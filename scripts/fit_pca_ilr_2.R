# version mit ESS
fit_pca_ilr_vs_5 <- function(x_data,
                             max_iter = 50,
                             r = 10,
                             lambda = 1,
                             eps = 0.01,
                             sc_factor = 1,
                             sum_exp = TRUE) {
  start_time <- Sys.time()
  if (!is.list(x_data) && !is.matrix(x_data)) {
    stop("Input x_data must be a list or a matrix")
  }

  if (is.data.frame(x_data) || is.matrix(x_data)) {
    x_data <- apply(x_data, 1, function(x) x, simplify = FALSE)
  }

  lengths <- unique(sapply(x_data, length))
  if (length(lengths) != 1) {
    stop("All observations must have the same number of components")
  }
  D <- lengths

  basis_vectors <- lapply(1:(D - 1), generate_orthonormal_basis, D)
  basis_matrix <- do.call(rbind, basis_vectors)

  # initial estimates
  nu <- rep(0, D - 1)
  Sigma <- diag(D - 1)
  pca <- prcomp(Sigma, center = FALSE)
  pca$center <- nu

  proposal_scores <- list(length(x_data))
  weights <- list(length(x_data))
  sdev_list <- list(length(max_iter))
  center_list <- list(length(max_iter))
  conditional_scores_list <- list(length(x_data))
  scores_median_list <- list(length(x_data))


  if (max_iter > 0) {
    for (k in 1:max_iter) {
      cat("Iteration:", k, "\n")
      for(i in 1:length(x_data)){

        optim_result <- optim(rep(0, length = length(pca$sdev)),
                              conditional_scores_log_ilr,
                              gr = gradient_cslc_vs1,
                              x_data_i = x_data[[i]],
                              pca = pca,
                              basis_matrix = basis_matrix,
                              sc_factor = sc_factor,
                              control = list(fnscale = -1),
                              method = "BFGS")
    
        scores_median <- as.vector(optim_result$par)
        proposal_scores[[i]] <- sapply(1:(r * k), function(t) {
          matrix(rnorm(length(scores_median),
                       mean = scores_median,
                       sd = lambda * pca$sdev))
        })
    
        log_weights <- apply(proposal_scores[[i]], 2, function(scores) {
          conditional_scores_log_ilr(scores,
                                          x_data[[i]],
                                          pca,
                                          basis_matrix,
                                          sc_factor) -
            sum(dnorm(scores,
                      mean = scores_median,
                      sd = lambda * pca$sdev,
                      log = TRUE))
        })
        if (sum_exp == TRUE) {
          weights[[i]] <- stabilize_weights(log_weights)
        } else {
          max_log_weight <- max(log_weights)
          weights <- exp(log_weights - max_log_weight)
          normalized_weights <- weights / sum(weights)
        }
      }

      monitor_global_ess(weights, k)
      mean_ess <- mean(sapply(all_weights, function(w) 1/sum(w^2)))

      # M-Step ###################
      scores_matrix <- sapply(seq_along(weights), function(i){
        proposal_scores[[i]] %*% weights[[i]]
      })
      na_count <- sum(is.na(scores_matrix))  
      mu_scores <- rowMeans(scores_matrix, na.rm = TRUE)
      # update parameters
      pca_old <- pca
      pca$center <- pca$center + pca$rotation %*% mu_scores
      cat("center:", pca$center, "\n")
      center_list[[k]] <- pca$center
      Sigma <- Reduce("+", lapply(seq_along(weights), function(i) {
          Reduce("+", lapply(1:(r * k), function(t) {
          C_it <- weights[[i]][t] * (proposal_scores[[i]][, t] - mu_scores) %*%
              t((proposal_scores[[i]][, t] - mu_scores))
          }))
      })) / length(weights)
      eigen_decomp <-  tryCatch({eigen(Sigma)}, error = function(e) {
          cat("error eigen() in iteration", k, "for observation", i, "\n")
          cat("error message:", e$message, "\n")
          print("pca$sdev:")
          print(pca$sdev)
      })
      negative_eigenvalues <- eigen_decomp$values < 0
      if (any(negative_eigenvalues)) {
          warning(sprintf("Warning: %d eigenvalues are negative.\n
          They have been set to zero.",
                          sum(negative_eigenvalues)))
      }
      pca$sdev <- sqrt(pmax(eigen_decomp$values, 0))
      sdev_list[[k]] <- pca$sdev
      pca$rotation <- pca$rotation %*% eigen_decomp$vectors
      # check convergence
      critical_value_1 <- sqrt(sum((pca_old$center - pca$center)^2))
      cat("critical value center_diff:", critical_value_1, "\n")
      Sigma_old <- Reduce("+", lapply(seq_along(pca_old$sdev), function(k) {
          pca_old$rotation[, k] %*% t(pca_old$rotation[, k]) * (pca_old$sdev[k]^2)
      }))
      Sigma_new <- Reduce("+", lapply(seq_along(pca$sdev), function(k) {
        pca$rotation[, k] %*% t(pca$rotation[, k]) * (pca$sdev[k]^2)
      }))
      Sigma_diff <- Sigma_old - Sigma_new
      critical_value_2 <- norm(Sigma_diff, type = "F")
      cat("critical value Sigma_diff:", critical_value_2, "\n")
      clr_rotation <- t(basis_matrix) %*% pca$rotation %*% basis_matrix
      cat("clr-coordinates PCA1:", clr_rotation[ , 1], "\n")
      # pca$sdev <- pca$sdev / max(pca$sdev) # change log
      
      if (max(critical_value_1, critical_value_2) < eps) {
        constant <- apply(pca$rotation, 2, function(g) {
          sqrt(sum(g^2))
        })
        pca$rotation <- t(t(pca$rotation) / constant)
        pca$sdev <- pca$sdev * constant  
        end_time <- Sys.time()
        elapsed_time <- end_time - start_time
        print(paste("The algorithm converged after:", elapsed_time, "minutes"))
        return(list("iteration" = k,
                    "pca" = pca,
                    "x_data" = x_data,
                    "list_center" = center_list,
                    "list_sdev" = sdev_list,
                    "time" = elapsed_time,
                    "ESS" = mean_ess))
      }
    }
  }
  constant <- apply(pca$rotation, 2, function(g) {
    sqrt(sum(g^2))
  })
  pca$rotation <- t(t(pca$rotation) / constant) # unit 1 for each column
  clr_rotation <- t(basis_matrix) %*% pca$rotation %*% basis_matrix
  cat("clr-coordinates PCA1:", clr_rotation[ , 1], "\n")
  pca$sdev <- pca$sdev * constant
  cat("square root Eigenvalues:", pca$sdev, "\n")
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  return(list("iteration" = max_iter,
              "pca" = pca,
              "x_data" = x_data,
              "list_center" = center_list,
              "list_sdev" = sdev_list,
              "time" = elapsed_time,
              "ESS" = mean_ess))
}

