

fit_compositional_pca_ilr_sc <- function(x_data,
                                         max_iter = 50,
                                         r = 10,
                                         lambda = 1,
                                         eps = 0.01,
                                         sc_factor = 0.001,
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
      # E-Step ###################
      for (i in seq_along(x_data)) {
        optim_result <- tryCatch({
          optim(rep(0, length = length(pca$sdev)),
                conditional_scores_log_ilr_sc,
                gr = gradient_cslc_ilr_sc,
                x_data_i = x_data[[i]],
                pca = pca,
                basis_matrix = basis_matrix,
                sc_factor = sc_factor,
                control = list(fnscale = -1),
                method = "BFGS")
        }, error = function(e) {
            cat("error optim() in iteration", k, "for observation", i, "\n")
            cat("error message:", e$message, "\n")
            print("pca$sdev:")
            print(pca$sdev)
            print("optim result:")
            print(optim_result)
            print("pca center:")
            print(pca$center)
        })
        scores_median <- as.vector(optim_result$par)
        scores_median_list[[i]] <- scores_median
        # importance sampling
        proposal_scores[[i]] <- sapply(1:(r * k), function(t) {
          matrix(rnorm(length(scores_median),
                       mean = scores_median,
                       sd = lambda * pca$sdev))
        })

        conditional_scores_list[[i]] <-
          apply(proposal_scores[[i]], 2, function(scores) {
            conditional_scores_log_ilr_sc(
                                          scores,
                                          x_data[[i]],
                                          pca,
                                          basis_matrix,
                                          sc_factor)
          })

        log_weights <- apply(proposal_scores[[i]], 2, function(scores) {
          conditional_scores_log_ilr_sc(
                                        scores,
                                        x_data[[i]],
                                        pca,
                                        basis_matrix,
                                        sc_factor) -
            sum(dnorm(
                      scores,
                      mean = scores_median,
                      sd = lambda * pca$sdev,
                      log = TRUE))
        })
        # increase numerical stability
        if (sum_exp) {
          weights[[i]] <- stabilize_weights(log_weights)
        } else {
          log_weights <- log_weights - mean(log_weights, na.rm = TRUE)
          weights[[i]] <- exp(log_weights)/sum(exp(log_weights))
        }
        # End of E-Step ###################
      }
      monitor_global_ess(weights, k)
      mean_conditional <- mean(unlist(conditional_scores_list), na.rm = TRUE)

      cat(sprintf("Conditional score mean value %.2f:\n",
                  mean_conditional))

      # M-Step ###################
      mu_scores <- rowMeans(sapply(seq_along(weights), function(i){
        proposal_scores[[i]] %*% weights[[i]]
      }), na.rm = TRUE)
      cat("Mean scores:", mu_scores, "\n")
      # update parameters
      pca_old <- pca
      pca$center <- pca$rotation %*% mu_scores
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
      cat("Eigenvalues:", pca$sdev, "\n")
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
                    "time" = elapsed_time))
      }
    }
  }
  constant <- apply(pca$rotation, 2, function(g) {
    sqrt(sum(g^2))
  })
  pca$rotation <- t(t(pca$rotation)/constant)
  pca$sdev <- pca$sdev * constant
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  return(list("iteration" = max_iter,
              "pca" = pca,
              "x_data" = x_data,
              "list_center" = center_list,
              "list_sdev" = sdev_list,
              "time" = elapsed_time))
}

conditional_scores_log_ilr_sc <- function(scores,
                                          x_data_i,
                                          pca,
                                          basis_matrix,
                                          sc_factor) {
  scaling_factor <- sc_factor
  ilr_comp <- as.vector(pca$center + pca$rotation %*% scores)
  clr_comp <- ilr2clr(ilr_comp)
  norm_constant <- sum(exp(clr_comp))

  # Compute scaled log likelihood
  log_likelihood <- scaling_factor * (sum(x_data_i * clr_comp)
                                      - sum(x_data_i) * log(norm_constant))

  # Prior remains unchanged as it's already well-scaled
  log_prior <- - sum(0.5 * scores^2 / (pca$sdev^2))

  return(log_likelihood + log_prior)
}

gradient_cslc_ilr_sc <- function(scores,
                                 x_data_i,
                                 pca,
                                 basis_matrix,
                                 sc_factor) {
  scaling_factor <- sc_factor
  m_i <- sum(x_data_i)
  ilr_comp <- as.vector(pca$center + pca$rotation %*% scores)
  clr_comp <- ilr2clr(ilr_comp)
  composition <- clrInv_long(clr_comp)

  grad <- sapply(seq_along(scores), function(k) {
    e_k <- basis_matrix[k, ] * (-1)
    term1 <- sum(x_data_i * e_k)
    term2 <- m_i * sum(composition * e_k)

    grad_k <- scaling_factor * (term1 - term2) - scores[k] / (pca$sdev[k]^2)

    return(grad_k)
  })

  return(grad)
}

clrInv_long <- function(clr_coords) {
    # Exponentiate the clr coordinates
    exp_coords <- exp(clr_coords)
    
    # Calculate the geometric mean normalization constant
    norm_const <- sum(exp_coords)
    
    # Return normalized compositions
    exp_coords / norm_const
}
