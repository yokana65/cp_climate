# version mit ESS




fit_pca_vs_5 <- function(x_data,
                         max_iter = 50,
                         r = 10,
                         lambda = 1,
                         eps = 0.01,
                         sc_factor = 1,
                         scores = TRUE,
                         sum_exp = TRUE,
                         fix_sign = TRUE) {
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
  V <- get_helmert(x_data)
  nu <- rep(0, D - 1)
  Sigma <- diag(D - 1)
  pca <- prcomp(Sigma, center = FALSE)
  pca$center <- nu

  proposal_scores <- list(length(x_data))
  weights <- list(length(x_data))

  if (max_iter > 0) {
    for (k in 1L:max_iter) {
      for (i in seq_along(x_data)){
        optim_result <- optim(rep(0L, length = length(pca$sdev)),
                              conditional_scores_log_ilr,
                              gr = gradient_cslc_vs1,
                              x_data_i = x_data[[i]],
                              pca = pca,
                              basis_matrix = t(V),
                              sc_factor = sc_factor,
                              control = list(fnscale = -1),
                              method = "BFGS")

        scores_median <- as.vector(optim_result$par)
        proposal_scores[[i]] <- sapply(1L:(r * k), function(t) {
          matrix(rnorm(length(scores_median),
                       mean = scores_median,
                       sd = lambda * pca$sdev))
        })
        log_weights <- apply(proposal_scores[[i]], 2L, function(scores) {
          conditional_scores_log_ilr(scores,
                                     x_data[[i]],
                                     pca,
                                     basis_matrix = t(V),
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
          weights[[i]] <- weights / sum(weights)
        }
      }
      # M-Step ###################
      scores_matrix <- sapply(seq_along(weights), function(i){
        proposal_scores[[i]] %*% weights[[i]]
      })
      mu_scores <- rowMeans(scores_matrix, na.rm = TRUE)
      # update parameters
      pca_old <- pca
      pca$center <- pca$center + pca$rotation %*% mu_scores
      Sigma <- Reduce("+", lapply(seq_along(weights), function(i) {
        Reduce("+", lapply(1L:(r * k), function(t) {
          C_it <- weights[[i]][t] * (proposal_scores[[i]][, t] - mu_scores) %*%
            t((proposal_scores[[i]][, t] - mu_scores))
        }))
      })) / length(weights)
      edc <-  tryCatch({eigen(Sigma, symmetric = TRUE)}, error = function(e) {
        cat("error eigen() in iteration", k, "for observation", i, "\n")
        cat("error message:", e$message, "\n")
      })
      edc$vectors <- fix_sign(edc$vectors)
      ev <- edc$values
      if (any(neg <- ev < 0)) {
        if (any(ev[neg] < -9 * .Machine$double.eps * ev[1L])) 
          stop("covariance matrix is not non-negative definite")
        else ev[neg] <- 0
      }
      pca$sdev <- sqrt(edc$values)
      sdev_list[[k]] <- pca$sdev
      # TODO: Das ist der entscheidende Teil für die Konvergenz: herausarbeiten! -> wird hier eigentlich auch eine Vorzeichenkorrektur benötigt?
      pca$rotation <- pca$rotation %*% edc$vectors
      critical_value_1 <- sqrt(sum((pca_old$center - pca$center)^2))
      Sigma_old <- Reduce("+", lapply(seq_along(pca_old$sdev), function(k) {
        pca_old$rotation[, k] %*% t(pca_old$rotation[, k]) * (pca_old$sdev[k]^2)
      }))
      Sigma_new <- Reduce("+", lapply(seq_along(pca$sdev), function(k) {
        pca$rotation[, k] %*% t(pca$rotation[, k]) * (pca$sdev[k]^2)
      }))
      Sigma_diff <- Sigma_old - Sigma_new
      critical_value_2 <- norm(Sigma_diff, type = "F")
      if (max(critical_value_1, critical_value_2) < eps) {
        constant <- apply(pca$rotation, 2, function(g) {
          sqrt(sum(g^2))
        })
        pca$rotation <- t(t(pca$rotation) / constant)
        pca$rotation_ilr <- pca$rotation
        pca$rotation <- V %*% pca$rotation_ilr
        pca$sdev <- pca$sdev / sum(pca$sdev)
        pca$scores_clr <- if (scores) {
          scale(clr(x_data),
                center = pca$center, scale = FALSE) %*% pca$rotation
        }
        rownames(pca$rotation) <- names(x_data[[1L]])
        cn <- paste0("Comp.", seq_len(ncol(pca$rotation)))
        colnames(pca$rotation) <- cn
        mean_ess <- mean(sapply(weights, function(w) 1 / sum(w^2)))
        end_time <- Sys.time()
        elapsed_time <- end_time - start_time
        print(paste("The algorithm converged after:", elapsed_time, attr(elapsed_time, "units")))
        return(list("iteration" = k,
                    "pca" = pca,
                    "x_data" = x_data,
                    "time" = elapsed_time,
                    "ESS" = mean_ess))
      }
    }
  }
  constant <- apply(pca$rotation, 2, function(g) {
    sqrt(sum(g^2))
  })
  pca$rotation <- t(t(pca$rotation) / constant)
  pca$rotation_ilr <- pca$rotation
  pca$rotation <- V %*% pca$rotation_ilr
  pca$sdev <- pca$sdev / sum(pca$sdev)
  pca$scores_clr <- if (scores) {
    scale(clr(x_data),
          center = pca$center, scale = FALSE) %*% pca$rotation
  }
  rownames(pca$rotation) <- names(x_data[[1L]])
  cn <- paste0("Comp.", seq_len(ncol(pca$rotation)))
  colnames(pca$rotation) <- cn
  mean_ess <- mean(sapply(weights, function(w) 1 / sum(w^2)))
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  return(list("iteration" = max_iter,
              "pca" = pca,
              "x_data" = x_data,
              "time" = elapsed_time,
              "ESS" = mean_ess))
}

