CoDaPcaMcEm <- function(x_data,
                        max_iter = 50,
                        r = 10,
                        lambda = 1,
                        eps = 0.01,
                        scores = FALSE,
                        sum_exp = TRUE,
                        fix_sign = TRUE,
                        dim = NULL,
                        coordinates = TRUE) {
  start_time <- Sys.time()
  
  prepared_data <- prepare_data(x_data, if(dim != FALSE) dim else NULL)
  
  nu <- rep(0, prepared_data$D - 1)
  Sigma <- diag(prepared_data$D - 1)
  pca <- prcomp(Sigma, center = FALSE)
  pca$center <- nu

  proposal_scores <- vector("list", length(prepared_data$x_data))
  weights <- vector("list", length(prepared_data$x_data))

  if (max_iter > 0) {
    for (k in 1L:max_iter) {
      for (i in seq_along(prepared_data$x_data)){
        optim_result <- optim(rep(0L, length = length(pca$sdev)),
                              conditional_scores_log_ilr_vs3b,
                              gr = gradient_cslc_vs1c,
                              x_data_i = prepared_data$x_data[[i]],
                              pca = pca,
                              basis_matrix = t(prepared_data$H),
                              control = list(fnscale = -1),
                              method = "BFGS")
    
        scores_median <- as.vector(optim_result$par)
        proposal_scores[[i]] <- sapply(1L:(r * k), function(t) {
          matrix(rnorm(length(scores_median),
                       mean = scores_median,
                       sd = lambda * pca$sdev))
        })

        log_weights <- apply(proposal_scores[[i]], 2L, function(scores) {
          conditional_scores_log_ilr_vs3b(scores,
                                          prepared_data$x_data[[i]],
                                          pca,
                                          t(prepared_data$H)) -
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

      pca_old <- pca
      pca$center <- as.vector(pca$center + pca$rotation %*% mu_scores)
      Sigma <- Reduce("+", lapply(seq_along(weights), function(i) {
          Reduce("+", lapply(1:(r * k), function(t) {
          C_it <- weights[[i]][t] * (proposal_scores[[i]][, t] - mu_scores) %*%
              t((proposal_scores[[i]][, t] - mu_scores))
          }))
      })) / length(weights)
      eigen_decomp <-  eigen(Sigma)
      pca$sdev <- sqrt(pmax(eigen_decomp$values, 0))
      pca$rotation <- pca$rotation %*% eigen_decomp$vectors

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
        pca$rotation <- prepared_data$H %*% pca$rotation_ilr
        pca$center <- prepared_data$H %*% pca$center
        pca$sdev <- pca$sdev / sum(pca$sdev)
        pca$scores_clr <- if (scores) {
          scale(clr(prepared_data$x_df),
                center = pca$center, scale = FALSE) %*% pca$rotation
        }
        # pca$coordinates <- if (coordinates) {
        #   predict_coordinates(comp_pca = pca, x_df = prepared_data$x_df, basis_matrix = t(H))
        # }
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
  pca$rotation <- t(t(pca$rotation) / constant) # unit 1 fom each column
  pca$sdev <- pca$sdev * constant
  pca$sdev <- pca$sdev / max(pca$sdev) # change log
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  return(list("iteration" = max_iter,
              "pca" = pca,
              "x_data" = x_data,
              "time" = elapsed_time))
}

CoDaPcaMcEm_num <- function(x_data,
                        max_iter = 50,
                        r = 10,
                        lambda = 1,
                        eps = 0.01,
                        scores = FALSE,
                        sum_exp = TRUE,
                        fix_sign = TRUE,
                        dim = NULL) {
  start_time <- Sys.time()
  
  prepared_data <- prepare_data(x_data, if(dim != FALSE) dim else NULL)
  
  nu <- rep(0, prepared_data$D - 1)
  Sigma <- diag(prepared_data$D - 1)
  pca <- prcomp(Sigma, center = FALSE)
  pca$center <- nu

  proposal_scores <- list(length(prepared_data$x_data))
  weights <- list(length(prepared_data$x_data))

  cat(sprintf("MCEM starts"))

  if (max_iter > 0) {
    for (k in 1L:max_iter) {
      cat(sprintf("Iteration %d", k))
      for (i in seq_along(prepared_data$x_data)){
        optim_result <- optim(rep(0L, length = length(pca$sdev)),
                              conditional_scores_log_ilr_vs3b,
                              x_data_i = prepared_data$x_data[[i]],
                              pca = pca,
                              basis_matrix = t(prepared_data$H),
                              control = list(fnscale = -1),
                              method = "Nelder-Mead")
    
        scores_median <- as.vector(optim_result$par)
        proposal_scores[[i]] <- sapply(1L:(r * k), function(t) {
          matrix(rnorm(length(scores_median),
                       mean = scores_median,
                       sd = lambda * pca$sdev))
        })

        log_weights <- apply(proposal_scores[[i]], 2L, function(scores) {
          conditional_scores_log_ilr_vs3b(scores,
                                          prepared_data$x_data[[i]],
                                          pca,
                                          t(prepared_data$H)) -
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
      # M-Step ###################
      scores_matrix <- sapply(seq_along(weights), function(i){
        proposal_scores[[i]] %*% weights[[i]]
      })
      mu_scores <- rowMeans(scores_matrix, na.rm = TRUE)
      # update parameters
      pca_old <- pca
      pca$center <- as.vector(pca$center + pca$rotation %*% mu_scores)
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
      pca$rotation <- pca$rotation %*% eigen_decomp$vectors
      # check convergence
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
        pca$rotation <- prepared_data$H %*% pca$rotation_ilr
        pca$center <- prepared_data$H %*% pca$center
        pca$sdev <- pca$sdev / sum(pca$sdev)
        pca$scores_clr <- if (scores) {
          scale(clr(x_df),
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
  pca$rotation <- t(t(pca$rotation) / constant) # unit 1 fom each column
  pca$sdev <- pca$sdev * constant
  pca$sdev <- pca$sdev / max(pca$sdev) # change log
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  return(list("iteration" = max_iter,
              "pca" = pca,
              "x_data" = x_data,
              "time" = elapsed_time))
}

CoDaPcaMcEm_princomp <- function(x_data,
                        max_iter = 50,
                        r = 10,
                        lambda = 1,
                        eps = 0.01,
                        scores = FALSE,
                        sum_exp = TRUE,
                        fix_sign = TRUE,
                        dim = NULL) {
  start_time <- Sys.time()
  
  prepared_data <- prepare_data(x_data, if(dim != FALSE) dim else NULL)
  
  nu <- rep(0, prepared_data$D - 1)
  Sigma <- diag(prepared_data$D - 1)
  pca <- princomp(cov = Sigma)
  pca$center <- nu

  proposal_scores <- list(length(prepared_data$x_data))
  weights <- list(length(prepared_data$x_data))

  if (max_iter > 0) {
    for (k in 1L:max_iter) {
      cat(sprintf("Iteration %d", k))
      for (i in seq_along(prepared_data$x_data)){
        optim_result <- optim(rep(0L, length = length(pca$sdev)),
                              conditional_scores_log_ilr_vs4,
                              gr = gradient_cslc_vs4,
                              x_data_i = prepared_data$x_data[[i]],
                              pca = pca,
                              basis_matrix = t(prepared_data$H),
                              control = list(fnscale = -1),
                              method = "BFGS")
    
        scores_median <- as.vector(optim_result$par)
        proposal_scores[[i]] <- sapply(1L:(r * k), function(t) {
          matrix(rnorm(length(scores_median),
                       mean = scores_median,
                       sd = lambda * pca$sdev))
        })
    
        log_weights <- apply(proposal_scores[[i]], 2L, function(scores) {
          conditional_scores_log_ilr_vs3(scores,
                                          prepared_data$x_data[[i]],
                                          pca,
                                          t(prepared_data$H)) -
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

      # M-Step ###################
      scores_matrix <- sapply(seq_along(weights), function(i){
        proposal_scores[[i]] %*% weights[[i]]
      })
      na_count <- sum(is.na(scores_matrix))  
      mu_scores <- rowMeans(scores_matrix, na.rm = TRUE)
      # update parameters
      pca_old <- pca
      pca$center <- pca$center + pca$loading %*% mu_scores
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
      pca$loading <- pca$loading %*% eigen_decomp$vectors
      # check convergence
      critical_value_1 <- sqrt(sum((pca_old$center - pca$center)^2))
      Sigma_old <- Reduce("+", lapply(seq_along(pca_old$sdev), function(k) {
          pca_old$loading[, k] %*% t(pca_old$loading[, k]) * (pca_old$sdev[k]^2)
      }))
      Sigma_new <- Reduce("+", lapply(seq_along(pca$sdev), function(k) {
        pca$loading[, k] %*% t(pca$loading[, k]) * (pca$sdev[k]^2)
      }))
      Sigma_diff <- Sigma_old - Sigma_new
      critical_value_2 <- norm(Sigma_diff, type = "F")
      
      if (max(critical_value_1, critical_value_2) < eps) {
        constant <- apply(pca$loading, 2, function(g) {
          sqrt(sum(g^2))
        })
        pca$loading <- t(t(pca$loading) / constant)
        pca$loading_ilr <- pca$loading
        pca$loading <- prepared_data$H %*% pca$loading_ilr
        pca$center <- prepared_data$H %*% pca$center
        pca$sdev <- pca$sdev / sum(pca$sdev)
        pca$scores_clr <- if (scores) {
          scale(clr(x_df),
                center = pca$center, scale = FALSE) %*% pca$loading
        }
        rownames(pca$loading) <- names(x_data[[1L]])
        cn <- paste0("Comp.", seq_len(ncol(pca$loading)))
        colnames(pca$loading) <- cn
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
  constant <- apply(pca$loading, 2, function(g) {
    sqrt(sum(g^2))
  })
  pca$loading <- t(t(pca$loading) / constant) # unit 1 fom each column
  pca$sdev <- pca$sdev * constant
  pca$sdev <- pca$sdev / max(pca$sdev) # change log
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  return(list("iteration" = max_iter,
              "pca" = pca,
              "x_data" = x_data,
              "time" = elapsed_time))
}


validate_dimensions <- function(x_data, dim = NULL) {
    lengths <- unique(sapply(x_data, length))
    if (length(lengths) != 1) {
        stop("All observations must have the same number of components")
    }
    
    D <- lengths[1]
    
    if (!is.null(dim)) {
        if (!is.numeric(dim) || dim != as.integer(dim)) {
            stop("dim must be an integer")
        }
        if (dim > D) {
            stop(sprintf("dim must be less than or equal to the original dimension %d", D))
        }
        return(dim)
    }
    
    return(D)
}

prepare_data <- function(x_data, dim = NULL) {
    # Input validation
    if (!is.list(x_data) && !is.matrix(x_data)) {
        stop("Input x_data must be a list or a matrix")
    }
    
    # Convert to consistent format
    if (is.list(x_data)) {
        x_df <- do.call(rbind, x_data)
    } else if (is.data.frame(x_data) || is.matrix(x_data)) {
        x_df <- x_data
        x_data <- apply(x_data, 1, function(x) x, simplify = FALSE)
    }
    
    # Validate dimensions
    D <- validate_dimensions(x_data, dim)
    
    # Create Helmert matrix
    H <- get_helmert(if(!is.null(dim)) dim else D)
    
    # Return as named list
    list(
        x_df = x_df,
        x_data = x_data,
        H = H,
        D = D
    )
}
