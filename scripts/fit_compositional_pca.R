co_pca_mcem <- function(x_data,
                        max_iter = 50,
                        r = 10,
                        lambda = 1,
                        eps = 0.01,
                        sum_exp = FALSE) {
  start_time <- Sys.time()

  prepared_data <- prepare_data(x_data)

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
                              log_conditional_scores,
                              gr = gradient_lcs,
                              x_data_i = prepared_data$x_data[[i]],
                              pca = pca,
                              basis_matrix = prepared_data$H,
                              control = list(fnscale = -1),
                              method = "BFGS")

        scores_median <- as.vector(optim_result$par)
        proposal_scores[[i]] <- sapply(1L:(r * k), function(t) {
          matrix(rnorm(length(scores_median),
                       mean = scores_median,
                       sd = lambda * pca$sdev))
        })

        log_weights <- apply(proposal_scores[[i]], 2L, function(scores) {
          log_conditional_scores(scores,
                                 prepared_data$x_data[[i]],
                                 pca,
                                 prepared_data$H) -
            sum(dnorm(scores,
                      mean = scores_median,
                      sd = lambda * pca$sdev,
                      log = TRUE))
        })
        if (sum_exp == TRUE) {
          weights[[i]] <- stabilize_weights(log_weights)
        } else {
          max_log_weight <- max(log_weights)
          weights_i <- exp(log_weights - max_log_weight)
          weights[[i]] <- weights_i / sum(weights_i)
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
        results <- finalize_pca(
            pca          = pca,
            prepared_data = prepared_data,
            x_data       = x_data,
            weights      = weights,
            start_time   = start_time
          )
        print(paste("The algorithm converged after:",
                    results$elapsed_time, attr(results$elapsed_time, "units")))
        return(list("iteration" = k,
                    "pca" = results$pca,
                    "x_data" = x_data,
                    "time" = results$elapsed_time,
                    "ESS" = results$ess))
      }
    }
  }
  results <- finalize_pca(
    pca          = pca,
    prepared_data = prepared_data,
    x_data       = x_data,
    weights      = weights,
    start_time   = start_time
  )
  print(paste("The algorithm did not converge after:",
              elapsed_time, 
              attr(elapsed_time, "units")))
  return(list("iteration" = k,
              "pca" = results$pca,
              "x_data" = x_data,
              "time" = results$elapsed_time,
              "ESS" = results$ess))
}

co_pca_mcem_nograd <- function(x_data,
                               max_iter = 50,
                               r = 10,
                               lambda = 1,
                               eps = 0.01,
                               sum_exp = FALSE) {
  start_time <- Sys.time()

  prepared_data <- prepare_data(x_data)

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
                              log_conditional_scores,
                              gr = gradient_lcs,
                              x_data_i = prepared_data$x_data[[i]],
                              pca = pca,
                              basis_matrix = prepared_data$H,
                              control = list(fnscale = -1),
                              method = "BFGS")

        scores_median <- as.vector(optim_result$par)
        proposal_scores[[i]] <- sapply(1L:(r * k), function(t) {
          matrix(rnorm(length(scores_median),
                       mean = scores_median,
                       sd = lambda * pca$sdev))
        })

        log_weights <- apply(proposal_scores[[i]], 2L, function(scores) {
          log_conditional_scores(scores,
                                 prepared_data$x_data[[i]],
                                 pca,
                                 prepared_data$H) -
            sum(dnorm(scores,
                      mean = scores_median,
                      sd = lambda * pca$sdev,
                      log = TRUE))
        })
        if (sum_exp == TRUE) {
          weights[[i]] <- stabilize_weights(log_weights)
        } else {
          max_log_weight <- max(log_weights)
          weights_i <- exp(log_weights - max_log_weight)
          weights[[i]] <- weights_i / sum(weights_i)
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
        results <- finalize_pca(
          pca          = pca,
          prepared_data = prepared_data,
          x_data       = x_data,
          weights      = weights,
          start_time   = start_time
        )
        print(paste("The algorithm converged after:",
                    results$elapsed_time,
                    attr(results$elapsed_time, "units")))
        return(list("iteration" = k,
                    "pca" = results$pca,
                    "x_data" = x_data,
                    "time" = results$elapsed_time,
                    "ESS" = results$ess))
      }
    }
  }
  results <- finalize_pca(
    pca          = pca,
    prepared_data = prepared_data,
    x_data       = x_data,
    weights      = weights,
    start_time   = start_time
  )
  print(paste("The algorithm did not converge after:",
              elapsed_time,
              attr(elapsed_time, "units")))

  return(list("iteration" = k,
              "pca" = results$pca,
              "x_data" = x_data,
              "time" = results$elapsed_time,
              "ESS" = results$ess))
}

