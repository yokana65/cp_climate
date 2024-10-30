fit_compositional_pca_vs1_0 <- function(x_data,
                            max_iter = 50, r = 10, lambda = 1, dim_reduction = 0.001,
                            eps = 0.01){
  start_time <- Sys.time()
  # TODO: error checks for structure of x_data
  # initial estimates
  D <- length(x_data[[1]])
  nu <- rep(0, D)
  Sigma <- diag(D)
  # compute initial pca
  pca <- prcomp(Sigma)
  pca$rotation <- apply(pca$rotation, 2, function(g) g - mean(g))
  pca$center <- nu
  
  proposal_scores <- list(length(x_data))
  weights <- list(length(x_data))
  if(max_iter > 0){
    for(k in 1:max_iter){
      # E-Step ###################
      for(i in 1:length(x_data)){
        # error check to analyse `vmmin` is not finite 
        optim_result <- optim(rep(0, length = length(pca$sdev)), conditional_scores_log_clr_composition, gr = gradient_cslc,
                              x_data_i = x_data[[i]], pca = pca,
                              control = list(fnscale = -1), method = "BFGS")
        scores_median <- as.vector(optim_result$par)
        # importance sampling
        proposal_scores[[i]] <- sapply(1:(r*k), function(t){
          matrix(rnorm(length(scores_median), mean = scores_median, sd = lambda*pca$sdev))
        })
        log_weights <- apply(proposal_scores[[i]], 2, function(scores){
          conditional_scores_log_clr_composition(scores, x_data[[i]], pca) -
            sum(dnorm(scores, mean = scores_median, sd = lambda*pca$sdev, log = TRUE))
        })
        # increase numerical stability
        log_weights <- log_weights - mean(log_weights, na.rm = TRUE)
        weights[[i]] <- exp(log_weights)/sum(exp(log_weights))
      }
      if (any(!is.finite(weights[[i]]))) {
         stop(paste("Infinite or NaN values found in weights at iteration", k, "for observation", i))
      }
      # M-Step ###################
      mu_scores <- rowMeans(sapply(seq_along(weights), function(i){
        proposal_scores[[i]]%*%weights[[i]]
      }))
      # update parameters
      pca_old <- pca
      pca$center <- pca$rotation%*%mu_scores
      Sigma <- Reduce("+", lapply(seq_along(weights), function(i){
        Reduce("+", lapply(1:(r*k), function(t){
          C_it <- weights[[i]][t]*(proposal_scores[[i]][,t] - mu_scores)%*%
            t((proposal_scores[[i]][,t] - mu_scores))
        }))
      }))/length(weights)

      eigen_decomp <- eigen(Sigma)
      # error check eigenvalues > 0
      negative_eigenvalues <- eigen_decomp$values < 0
        if (any(negative_eigenvalues)) {
        warning(sprintf("Warning: %d eigenvalues are negative. They have been set to zero.", sum(negative_eigenvalues)))
        }
      pca$sdev <- sqrt(pmax(eigen_decomp$values, 0))
      pca$rotation <- pca$rotation%*%eigen_decomp$vectors
      pca$rotation <- apply(pca$rotation, 2, function(g) g - mean(g))
      
      # check convergence
      critical_value_1 <- sqrt(sum((pca_old$center - pca$center)^2))
      K_old <- Reduce("+", lapply(seq_along(pca_old$sdev), function(k){
        # standard decomposition VLV
        pca_old$rotation[,k]%*%t(pca_old$rotation[,k])*(pca_old$sdev[k]^2)
      }))
      K_new <- Reduce("+", lapply(seq_along(pca$sdev), function(k){
        pca$rotation[,k]%*%t(pca$rotation[,k])*(pca$sdev[k]^2)
      }))
      K_diff <- K_old - K_new
      # Berechne die Frobenius-Norm der Differenzmatrix
      critical_value_2 <- norm(K_diff, type = "F")

      if(max(critical_value_1, critical_value_2) < eps){
        constant <- apply(pca$rotation, 2, function(g) { sqrt(sum(g^2)) })
        pca$rotation <- t(t(pca$rotation)/constant)
        pca$sdev <- pca$sdev*constant
        
        end_time <- Sys.time()
        elapsed_time <- end_time - start_time
        print(paste("The algorithm converged after:", elapsed_time, "seconds"))
        return(list("iteration" = k, "pca" = pca, "x_data" = x_data))
      }
    }
  }
  constant <- apply(pca$rotation, 2, function(g) { sqrt(sum(g^2)) })
  pca$rotation <- t(t(pca$rotation)/constant)
  pca$sdev <- pca$sdev*constant
  return(list("iteration" = max_iter, "pca" = pca, "x_data" = x_data))
}

conditional_scores_log_clr_composition <- function(scores, x_data_i, pca){
clr_comp <- pca$center + pca$rotation%*%scores
norm_constant <- sum(exp(clr_comp))
log_likelihood <- sum(x_data_i * clr_comp) - sum(x_data_i)*log(norm_constant)
log_prior <- - sum(0.5*scores^2/(pca$sdev^2))
log_posterior <- log_likelihood + log_prior

return(log_posterior) 
}

gradient_cslc <- function(scores, x_data_i, pca){

  m_i <- sum(x_data_i)

  composition <- clrInv(pca$center + pca$rotation%*%scores)

  grad <- sapply(seq_along(scores), function(k) {
    e_k <- pca$rotation[, k]

    term1 <- sum(x_data_i * e_k)

    term2 <- m_i * sum(composition * e_k)

    grad_k <- term1 - term2 - scores[k] / (pca$sdev[k]^2)

    return(grad_k)
})

return(grad)
}


############ Simulate composition data ############ 
simulate_composition_1 <- function(n_components, n_data, n_counts, n_samples, lambda_1, lambda_2){
set.seed(123)
x_grid <- seq(1, n_components)
raw_data <- data.frame("x" = x_grid, "y" = rep(0, n_components))

clr_mean <- center_function_comp(raw_data)

pc_1 <- data.frame("x" = x_grid, "y" = c(rep(sqrt(5/(5+1)*5^(-1)), 5),-1, rep(0, 7)))
pc_1 <- center_function_comp(pc_1)
pc_1[, 2] <- pc_1[, 2] / norm(pc_1[, 2], type = "2")

pc_2 <- data.frame("x" = x_grid, "y" = c(rep(sqrt(12/(12+1)*12^(-1)), 12),0))
pc_2 <- center_function_comp(pc_2)
pc_2[, 2] <- pc_2[,2]/norm(pc_2[,2], type="2")

true_observed_clr_comp <- sapply(1:n_data, function(i){
  clr_mean[,2] + rnorm(1, 0, lambda_1)*pc_1[,2] + rnorm(1, 0, lambda_2)*pc_2[,2]
})

check_columns_sum_to(true_observed_clr_comp, 0.001)

true_observed_comp <- lapply(1:n_data, function(i){
  # TODO: remove x_grid
  clr_density <- data.frame(x_grid, true_observed_clr_comp[,i])
  inverse_clr_trafo(clr_density)
})

x_data <- unlist(lapply(1:n_data, function(i) {
  probs <- true_observed_comp[[i]][,2]
  samples <- t(rmultinom(n_samples, n_counts, probs))
  lapply(1:nrow(samples), function(j) samples[j,])
}), recursive = FALSE)

x_data_matrix <- do.call(rbind, x_data)

return(list(x_data = x_data, x_data_matrix = x_data_matrix, true_observed_comp = true_observed_comp))
}

center_function_comp <- function(comp_data){
  mean <- mean(comp_data[,2])
  comp_data[,2] <- comp_data[,2] - mean
  comp_data
}

check_columns_sum_to <- function(data, integer) {
  n_cols <- ncol(data)
  
  columns_sum_to_zero <- logical(n_cols)
  
  for (i in 1:n_cols) {
    column_sum <- sum(data[, i])
    columns_sum_to_zero[i] <- (column_sum < integer)
  }
  
  return(columns_sum_to_zero)
}

inverse_clr_trafo <- function(clr_density){
  f_integral <- sum(exp(clr_density[,2]))
  data.frame("x" = clr_density[,1], "y" = exp(clr_density[,2])/f_integral)
}


