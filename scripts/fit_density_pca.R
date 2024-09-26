fit_density_pca <- function(x_data, x_grid = seq(min(unlist(x_data)), max(unlist(x_data)), length = 200),
                            max_iter = 50, r = 10, lambda = 1, dim_reduction = 0.001,
                            bw = (max(x_grid) - min(x_grid))/10, eps = 0.01){
  # initial estimates
  # kernel density estimates
  densities_estimated <- lapply(1:length(x_data), function(i){
    density <- density(x_data[[i]], from = min(x_grid), to = max(x_grid), 
                       kernel = "gaussian", bw, 
                       n = length(x_grid))
    data.frame("x" = density$x, "y" = density$y)
  })
  # compute initial pca
  clr_densities_estimated <- lapply(densities_estimated, clr_trafo)
  clr_densities <- do.call("rbind", sapply(clr_densities_estimated, '[', 2))
  pca <- prcomp(na.omit(clr_densities))
  which_reduced <- rev(cumsum(rev(pca$sdev^2))/sum(pca$sdev^2) > dim_reduction)
  which_reduced <- which_reduced|c(TRUE, TRUE, rep(FALSE, length(which_reduced) - 2))
  pca$sdev <- pca$sdev[which_reduced]
  pca$rotation <- pca$rotation[,which_reduced, drop = FALSE]
  
  proposal_scores <- list(length(x_data))
  weights <- list(length(x_data))
  if(max_iter > 0){
    for(k in 1:max_iter){
      # E-Step ###################
      # draw densities conditional on observations and current pca
      for(i in 1:nrow(clr_densities)){
        # find median of the posterior score distribution
        optim_result <- optim(rep(0, length = length(pca$sdev)), conditional_scores_log_density, gr = gradient_csld,
                              x_grid = x_grid, x_data_i = x_data[[i]], pca = pca,
                              control = list(fnscale = -1), method = "BFGS")
        scores_median <- as.vector(optim_result$par)
        # importance sampling
        proposal_scores[[i]] <- sapply(1:(r*k), function(t){
          matrix(rnorm(length(scores_median), mean = scores_median, sd = lambda*pca$sdev))
        })
        log_weights <- apply(proposal_scores[[i]], 2, function(scores){
          conditional_scores_log_density(scores, x_grid, x_data[[i]], pca) -
            sum(dnorm(scores, mean = scores_median, sd = lambda*pca$sdev, log = TRUE))
        })
        # increase numerical stability
        log_weights <- log_weights - mean(log_weights, na.rm = TRUE)
        weights[[i]] <- exp(log_weights)/sum(exp(log_weights))
      }
      # M-Step ###################
      mu_scores <- rowMeans(sapply(seq_along(weights), function(i){
        proposal_scores[[i]]%*%weights[[i]]
      }))
      
      # update pca
      pca_old <- pca
      pca$center <- center_function(cbind(x_grid, pca$center + pca$rotation%*%mu_scores))[,2]
      
      Sigma <- Reduce("+", lapply(seq_along(weights), function(i){
        Reduce("+", lapply(1:(r*k), function(t){
          C_it <- weights[[i]][t]*(proposal_scores[[i]][,t] - mu_scores)%*%
            t((proposal_scores[[i]][,t] - mu_scores))
        }))
      }))/length(weights)
      eigen_decomp <- eigen(Sigma)
      pca$sdev <- sqrt(eigen_decomp$values)
      pca$rotation <- pca$rotation%*%eigen_decomp$vectors
      pca$rotation <- apply(pca$rotation, 2, function(g) center_function(cbind(x_grid, g))[,2])
      
      # check convergence
      critical_value_1 <- L_2_norm(cbind(x_grid, pca_old$center - pca$center))
      K_old <- Reduce("+", lapply(seq_along(pca_old$sdev), function(k){
        pca_old$rotation[,k]%*%t(pca_old$rotation[,k])*(pca_old$sdev[k]^2)
      }))
      K_new <- Reduce("+", lapply(seq_along(pca$sdev), function(k){
        pca$rotation[,k]%*%t(pca$rotation[,k])*(pca$sdev[k]^2)
      }))
      critical_value_2 <- L_2_norm(cbind(x_grid, sapply(1:nrow(K_old), function(k){
        L_2_norm(cbind(x_grid, K_old[k,] - K_new[k,]))
      })))
      
      if(max(critical_value_1, critical_value_2) < eps){
        # normalize result
        constant <- apply(pca$rotation, 2,  function(g){
          L_2_norm(cbind(x_grid, g))
        })
        pca$rotation <- t(t(pca$rotation)/constant)
        pca$sdev <- pca$sdev*constant
        return(list("iteration" = k, "pca" = pca, "x_grid" = x_grid, "x_data" = x_data))
      }
      which_reduced <- rev(cumsum(rev(pca$sdev^2))/sum(pca$sdev^2) > dim_reduction)
      which_reduced <- which_reduced|c(TRUE, TRUE, rep(FALSE, length(which_reduced) - 2))
      pca$sdev <- pca$sdev[which_reduced]
      pca$rotation <- pca$rotation[,which_reduced, drop = FALSE]
    }
  }
  # normalize result
  constant <- apply(pca$rotation, 2,  function(g){
    L_2_norm(cbind(x_grid, g))/(max(x_grid) - min(x_grid))
  })
  pca$rotation <- t(t(pca$rotation)/constant)
  pca$sdev <- pca$sdev*constant
  return(list("iteration" = max_iter, "pca" = pca, "x_grid" = x_grid, "x_data" = x_data))
}

################################################################################
# objective function and gradient
conditional_scores_log_density <- function(scores, x_grid, x_data_i, pca){
  clr_density <- cbind(x_grid, pca$center + pca$rotation%*%scores)
  idxs <- sapply(x_data_i, function(x) {
    which.min((x - x_grid) ^ 2)
  })
  mid_points <- c(x_grid[1], x_grid[-1] - 0.5*diff(x_grid), x_grid[length(x_grid)])
  f_integral <- sum(exp(clr_density[,2])*diff(mid_points))
  
  sum(clr_density[idxs, 2]) - length(idxs)*log(f_integral) - sum(0.5*scores^2/(pca$sdev^2))
}

gradient_csld <- function(scores, x_grid, x_data_i, pca){
  idxs <- sapply(x_data_i, function(x) {
    which.min((x - x_grid) ^ 2)
  })
  mid_points <- c(x_grid[1], x_grid[-1] - 0.5*diff(x_grid), x_grid[length(x_grid)])
  density <- inverse_clr_trafo(cbind(x_grid, pca$center + pca$rotation%*%scores))
  
  sapply(seq_along(scores), function(k){
    scalar_prod <- sum(density[,2]*pca$rotation[, k]*diff(mid_points))
    sum(pca$rotation[idxs, k]) - length(idxs)*scalar_prod - scores[k]/(pca$sdev[k]^2)
  })
}

################################################################################
# helper functions
center_function <- function(f_data){
  mid_points <- c(f_data[1,1], f_data[-1,1] - 0.5*diff(f_data[,1]), f_data[nrow(f_data),1])
  f_integral <- sum(f_data[,2]*diff(mid_points))
  f_data[,2] <- f_data[,2] - f_integral/(mid_points[length(mid_points)] - mid_points[1])
  f_data
}
L_2_norm <- function(f_data){
  mid_points <- c(f_data[1,1], f_data[-1,1] - 0.5*diff(f_data[,1]), f_data[nrow(f_data),1])
  sqrt(sum(f_data[,2]^2*diff(mid_points)))
}

clr_trafo <- function(f_data){
  f_data[,2] <- log(f_data[,2])
  center_function(f_data)
}

inverse_clr_trafo <- function(clr_density){
  mid_points <- c(clr_density[1, 1], clr_density[-1, 1] - 0.5*diff(clr_density[,1]),
                  clr_density[nrow(clr_density),1])
  f_integral <- sum(exp(clr_density[,2])*diff(mid_points))
  data.frame("x" = clr_density[,1], "y" = exp(clr_density[,2])/f_integral)
}

predict_latent_densities <- function(density_pca){
  predicted_scores <- sapply(seq_along(density_pca$x_data), function(i){
    optim_result <- optim(rep(0, length = length(density_pca$pca$sdev)), conditional_scores_log_density, gr = gradient_csld,
                          x_grid = density_pca$x_grid, x_data_i = density_pca$x_data[[i]], pca = density_pca$pca,
                          control = list(fnscale = -1), method = "BFGS")
    as.vector(optim_result$par)
  })
  clr_densities <- density_pca$pca$center + density_pca$pca$rotation%*%predicted_scores
  return(list("clr_densities" = clr_densities, "x_grid" = density_pca$x_grid, "predicted_scores" = predicted_scores))
}

################################################################################
# plot functions

# plot all functional principal components
plot_pca <- function(pca, x_grid, dim = 2){
  pca_mean <- data.frame("x" = x_grid, "y" = pca$center)
  pcs <- lapply(1:dim, function(i){
    data.frame("x" = x_grid, "y" = pca$rotation[,i])
  })
  y_lim <- range(rbind(pca_mean, do.call("rbind", pcs))[,2])
  plot_function <- function() {
    plot(pca_mean, ylim = y_lim, type = "l", main = "Principal component decomposition")
    invisible(lapply(1:dim, function(i) lines(pcs[[i]], col = rainbow(dim)[i])))
    legend("bottom", legend = c("mean", paste("pc", 1:dim)), col = c("black", rainbow(dim)),
          lty = 1)
  }

  return(plot_function)
}

# obtain data to plot effect of one principal component on mean density
get_predicted_densities <- function(density_pca, idx = 1, fac = 0.2){
  clr_mean <- density_pca$pca$center
  clr_mean_minus <- clr_mean - fac*density_pca$pca$rotation[,idx]
  clr_mean_plus <- clr_mean + fac*density_pca$pca$rotation[,idx]
  
}

################################################################################
# Simulation Function
simulate_densities_greven <- function(n_data, n_samples, x_grid, lambda_1, lambda_2){
set.seed(12)
f_data <- data.frame("x" = x_grid, "y" = -20*(x_grid-0.5)^2+5/3)
clr_mean <- center_function(f_data)

pc_1 <- data.frame("x" = x_grid, "y" = -0.2*sin(10*(x_grid-0.5)))
pc_1 <- center_function(pc_1)
pc_1[,2] <- pc_1[,2]/L_2_norm(pc_1)

pc_2 <- data.frame("x" = x_grid, "y" = 0.1*cos(2*pi*(x_grid-0.5)))
pc_2 <- center_function(pc_2)
pc_2[,2] <- pc_2[,2]/L_2_norm(pc_2)

true_observed_clr_densities <- sapply(1:n_data, function(i){
  clr_mean[,2] + rnorm(1, 0, lambda_1)*pc_1[,2] + rnorm(1, 0, lambda_2)*pc_2[,2]
})
true_observed_densities <- lapply(1:n_data, function(i){
  clr_density <- data.frame(x_grid, true_observed_clr_densities[,i])
  inverse_clr_trafo(clr_density)
})

x_data_densities <- lapply(1:n_data, function(i){
  probs <- true_observed_densities[[i]][,2]
  x_grid <- true_observed_densities[[i]][,1]
  sample(x_grid, n_samples, replace = TRUE, prob = probs)
})

return(list("x_data_densities" = x_data_densities, "true_observed_densities" = true_observed_densities, "clr_mean" = clr_mean, "pc_1" = pc_1, "pc_2" = pc_2, "x_grid" = x_grid))
}


