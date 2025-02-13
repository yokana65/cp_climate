set_1 <- c("#e0ecf4", "#9ebcda", "#8856a7")

normalize <- function(v) {
  v / sqrt(sum(v^2))
}

base_setting <- function(n_observations = 100,
                            eigenvalues = c(0.7, 0.3, 0.0, 0.0),
                            mean = c(0, 0.01, 0.02, -0.02, -0.01),
                            n_counts = 500) {
pc_1 <- c(0, sqrt(1/2), -sqrt(1/2), 0, 0)
pc_2 <- c(-sqrt(1/2), 0, 0, sqrt(1/2), 0)
lambda_1 <- eigenvalues[1]
lambda_2 <- eigenvalues[2]

clr_coords <- lapply(1:n_observations, function(i){
  mean + rnorm(1, 0, sqrt(lambda_1))*pc_1 + rnorm(1, 0, sqrt(lambda_2))*pc_2
})
composition_list <- lapply(clr_coords, clrInv)
x_data <- lapply(1:n_observations, function(i) {
  probs <- composition_list[[i]]
  rmultinom(1, n_counts, probs)[, 1]
})
x_data_matrix <- do.call(rbind, x_data)
return(list("x_data" = x_data,"x_data_matrix" = x_data_matrix))
}


base_setting_vs2 <- function(n_observations = 100,
                            eigenvalues = c(0.6, 0.3, 0.1, 0.1),
                            mean = c(0, 0.01, 0.02, -0.02, -0.01),
                            n_counts = 500) {
pc_1 <- c(0, sqrt(1/2), -sqrt(1/2), 0, 0)
pc_2 <- c(-sqrt(1/2), 0, 0, sqrt(1/2), 0)
pc_3 <- c(sqrt(1/2), 0, -sqrt(1/2), 0, 0)
pc_4 <- c(-sqrt(1/2), 0, 0, 0, sqrt(1/2))

lambda_1 <- eigenvalues[1]
lambda_2 <- eigenvalues[2]
lambda_3 <- eigenvalues[3]
lambda_4 <- eigenvalues[4]

clr_coords <- lapply(1:n_observations, function(i){
  mean + rnorm(1, 0, sqrt(lambda_1))*pc_1 +
      rnorm(1, 0, sqrt(lambda_2)) * pc_2 +
      rnorm(1, 0, sqrt(lambda_3)) * pc_3 +
      rnorm(1, 0, sqrt(lambda_4)) * pc_4
})
composition_list <- lapply(clr_coords, clrInv)
x_data <- lapply(1:n_observations, function(i) {
  probs <- composition_list[[i]]
  rmultinom(1, n_counts, probs)[, 1]
})
x_data_matrix <- do.call(rbind, x_data)
return(list("x_data" = x_data,"x_data_matrix" = x_data_matrix))
}

complex_setting <- function(data,
  eigenvalues = c(0.43, 0.26, 0.12, 0.08),
  mean = c(-1.21, -1.76, 1.70, -0.6,
           1.04, -3.62, -1.26, 1.04, -0.83,
           0.50, 4.0, -0.49, 1.51),
  n_observations = nrow(data)) {

  mean_data <- mean(data$aggregate)
  sd_data <- sd(data$aggregate)

  n_counts <- round(rnorm(n_observations, mean = mean_data, sd = sd_data))

  pc_1 <- c(0.4, 0, 0.3, -0.3, 0, 0, -0.2, -0.2, 0.3, -0.1, 0.3, -0.3, -0.2)
  pc_2 <- c(0.7, 0, -0.3, 0.2, 0, -0.2, 0, 0, 0, 0, -0.4, 0, 0)
  pc_3 <- c(0.2, 0, -0.3, 0, -0.4, 0.6, -0.3, 0.2, 0, 0, 0, 0, 0)
  pc_4 <- c(0.3, 0, 0, 0, 0, 0.2, 0, 0, -0.8, 0, 0.3, 0, 0)

  pc_1_norm <- normalize(pc_1)
  pc_2_norm <- normalize(pc_2)
  pc_3_norm <- normalize(pc_3)
  pc_4_norm <- normalize(pc_4)

  lambda_1 <- eigenvalues[1]
  lambda_2 <- eigenvalues[2]
  lambda_3 <- eigenvalues[3]
  lambda_4 <- eigenvalues[4]
  
  clr_coords <- lapply(1:n_observations, function(i){
    mean + rnorm(1, 0, sqrt(lambda_1)) * pc_1_norm +
      rnorm(1, 0, sqrt(lambda_2)) * pc_2_norm +
      rnorm(1, 0, sqrt(lambda_3)) * pc_3_norm +
      rnorm(1, 0, sqrt(lambda_4)) * pc_4_norm
  })

  composition_list <- lapply(clr_coords, clrInv)

  x_data <- lapply(1:n_observations, function(i) {
    probs <- composition_list[[i]]
    rmultinom(1, n_counts[i] , probs)[, 1]
  })
  x_data_matrix <- do.call(rbind, x_data)

  return(list("x_data" = x_data,"x_data_matrix" = x_data_matrix))
}

complex_setting_auto <- function(data,
  eigenvalues = c(0.43, 0.26, 0.12, 0.08),
  mean = c(-1.21, -1.76, 1.70, -0.6,
           1.04, -3.62, -1.26, 1.04, -0.83,
           0.50, 4.0, -0.49, 1.51),
  n_observations = nrow(data),
  n_periods = 10) {

  mean_data <- mean(data$aggregate)
  sd_data <- sd(data$aggregate)

  n_counts <- round(rnorm(n_observations, mean = mean_data, sd = sd_data))

  pc_1 <- c(0.4, 0, 0.3, -0.3, 0, 0, -0.2, -0.2, 0.3, -0.1, 0.3, -0.3, -0.2)
  pc_2 <- c(0.7, 0, -0.3, 0.2, 0, -0.2, 0, 0, 0, 0, -0.4, 0, 0)
  pc_3 <- c(0.2, 0, -0.3, 0, -0.4, 0.6, -0.3, 0.2, 0, 0, 0, 0, 0)
  pc_4 <- c(0.3, 0, 0, 0, 0, 0.2, 0, 0, -0.8, 0, 0.3, 0, 0)

  pc_1_norm <- normalize(pc_1)
  pc_2_norm <- normalize(pc_2)
  pc_3_norm <- normalize(pc_3)
  pc_4_norm <- normalize(pc_4)

  lambda_1 <- eigenvalues[1]
  lambda_2 <- eigenvalues[2]
  lambda_3 <- eigenvalues[3]
  lambda_4 <- eigenvalues[4]

  coef1 <- generate_temporal_scores(n_observations, n_periods, 
                                    eigenvalues[1], phase = 0)
  coef2 <- generate_temporal_scores(n_observations, n_periods, 
                                    eigenvalues[2], phase = pi / 4)
  coef3 <- generate_temporal_scores(n_observations, n_periods, 
                                    eigenvalues[3], phase = pi / 2)
  coef4 <- generate_temporal_scores(n_observations, n_periods, 
                                    eigenvalues[4], phase = pi / 3)
  
  clr_coords <- lapply(1:n_observations, function(i){
    mean + coef1[i] * pc_1_norm + coef2[i]*pc_2_norm +
      coef3[i] * pc_3_norm + coef4[i] * pc_4_norm
  })

  composition_list <- lapply(clr_coords, clrInv)

  x_data <- lapply(1:n_observations, function(i) {
    probs <- composition_list[[i]]
    rmultinom(1, n_counts[i] , probs)[, 1]
  })
  x_data_matrix <- do.call(rbind, x_data)

  return(list("x_data" = x_data,"x_data_matrix" = x_data_matrix))
}

generate_temporal_scores <- function(n_observations, n_periods, eigenvalue, phase = 0) {
  base_scores <- rnorm(n_observations, mean = 0, sd = sqrt(eigenvalue))
  t <- seq(0, 2*pi*n_periods, length.out = n_observations)
  temporal_component <- sin(t + phase)
  scores <- base_scores + temporal_component * sqrt(eigenvalue)
  return(scores)
}

run_composition_simulation <- function(n_simulations, n_observations, n_counts, balanced_setting_fn) {
  sim_composition_results <- vector("list", n_simulations)

  for (i in 1:n_simulations) {
    set.seed(i)
    sim_composition_results[[i]] <- balanced_setting_fn(
      n_observations = n_observations,
      n_counts = n_counts
    )
  }

  pca_results_list <- vector("list", length(sim_composition_results))
  
  for (l in seq_along(sim_composition_results)) {
    sim <- sim_composition_results[[l]]
    x_data <- sim$x_data
    
    set.seed(l)
    result <- co_pca_mcem_nograd(
      x_data,
      lambda = 1,
      max_iter = 40,
      eps = 0.05,
      sum_exp = TRUE
    )
    pca_results_list[[l]] <- result
  }
  
  return(list(
    simulations = sim_composition_results,
    pca_results = pca_results_list
  ))
}

run_complex_simulation <- function(n_simulations, n_observations, data, complex_setting_fn) {
  sim_composition_results <- vector("list", n_simulations)

  for (i in 1:n_simulations) {
    set.seed(i)
    sim_composition_results[[i]] <- complex_setting_fn(
      data = data,
      n_observations = n_observations
    )
  }

  pca_results_list <- vector("list", length(sim_composition_results))
  for (l in seq_along(sim_composition_results)) {
    sim <- sim_composition_results[[l]]
    x_data <- sim$x_data
    set.seed(l)
    result <- co_pca_mcem_nograd(
      x_data,
      lambda = 1,
      max_iter = 40,
      eps = 0.12,
      sum_exp = TRUE
    )
    pca_results_list[[l]] <- result
  }

  return(list(
    simulations = sim_composition_results,
    pca_results = pca_results_list
  ))
}

calculate_pca <- function(simulation_data) {

  pca_results <- lapply(simulation_data$simulations, function(sim) {
    x_data <- do.call(rbind, sim$x_data)
    
    if (any(x_data == 0)) {
      x_data_repl <- cmultRepl(x_data, method = "CZM", suppress.print = TRUE)
    } else {
      x_data_repl <- x_data
    }
    
    list(
      standard = prcomp(clr(x_data_repl)),
      robust = pcaCoDa(x_data_repl, method = "robust")
    )
  })
  
  return(pca_results)
}

calculate_distances <- function(sim_results, pca_results, true_v, true_sigma) {
  n_sims <- length(sim_results$simulations)

  diff_mean_mcem <- sapply(sim_results$pca_results, function(x) 
    sqrt(sum((true_mu - x$pca$center)^2)))

  diff_mean_std <- sapply(pca_results, function(x) 
    sqrt(sum((true_mu - x$standard$center)^2)))
  
  sigma_dist_mcem <- sapply(sim_results$pca_results, function(x) {
    sigma_hat <- with(x$pca, rotation %*% diag(eigenvalues) %*% t(rotation))
    norm(true_sigma - sigma_hat, type = "F")
  })
  
  sigma_dist_std <- sapply(pca_results, function(x) {
    sigma_hat <- with(x$standard, rotation %*% diag(sdev^2) %*% t(rotation))
    norm(true_sigma - sigma_hat, type = "F")
  })

  sigma_dist_rob <- sapply(pca_results, function(x) {
    sigma_hat <- with(x$robust, loadings %*% diag(eigenvalues) %*% t(loadings))
    norm(true_sigma - sigma_hat, type = "F")
  })
  
  return(list(
    mean_mcem = diff_mean_mcem,
    mean_std = diff_mean_std,
    sigma_mcem = sigma_dist_mcem,
    sigma_std = sigma_dist_std,
    sigma_rob = sigma_dist_rob
  ))
}

restructure_results <- function(simulation_data_list, simulation_results_list) {
  if (length(simulation_data_list) != length(simulation_results_list)) {
    stop("The lengths of simulation_data_list and simulation_results_list must be equal.")
  }

  combined_list <- lapply(seq_along(simulation_data_list), function(i) {
    list(
      simulations = simulation_data_list[[i]],
      pca_results = simulation_results_list[[i]]
    )
  })

  return(combined_list)
}