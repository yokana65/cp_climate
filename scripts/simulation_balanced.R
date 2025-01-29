balanced_setting <- function(n_observations = 100,
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


balanced_setting_vs2 <- function(n_observations = 100,
                            eigenvalues = c(0.6, 0.3, 0.1, 0.1),
                            mean = c(0, 0.01, 0.02, -0.02, -0.01),
                            n_counts = 500) {
pc_1 <- c(0, sqrt(1/2), -sqrt(1/2), 0, 0)
pc_2 <- c(-sqrt(1/2), 0, 0, sqrt(1/2), 0)
pc_3 <- c(sqrt(1/2), 0, -sqrt(1/2), 0, 0)
pc_4 <- c(-sqrt(1/2), 0, 0, 0, sqrt(1/2))
lambda_1 <- eigenvalues[1]
lambda_2 <- eigenvalues[2]
lambda_1 <- eigenvalues[3]
lambda_2 <- eigenvalues[4]

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

# alternative functional version
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
    result <- co_pca_mcem(
      x_data,
      lambda = 1,
      max_iter = 40,
      eps = 0.02,
      sum_exp = TRUE
    )
    pca_results_list[[l]] <- result
  }
  
  return(list(
    simulations = sim_composition_results,
    pca_results = pca_results_list
  ))
}

results <- run_composition_simulation(
  n_simulations = 100,
  n_observations = 100,
  n_counts = 20,
  balanced_setting
)