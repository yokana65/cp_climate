library(ggplot2)
library(gridExtra)
library(compositions)
library(robCompositions)
library(zCompositions)

source("scripts/fit_composition_pca.R")
source("scripts/helper_functions.R")
source("scripts/conditional_scores_function.R")
source("scripts/gradient.R")

set_1 <- c("#e0ecf4", "#9ebcda", "#8856a7")

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

calculate_pca <- function(simulation_data) {
  pca_clr_std <- list()
  pca_clr_rob <- list()
  
  for(i in 1:length(simulation_data$simulations)) {
    x_data <- do.call(rbind, simulation_data$simulations[[i]]$x_data)
    
    if (any(x_data == 0)) {
      x_data_repl <- cmultRepl(x_data, method = "CZM", suppress.print=TRUE)
    } else {
      x_data_repl <- x_data
    }
    
    pca_clr_std[[i]] <- prcomp(clr(x_data_repl))
    pca_clr_rob[[i]] <- pcaCoDa(x_data_repl, method = "robust")
  }
  
  return(list(standard = pca_clr_std, robust = pca_clr_rob))
}

calculate_distances <- function(sim_results, pca_results, true_v, true_sigma) {
  n_sims <- length(sim_results$simulations)
  
  # Mean distances
  diff_mean_mcem <- sapply(sim_results$pca_results, function(x) 
    sqrt(sum((true_v - x$center)^2)))
  
  diff_mean_std <- sapply(pca_results$standard, function(x) 
    sqrt(sum((true_v - x$center)^2)))
  
  # Covariance distances
  sigma_dist_mcem <- sapply(sim_results$pca_results, function(x) {
    sigma_hat <- with(x, rotation %*% diag(eigenvalues) %*% t(rotation))
    norm(true_sigma - sigma_hat, type = "F")
  })
  
  sigma_dist_std <- sapply(pca_results$standard, function(x) {
    sigma_hat <- with(x, rotation %*% diag(sdev^2) %*% t(rotation))
    norm(true_sigma - sigma_hat, type = "F")
  })
  
  return(list(
    mean_mcem = diff_mean_mcem,
    mean_std = diff_mean_std,
    sigma_mcem = sigma_dist_mcem,
    sigma_std = sigma_dist_std
  ))
}

#*******reproduction plo1***********#
n_simulations <- 100
n_observations <- 100
counts_list <- c(20, 40, 80, 160)

true_mu <- c(0, 0.01, 0.02, -0.02, -0.01) 
eigenvalues <- c(0.7, 0.3, 0, 0)
V <- cbind(c(0, sqrt(1/2), -sqrt(1/2), 0, 0), c(-sqrt(1/2), 0, 0, sqrt(1/2), 0), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
true_sigma <- V %*% diag(eigenvalues) %*% t(V) 

simulation_results <- lapply(counts_list, function(n_counts) {
  run_composition_simulation(
    n_simulations = n_simulations,
    n_observations = n_observations,
    n_counts = n_counts,
    balanced_setting_vs2
  )
})

pca_results <- lapply(simulation_results, calculate_pca)

distances <- lapply(1:length(counts_list), function(i) {
  calculate_distances(
    simulation_results[[i]], 
    pca_results[[i]], 
    true_mu, 
    true_sigma
  )
})

mean_df <- data.frame(
  differences = c(
    unlist(lapply(distances, `[[`, "mean_mcem")),
    unlist(lapply(distances, `[[`, "mean_std"))
  ),
  method = factor(rep(c("MCEM", "sample mean"), 
                     each = n_simulations * length(counts_list))),
  size = factor(rep(rep(as.character(counts_list), 
                       each = n_simulations), 2))
)

plot1 <- ggplot(mean_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_grey() +
  scale_fill_manual(values = set_1) +
  labs(y = "dist mean",
       x = "Scale",
       fill = "method",
       title = "Distance of true mean to estimated mean") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 10))

cov_df <- data.frame(
  differences = c(
    unlist(lapply(distances, `[[`, "sigma_mcem")),
    unlist(lapply(distances, `[[`, "sigma_std")),
    unlist(lapply(distances, `[[`, "sigma_rob"))  # Robuste Schätzungen hinzugefügt
  ),
  method = factor(rep(c("MCEM", "standard", "robust"), 
                     each = n_simulations * length(counts_list))),
  size = factor(rep(rep(as.character(counts_list), 
                       each = n_simulations), 3))
)

plot2 <- ggplot(cov_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_grey() +
  scale_fill_manual(values = set_1) +
  labs(y = "dist covariance",
       x = "Scale",
       fill = "method",
       title = "Distance of true covariance to estimated covariance") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 10))

grid.arrange(plot1, plot2, ncol = 2)