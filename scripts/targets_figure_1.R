source("scripts/helper_functions.R")

load_required_packages()
set_1 <- get_color_palette()

# Reproduction of Plot 1 with targets
true_mu <- c(0, 0.01, 0.02, -0.02, -0.01) 

eigenvalues <- c(0.7, 0.3, 0, 0)
V <- cbind(c(0, sqrt(1/2), -sqrt(1/2), 0, 0), c(-sqrt(1/2), 0, 0, sqrt(1/2), 0), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
true_sigma <- V %*% diag(eigenvalues) %*% t(V) 

sample_sizes <- c("20", "40", "80", "160")
methods <- c("MCEM", "standard", "robust")

simulation_data_list <- list(
  tar_read(sim_balanced_m20),
  tar_read(sim_balanced_m40),
  tar_read(sim_balanced_m80),
  tar_read(sim_balanced_m160)
  )

n_simulations <- length(simulation_data_list[[1]])

simulation_results_list <- list(
  tar_read(results_balanced_m20),
  tar_read(results_balanced_m40),
  tar_read(results_balanced_m80),
  tar_read(results_balanced_m160)  
)

pca_clr_rob <- lapply(1:4, function(x) list(n_simulations))
pca_clr_std <- lapply(1:4, function(x) list(n_simulations))

for(j in 1:4) {
  simulation_data <- simulation_data_list[[j]]
  for(i in 1:length(simulation_data)) {
    sim <- simulation_data[[i]]
    x_data_list <- sim$x_data
    x_data <- do.call(rbind, x_data_list)
    if (any(x_data == 0)) {
      x_data_repl <- cmultRepl(x_data, method = "CZM", suppress.print=TRUE)
    } else {
      x_data_repl <- x_data
    }

    pca_clr_std[[j]][[i]] <- prcomp(clr(x_data_repl))
    pca_clr_rob[[j]][[i]] <- pcaCoDa(x_data_repl, method = "robust")
  }
}

diff_mean <- lapply(1:4, function(x) list(n_simulations))
diff_mean_std_clr <- lapply(1:4, function(x) list(n_simulations))

for(j in 1:4) {
  for(i in 1:n_simulations) {
    diff_mean[[j]][[i]] <- sqrt(sum((true_mu - simulation_results_list[[j]][[i]]$pca$center)^2))
    diff_mean_std_clr[[j]][[i]] <- sqrt(sum((true_mu - pca_clr_std[[j]][[i]]$center)^2))
  }
}

diff_df <- data.frame(
  differences = c(
    unlist(diff_mean[[1]]), unlist(diff_mean[[2]]), 
    unlist(diff_mean[[3]]), unlist(diff_mean[[4]]),
    unlist(diff_mean_std_clr[[1]]), unlist(diff_mean_std_clr[[2]]),
    unlist(diff_mean_std_clr[[3]]), unlist(diff_mean_std_clr[[4]])
  ),
  method = factor(rep(c("MCEM", "sample mean"), each = 4 * n_simulations)),
  size = factor(rep(rep(c("20", "40", "80", "160"), 
                         each = n_simulations), 2))
)

diff_df$method <- factor(diff_df$method, 
                        levels = c("MCEM", "sample mean"))

diff_df$size <- factor(diff_df$size, 
                        levels = c("20", "40", "80", "160"))

plot1 <- ggplot(diff_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = set_1) + 
  labs(y = "dist mean",
       x = expression(m),
       fill = "method",
       title = "Distance of true mean to estimated mean") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16))

sigma_distances <- lapply(1:4, function(x) list(n_simulations))
sigma_distances_std_clr <- lapply(1:4, function(x) list(n_simulations))
sigma_distances_rob_clr <- lapply(1:4, function(x) list(n_simulations))


for(j in 1:4) {
  for(i in 1:n_simulations) {
    sigma_hat_mcem <- with(simulation_results_list[[j]][[i]]$pca,
      rotation %*% diag(eigenvalues) %*% t(rotation))
    sigma_distances[[j]][[i]] <- norm(true_sigma - sigma_hat_mcem, type = "F")
    
    sigma_hat_clr <- with(pca_clr_std[[j]][[i]],
      rotation %*% diag(sdev^2) %*% t(rotation))
    sigma_distances_std_clr[[j]][[i]] <- norm(true_sigma - sigma_hat_clr, type = "F")
    
    sigma_hat_rob <- with(pca_clr_rob[[j]][[i]],
    loadings %*% diag(eigenvalues) %*% t(loadings))
    sigma_distances_rob_clr[[j]][[i]] <- norm(true_sigma - sigma_hat_rob, type = "F")
  }
}

sigma_diff_df <- data.frame(
  differences = c(
    unlist(sigma_distances[[1]]), unlist(sigma_distances[[2]]), 
    unlist(sigma_distances[[3]]), unlist(sigma_distances[[4]]),
    unlist(sigma_distances_std_clr[[1]]), unlist(sigma_distances_std_clr[[2]]),
    unlist(sigma_distances_std_clr[[3]]), unlist(sigma_distances_std_clr[[4]]),
    unlist(sigma_distances_rob_clr[[1]]), unlist(sigma_distances_rob_clr[[2]]),
    unlist(sigma_distances_rob_clr[[3]]), unlist(sigma_distances_rob_clr[[4]])
  ),
  method = factor(rep(c("MCEM", "standard", "robust"), each = 4 * n_simulations)),
  size = factor(rep(rep(c("20", "40", "80", "160"), 
                         each = n_simulations), 3))
)

sigma_diff_df$method <- factor(sigma_diff_df$method, 
                        levels = c("MCEM", "standard", "robust"))

sigma_diff_df$size <- factor(sigma_diff_df$size, 
                        levels = c("20", "40", "80", "160"))

plot2 <- ggplot(sigma_diff_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = set_1) + 
  labs(y = "dist covariance",
       x = expression(m),
       fill = "method",
       title = "Distance of true covariance to estimated covariance") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16))

png("./scripts/figures/figure_1.png", width = 12, height = 5, units = "in", res = 300)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()