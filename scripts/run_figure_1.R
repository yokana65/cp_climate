required_packages <- c("ggplot2", "gridExtra", "compositions", 
                      "robCompositions", "zCompositions")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

source("scripts/fit_comp_pca.R")
source("scripts/help_functions.R")
source("scripts/cond_scores_function.R")
source("scripts/grad_function.R")
source("scripts/simulation_functions.R")

#*******reproduction of figure 1***********#
n_simulations <- 5
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
    base_setting
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

mean_df$method <- factor(mean_df$method, 
                        levels = c("MCEM", "sample mean"))

mean_df$size <- factor(mean_df$size, 
                        levels = c("20", "40", "80", "160"))

plot1 <- ggplot(mean_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = set_1) +
  labs(y = "dist mean",
       x = expression(m),
       fill = "method",
       title = "Distance of true mean to estimated mean") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16))

cov_df <- data.frame(
  differences = c(
    unlist(lapply(distances, `[[`, "sigma_mcem")),
    unlist(lapply(distances, `[[`, "sigma_std")),
    unlist(lapply(distances, `[[`, "sigma_rob")) 
  ),
  method = factor(rep(c("MCEM", "standard", "robust"), 
                     each = n_simulations * length(counts_list))),
  size = factor(rep(rep(as.character(counts_list), 
                       each = n_simulations), 3))
)


cov_df$method <- factor(cov_df$method, 
                        levels = c("MCEM", "standard", "robust"))

cov_df$size <- factor(cov_df$size, 
                        levels = c("20", "40", "80", "160"))

plot2 <- ggplot(cov_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = set_1) +
  labs(y = "dist covariance",
       x = expression(m),
       fill = "method",
       title = "Distance of true covariance to estimated covariance") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16))

png("./scripts/figures/figure_1_upd.png", width = 12, height = 5, units = "in", res = 300)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()