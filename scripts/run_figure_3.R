required_packages <- c("ggplot2", "gridExtra", "grid", "compositions", 
                      "robCompositions", "zCompositions", "targets", "dplyr")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

set_1 <- get_color_palette()

source("scripts/fit_compositional_pca.R")
source("scripts/helper_functions.R")
source("scripts/conditional_scores_function.R")
source("scripts/gradient.R")
source("scripts/simulation_functions.R")
source("scripts/read_data_KL15_XRF.R")

# #*******reproduction plot3***********#
dir <- paste0(getwd(), "/data/Africa_NE_200/data/")
data_kl15_xrf <- read.table(paste0(dir,'data_KL15_XRF.txt'), header = TRUE, sep = "\t")
data_kl15_xrf <- rename(data_kl15_xrf, depth = User_ID)
data_kl15_agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
data_kl15_agem <- rename(data_kl15_agem, age = best)
results <- read_data_kl15_xrf(data_kl15_xrf, data_kl15_agem)
data <- results$data_kl15

# data <- tar_read(data_kl15)

n_simulations <- 5
n_observations <- 100
scales <- c(1, 0.1, 0.01)

eigenvalues <- c(0.43, 0.26, 0.12, 0.08)
true_mu <- c(-1.21, -1.76, 1.70, -0.6,
                                 1.04, -3.62, -1.26, 1.04, -0.83,
                                 0.50, -4.0, -0.49, 1.51)

pc_1 <- c(0.4, 0, 0.3, -0.3, 0, 0, -0.2, -0.2, 0.3, -0.1, 0.3, -0.3, -0.2)
pc_2 <- c(0.7, 0, -0.3, 0.2, 0, -0.2, 0, 0, 0, 0, -0.4, 0, 0)
pc_3 <- c(0.2, 0, -0.3, 0, -0.4, 0.6, -0.3, 0.2, 0, 0, 0, 0, 0)
pc_4 <- c(0.3, 0, 0, 0, 0, 0.2, 0, 0, -0.8, 0, 0.3, 0, 0)

# normalize <- function(v) {
#   v / sqrt(sum(v^2))
# }

pc_1_norm <- normalize(pc_1)
pc_2_norm <- normalize(pc_2)
pc_3_norm <- normalize(pc_3)
pc_4_norm <- normalize(pc_4)

V <- cbind(pc_1_norm, pc_2_norm, pc_3_norm, pc_4_norm)
true_sigma <- V %*% diag(eigenvalues) %*% t(V) 

simulation_results <- lapply(scales, function(scale) {
  run_complex_simulation(
    n_simulations = n_simulations,
    n_observations = n_observations,
    data = data * scale,
    complex_setting
  )
})

pca_results <- lapply(simulation_results, calculate_pca)

distances <- lapply(1:length(scales), function(i) {
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
                     each = n_simulations * length(scales))),
  size = factor(rep(rep(c("cts per sec", "cts per dsec", "cts per csec"), 
                       each = n_simulations), 2))
)

mean_df$method <- factor(mean_df$method, 
                        levels = c("MCEM", "sample mean"))

mean_df$size <- factor(mean_df$size, 
                      levels = c("cts per sec", "cts per dsec", "cts per csec"))

plot1 <- ggplot(mean_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = set_1) +
  labs(y = "dist mean",
       x = "Scale",
       fill = "method") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16))

cov_df <- data.frame(
  differences = c(
    unlist(lapply(distances, `[[`, "sigma_mcem")),
    unlist(lapply(distances, `[[`, "sigma_std")),
    unlist(lapply(distances, `[[`, "sigma_rob")) 
  ),
  method = factor(rep(c("MCEM", "standard", "robust"), 
                     each = n_simulations * length(scales))),
  size = factor(rep(rep(c("cts per sec", "cts per dsec", "cts per csec"), 
                       each = n_simulations), 3))
)


cov_df$method <- factor(cov_df$method, 
                        levels = c("MCEM", "standard", "robust"))

cov_df$size <- factor(cov_df$size, 
                        levels = c("cts per sec", "cts per dsec", "cts per csec"))

plot2 <- ggplot(cov_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = set_1) +
  labs(y = "dist covariance",
       x = "Scale",
       fill = "method") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16))

#*************computation with autocorrelated scores***************#
simulation_results_auto <- lapply(scales, function(scale) {
  run_complex_simulation(
    n_simulations = n_simulations,
    n_observations = n_observations,
    data = data * scale,
    complex_setting_auto
  )
})

pca_results_auto <- lapply(simulation_results_auto, calculate_pca)

distances_auto <- lapply(1:length(scales), function(i) {
  calculate_distances(
    simulation_results_auto[[i]],
    pca_results_auto[[i]],
    true_mu,
    true_sigma
  )
})

cov_df_auto <- data.frame(
  differences = c(
    unlist(lapply(distances_auto, `[[`, "sigma_mcem")),
    unlist(lapply(distances_auto, `[[`, "sigma_std")),
    unlist(lapply(distances_auto, `[[`, "sigma_rob")) 
  ),
  method = factor(rep(c("MCEM", "standard", "robust"), 
                     each = n_simulations * length(scales))),
  size = factor(rep(rep(c("cts per sec", "cts per dsec", "cts per csec"), 
                       each = n_simulations), 3))
)


cov_df_auto$method <- factor(cov_df$method, 
                        levels = c("MCEM", "standard", "robust"))

cov_df_auto$size <- factor(cov_df$size, 
                        levels = c("cts per sec", "cts per dsec", "cts per csec"))

plot3 <- ggplot(cov_df_auto, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = set_1) +
  labs(y = "dist covariance",
       x = "Scale",
       fill = "method") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16))

grid.arrange(plot2, plot3, ncol = 2, top = textGrob("Distance of true covariance to estimated covariance", gp = gpar(fontsize = 16)))

png("./scripts/figures/figure_3_upd.png", width = 12, height = 5, units = "in", res = 300)
grid.arrange(plot2, plot3, ncol = 2, top = textGrob("Distance of true covariance to estimated covariance", gp = gpar(fontsize = 16)))
dev.off()