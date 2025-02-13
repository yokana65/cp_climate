source("scripts/helper_functions.R")

required_packages <- c("ggplot2", "gridExtra", "grid", "compositions", 
                      "robCompositions", "zCompositions", "targets", "dplyr")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

set_1 <- get_color_palette()

#*******reproduction figure A1 with targets***********#
eigenvalues <- c(0.43, 0.26, 0.12, 0.08)
true_mu <- c(-1.21, -1.76, 1.70, -0.6,
                                 1.04, -3.62, -1.26, 1.04, -0.83,
                                 0.50, -4.0, -0.49, 1.51)

pc_1 <- c(0.4, 0, 0.3, -0.3, 0, 0, -0.2, -0.2, 0.3, -0.1, 0.3, -0.3, -0.2)
pc_2 <- c(0.7, 0, -0.3, 0.2, 0, -0.2, 0, 0, 0, 0, -0.4, 0, 0)
pc_3 <- c(0.2, 0, -0.3, 0, -0.4, 0.6, -0.3, 0.2, 0, 0, 0, 0, 0)
pc_4 <- c(0.3, 0, 0, 0, 0, 0.2, 0, 0, -0.8, 0, 0.3, 0, 0)

pc_1_norm <- normalize(pc_1)
pc_2_norm <- normalize(pc_2)
pc_3_norm <- normalize(pc_3)
pc_4_norm <- normalize(pc_4)

V <- cbind(pc_1_norm, pc_2_norm, pc_3_norm, pc_4_norm)
true_sigma <- V %*% diag(eigenvalues) %*% t(V)

sample_sizes <- c("cts per sec", "cts per dsec", "cts per csec")
methods <- c("MCEM", "standard", "robust")

simulation_data_list <- list(
  tar_read(sim_complex_sec),
  tar_read(sim_complex_dsec),
  tar_read(sim_complex_csec)
)

n_simulations <- length(simulation_data_list[[1]])

simulation_results_list <- list(
  tar_read(result_complex_sec),
  tar_read(result_complex_dsec),
  tar_read(result_complex_csec)
)

pca_clr_rob <- lapply(seq_along(methods), function(x) list(n_simulations))
pca_clr_std <- lapply(seq_along(methods), function(x) list(n_simulations))

for(j in seq_along(methods)) {
  simulation_data <- simulation_data_list[[j]]
  for(i in seq_along(simulation_data)) {
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

diff_mean <- lapply(seq_along(methods), function(x) list(n_simulations))
diff_mean_std_clr <- lapply(seq_along(methods), function(x) list(n_simulations))

for(j in seq_along(methods)) {
  for(i in 1:n_simulations) {
    diff_mean[[j]][[i]] <- sqrt(sum((true_mu - simulation_results_list[[j]][[i]]$pca$center)^2))
    diff_mean_std_clr[[j]][[i]] <- sqrt(sum((true_mu - pca_clr_std[[j]][[i]]$center)^2))
  }
}

diff_df <- data.frame(
  differences = c(
    unlist(diff_mean[[1]]), unlist(diff_mean[[2]]), 
    unlist(diff_mean[[3]]),
    unlist(diff_mean_std_clr[[1]]), unlist(diff_mean_std_clr[[2]]),
    unlist(diff_mean_std_clr[[3]])
  ),
  method = factor(rep(c("MCEM", "sample mean"), each = length(sample_sizes) * n_simulations)),
  size = factor(rep(rep(sample_sizes, 
                         each = n_simulations), 2))
)

diff_df$method <- factor(diff_df$method, 
                        levels = c("MCEM", "sample mean"))

diff_df$size <- factor(diff_df$size, 
                        levels = sample_sizes)

plot1 <- ggplot(diff_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = set_1) + 
  labs(y = "dist mean",
       x = "Scale",
       fill = "method",
       title = "Distance of true mean to estimated mean") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16))

sigma_distances <- lapply(seq_along(methods), function(x) list(n_simulations))
sigma_distances_std_clr <- lapply(seq_along(methods), function(x) list(n_simulations))
sigma_distances_rob_clr <- lapply(seq_along(methods), function(x) list(n_simulations))

for(j in seq_along(methods)) {
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
    unlist(sigma_distances[[3]]), 
    unlist(sigma_distances_std_clr[[1]]), unlist(sigma_distances_std_clr[[2]]),
    unlist(sigma_distances_std_clr[[3]]), 
    unlist(sigma_distances_rob_clr[[1]]), unlist(sigma_distances_rob_clr[[2]]),
    unlist(sigma_distances_rob_clr[[3]])
  ),
  method = factor(rep(c("MCEM", "standard with CZM", "robust with CZM"), each = length(methods) * n_simulations)),
  size = factor(rep(rep(sample_sizes, 
                         each = n_simulations), length(sample_sizes)))
)

sigma_diff_df$method <- factor(sigma_diff_df$method, 
                        levels = c("MCEM", "standard with CZM", "robust with CZM"))

sigma_diff_df$size <- factor(sigma_diff_df$size, 
                        levels = sample_sizes)

plot2 <- ggplot(sigma_diff_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = set_1) + 
  labs(y = "dist covariance",
       x = "Scale",
       fill = "method") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16))

#**********zero imputation with GBM method**********#
simulation_data_list <- list(
  tar_read(sim_complex_sec),
  tar_read(sim_complex_dsec),
  tar_read(sim_complex_csec)
)

n_simulations <- length(simulation_data_list[[1]])

simulation_results_list <- list(
  tar_read(result_complex_sec),
  tar_read(result_complex_dsec),
  tar_read(result_complex_csec)
)

pca_clr_rob <- lapply(seq_along(methods), function(x) list(n_simulations))
pca_clr_std <- lapply(seq_along(methods), function(x) list(n_simulations))

for(j in seq_along(methods)) {
  simulation_data <- simulation_data_list[[j]]
  for(i in 1:length(simulation_data)) {
    sim <- simulation_data[[i]]
    x_data_list <- sim$x_data
    x_data <- do.call(rbind, x_data_list)
    if (any(x_data == 0)) {
      x_data_repl <- cmultRepl(x_data, method = "GBM", suppress.print=TRUE)
    } else {
      x_data_repl <- x_data
    }

    pca_clr_std[[j]][[i]] <- prcomp(clr(x_data_repl))
    pca_clr_rob[[j]][[i]] <- pcaCoDa(x_data_repl, method = "robust")
  }
}

sigma_distances <- lapply(seq_along(methods), function(x) list(n_simulations))
sigma_distances_std_clr <- lapply(seq_along(methods), function(x) list(n_simulations))
sigma_distances_rob_clr <- lapply(seq_along(methods), function(x) list(n_simulations))

for(j in seq_along(methods)) {
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
    unlist(sigma_distances[[3]]), 
    unlist(sigma_distances_std_clr[[1]]), unlist(sigma_distances_std_clr[[2]]),
    unlist(sigma_distances_std_clr[[3]]), 
    unlist(sigma_distances_rob_clr[[1]]), unlist(sigma_distances_rob_clr[[2]]),
    unlist(sigma_distances_rob_clr[[3]])
  ),
  method = factor(rep(c("MCEM", "standard with GBM", "robust with GBM"), each = length(methods) * n_simulations)),
  size = factor(rep(rep(sample_sizes, 
                         each = n_simulations), length(sample_sizes)))
)

sigma_diff_df$method <- factor(sigma_diff_df$method, 
                        levels = c("MCEM", "standard with GBM", "robust with GBM"))

sigma_diff_df$size <- factor(sigma_diff_df$size, 
                        levels = sample_sizes)

plot3 <- ggplot(sigma_diff_df, aes(x = size, y = differences, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = set_1) + 
  labs(y = "dist covariance",
       x = "Scale",
       fill = "method") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16))

png("./scripts/figures/figure_A.png", width = 12, height = 5, units = "in", res = 300)
grid.arrange(plot2, plot3, ncol = 2, top = textGrob("Distance of true covariance to estimated covariance for two different imputation methods", gp = gpar(fontsize = 16)))
dev.off()