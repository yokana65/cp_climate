source("scripts/helper_functions.R")

required_packages <- c("ggplot2", "targets")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

source("scripts/read_data_KL15_XRF.R")

simulation_results_list <- list(
  tar_read(result_complex_csec),
  tar_read(result_complex_csec_auto)
)

n_simulations <- length(simulation_results_list[[1]])

iteration_data <- data.frame(
  iteration = c(
    unlist(lapply(simulation_results_list[[1]], function(x) x$iteration)),
    unlist(lapply(simulation_results_list[[2]], function(x) x$iteration))
  ),
  setting = factor(
    c(
      rep("uncorrelated scores", length(unlist(lapply(simulation_results_list[[1]], function(x) x$iteration)))),
      rep("scores with temporal trend", length(unlist(lapply(simulation_results_list[[2]], function(x) x$iteration))))
    ),
    levels = c("uncorrelated scores", "scores with temporal trend")
  )
)

plot1 <- ggplot(iteration_data, aes(x = iteration, fill = setting)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "dodge", 
                 bins = 30, alpha = 0.7) +
  geom_density(aes(color = setting), fill = NA, linewidth = 0.5, adjust = 1.3) +
  scale_fill_manual(values = c("#e0ecf4", "#9ebcda")) +
  scale_color_manual(values = c("#e0ecf4", "#9ebcda")) +
  theme_minimal() +
  labs(
    title = "Number of iterations until convergence",
    x = "Iterations",
    y = "Distribution"
  )

png("./scripts/figures/figure_C.png", width = 12, height = 5, units = "in", res = 300)
print(plot1)
dev.off()