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

#*******reproduction of figure 4***********#
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

n_simulations <- length(simulation_results_list[[1]])

ess_data <- data.frame(
  ESS = c(
    simulation_results_list[[1]][[1]]$ESS / (simulation_results_list[[1]][[1]]$iteration * 10),
    simulation_results_list[[2]][[1]]$ESS / (simulation_results_list[[1]][[1]]$iteration * 10),
    simulation_results_list[[3]][[1]]$ESS / (simulation_results_list[[1]][[1]]$iteration * 10),
    simulation_results_list[[4]][[1]]$ESS / (simulation_results_list[[1]][[1]]$iteration * 10)
  ),
  size = factor(rep(c("20", "40", "80", "160"), each = 100))
)

ess_data$size <- factor(ess_data$size, 
                        levels = c("20", "40", "80", "160"))

plot1 <- ggplot(ess_data, aes(x = ESS, fill = size)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "dodge", 
                 bins = 30, alpha = 0.7) +
  geom_density(aes(color = size), fill = NA, linewidth = 0.5, adjust = 0.5) +
  scale_fill_manual(values = c("#440154FF", "#31688EFF", "#9ebcda", "#e0ecf4")) +
  scale_color_manual(values = c("#440154FF", "#31688EFF", "#9ebcda", "#e0ecf4")) +
  theme_minimal() +
  labs(
    title = "Effective sample sizes for important sampling",
    x = "Percentage of effective weights over the number of replicates",
    y = ""
  )

png("./scripts/figures/figure_4_upd.png", width = 12, height = 5, units = "in", res = 300)
print(plot1)
dev.off()
