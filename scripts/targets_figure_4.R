required_packages <- c("ggplot2", "targets")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

source("scripts/helper_functions.R")
source("scripts/read_data_KL15_XRF.R")

simulation_results_list <- list(
  tar_read(result_balanced_m20),
  tar_read(result_balanced_m40),
  tar_read(result_balanced_m80),
  tar_read(result_balanced_m160) 
)

n_simulations <- length(simulation_results_list[[1]])

ess_data <- data.frame(
  ESS = c(
    simulation_results_list[[1]][[1]]$ESS / (simulation_results_list[[1]][[1]]$iteration * 10),
    simulation_results_list[[2]][[1]]$ESS / (simulation_results_list[[1]][[1]]$iteration * 10),
    simulation_results_list[[3]][[1]]$ESS / (simulation_results_list[[1]][[1]]$iteration * 10),
    simulation_results_list[[4]][[1]]$ESS / (simulation_results_list[[1]][[1]]$iteration * 10)
  ),
  m = factor(rep(c("20", "40", "80", "160"), each = 100))
)

ess_data$m <- factor(ess_data$m, 
                        levels = c("20", "40", "80", "160"))

plot1 <- ggplot(ess_data, aes(x = ESS, fill = m)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "dodge", 
                 bins = 30, alpha = 0.7) +
  geom_density(aes(color = m), fill = NA, linewidth = 0.5, adjust = 0.5) +
  scale_fill_manual(values = c("#440154FF", "#31688EFF", "#9ebcda", "#e0ecf4")) +
  scale_color_manual(values = c("#440154FF", "#31688EFF", "#9ebcda", "#e0ecf4")) +
  theme_minimal() +
  labs(
    title = "Effective sample sizes for important sampling",
    x = "Percentage of effective weights over the number of replicates",
    y = ""
  )

png("./scripts/figures/figure_4.png", width = 12, height = 5, units = "in", res = 300)
print(plot1)
dev.off()