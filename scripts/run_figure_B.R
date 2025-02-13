required_packages <- c("ggplot2", "dplyr")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

source("scripts/read_data_KL15_XRF.R")

dir <- paste0(getwd(), "/data/Africa_NE_200/data/")
data_kl15_xrf <- read.table(paste0(dir,'data_KL15_XRF.txt'), header = TRUE, sep = "\t")
data_kl15_xrf <- rename(data_kl15_xrf, depth = User_ID)
data_kl15_agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
data_kl15_agem <- rename(data_kl15_agem, age = best)
results <- read_data_kl15_xrf(data_kl15_xrf, data_kl15_agem)
data <- results$data_kl15

plot1 <- ggplot(data.frame(values = data$aggregate), aes(x = values)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, 
                 fill = "#9ebcda", color = "black", alpha = 0.7) +
  geom_density(color = "#8856a7", linewidth = 1) +
  theme_minimal() +
  labs(title = "Distribution of sample sizes for 2119 observations in the KL15 dataset",
       x = "Values",
       y = "Density") +
  theme(plot.title = element_text(hjust = 0.5))
  
png("./scripts/figures/figure_B.png", width = 12, height = 5, units = "in", res = 300)
print(plot1)
dev.off()