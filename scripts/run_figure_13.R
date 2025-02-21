required_packages <- c("ggplot2", "gridExtra", "grid", "compositions")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

source("scripts/read_data_KL15_XRF.R")
source("scripts/help_functions.R")

dir <- paste0(getwd(), "/data/Africa_NE_200/data/")
data_kl15_xrf <- read.table(paste0(dir,'data_KL15_XRF.txt'), header = TRUE, sep = "\t")
data_kl15_xrf <- rename(data_kl15_xrf, depth = User_ID)
data_kl15_agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
data_kl15_agem <- rename(data_kl15_agem, age = best)
results <- read_data_kl15_xrf(data_kl15_xrf, data_kl15_agem)
data_kl15 <- results$data_kl15

data_sel <- data_kl15[4:ncol(data_kl15)-1]
x_ilr <- ilr(data_sel)

colnames(x_ilr) <- c("coordinate 1", "coordinate 2", "coordinate 3", "coordinate 4", 
                        "coordinate 5", "coordinate 6", "coordinate 7", "coordinate 8", 
                        "coordinate 9", "coordinate 10", "coordinate 11", "coordinate 12")

x_ilr_sc <- scale(x_ilr, center=TRUE, scale=FALSE)

#*******reproduction figure 13***********#
create_qq_plot <- function(data, var_name) {
  ggplot(data.frame(x = data), aes(sample = x)) +
    stat_qq() +
    stat_qq_line() +
    ggtitle(paste("Q-Q Plot:", var_name)) +
    theme_minimal()
}

qq_plots <- lapply(1:ncol(x_ilr_sc), function(i) {
  create_qq_plot(x_ilr_sc[,i], colnames(x_ilr_sc)[i])
})

png("./scripts/figures/figure_13.png", width = 12, height = 8, units = "in", res = 300)
grid.arrange(grobs = qq_plots, ncol = 4)
dev.off()

#*******reproduction table ??***********#
descriptives <- summary_stats(data_sel)
descriptives <- summary_stats(round(data_sel * 0.01))