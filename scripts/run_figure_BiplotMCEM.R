source("scripts/read_data_KL15_XRF.R")
source("scripts/help_functions.R")
source("scripts/fit_comp_pca.R")
source("scripts/cond_scores_function.R")
source("scripts/grad_function.R")

set_1 <- get_color_palette()

dir <- paste0(getwd(), "/data/Africa_NE_200/data/")
data_kl15_xrf <- read.table(paste0(dir,'data_KL15_XRF.txt'), header = TRUE, sep = "\t")
data_kl15_xrf <- rename(data_kl15_xrf, depth = User_ID)
data_kl15_agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
data_kl15_agem <- rename(data_kl15_agem, age = best)
results <- read_data_kl15_xrf(data_kl15_xrf, data_kl15_agem)
data_kl15 <- results$data_kl15

#*******reproduction figure 6***********#
data_sel <- data_kl15[4:ncol(data_kl15)-1] * 0.01
colnames(data_sel) <- gsub("_cts", "", colnames(data_sel))

pca_results_mcem <-
  co_pca_mcem_nm(data_sel,
            max_iter = 50,
            r = 10,
            lambda = 1,
            eps = 0.25)

png("./scripts/figures/figure_BiplotMCEM.png", width = 8, height = 8, units = "in", res = 300)
plot_pca_rotations(pca_results_mcem$pca$rotation, components = c(1,2), main = "PC1 vs. PC2", fixed=FALSE, pos_vector = c(4,4,4,2,4,4,2,2,4,2,4,2,2))
dev.off()
