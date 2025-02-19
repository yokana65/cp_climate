required_packages <- c("ggplot2", "targets", "zCompositions", "gridExtra")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

source("scripts/read_data_KL15_XRF.R")
source("scripts/help_functions.R")
source("scripts/fit_comp_pca.R")
source("scripts/cond_scores_function.R")
source("scripts/grad_function.R")

dir <- paste0(getwd(), "/data/Africa_NE_200/data/")
data_kl15_xrf <- read.table(paste0(dir,'data_KL15_XRF.txt'), header = TRUE, sep = "\t")
data_kl15_xrf <- rename(data_kl15_xrf, depth = User_ID)
data_kl15_agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
data_kl15_agem <- rename(data_kl15_agem, age = best)
results <- read_data_kl15_xrf(data_kl15_xrf, data_kl15_agem)
data_kl15 <- results$data_kl15

data_sel <- round(data_kl15[4:ncol(data_kl15)-1] * 0.01)
colnames(data_sel) <- gsub("_cts", "", colnames(data_sel))
data_repl <- cmultRepl(data_sel, method = "CZM", suppress.print=TRUE)
x_clr <- clr(data_repl)

pca_classic <- prcomp(x_clr)

block_bootstrap_pca <- function(data, block_length = 50, R = 1000, reference_pca = NULL) {
  n <- nrow(data)
  results <- list()
  
  # Calculate reference PCA if not provided
  if (is.null(reference_pca)) {
    reference_pca <- prcomp(data, scale. = FALSE)
  }
  
  for(i in 1:R) {
    # Generate blocks
    blocks <- ceiling(n/block_length)
    start_points <- sample(1:(n - block_length + 1), blocks, replace = TRUE)
    
    # Create bootstrap sample
    boot_indices <- unlist(lapply(start_points, function(x) x:(x + block_length - 1)))
    boot_indices <- boot_indices[boot_indices <= n][1:n]
    
    # Perform PCA on bootstrap sample
    boot_data <- data[boot_indices, ]
    boot_pca <- prcomp(boot_data, scale. = FALSE)
    
    # Align signs with reference PCA
    for(j in 1:ncol(boot_pca$rotation)) {
      # Calculate correlation with reference
      cor_sign <- sign(sum(boot_pca$rotation[,j] * reference_pca$rotation[,j] ))
      # Flip sign if negative correlation
      if(cor_sign < 0) {
        boot_pca$rotation[,j] <- -boot_pca$rotation[,j]
        boot_pca$x[,j] <- -boot_pca$x[,j]
      }
    }
    
    results[[i]] <- list(
      rotation = boot_pca$rotation,
      sdev = boot_pca$sdev
    )
  }
  
  return(results)
}

# Anwendung:
bootstrap_results <- block_bootstrap_pca(x_clr, 
                                       block_length = 50, 
                                       R = 1000, 
                                       reference_pca = pca_classic)

# Calculate confidence intervals for loadings
loading_ci <- apply(simplify2array(lapply(bootstrap_results, function(x) x$rotation)), 1:2, 
                   quantile, probs = c(0.025, 0.975))

# direction is changed to be consistent with the mcem results
loading_df <- data.frame(
  Variable = rownames(pca_classic$rotation),
  Loading = pca_classic$rotation[,1] * (-1),
  Lower = loading_ci[1,,1] * (-1),
  Upper = loading_ci[2,,1] * (-1)
)

plot1 <- ggplot(loading_df, aes(x = Variable, y = Loading)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper)) +
  coord_flip() +
  labs(title = "PC1 loading intervals for the standard estimation")

btst_results <- tar_read(bbootstrap_01_3)
mcem_result <- tar_read(results_mcem)
# mcem_result <- co_pca_mcem_nograd(
#           data_sel, 
#           lambda = 1,
#           max_iter = 40,
#           eps = 0.15,
#           sum_exp = FALSE
#       )

# # Calculate confidence intervals for loadings
loading_ci <- apply(simplify2array(lapply(btst_results, function(x) 
                  x$rotation)), 1:2, 
                   quantile, probs = c(0.025, 0.975))

loading_df <- data.frame(
  Variable = rownames(pca_classic$rotation),
  Loading = mcem_result$pca$rotation[,1],
  Lower = loading_ci[1,,1],
  Upper = loading_ci[2,,1]
)

plot2 <- ggplot(loading_df, aes(x = Variable, y = Loading)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper)) +
  coord_flip() +
  labs(title = "PC1 loading intervals for mcem estimation")

png("./scripts/figures/figure_E.png", width = 12, height = 8, units = "in", res = 300)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()