required_packages <- c("ggplot2", "gridExtra", "grid", "compositions",
                        "robCompositions", "dplyr", "reshape2",
                        "scales", "mgcv")

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
colnames(data_sel) <- gsub("_cts", "", colnames(data_sel))
x_clr <- clr(data_sel)

pca_rob <- pcaCoDa(data_sel)

#*******reproduction figure 9***********#
pc_time_data_list <- lapply(1:4, function(i) {
  data.frame(
    age = data_kl15$age,
    score = calculate_pc_scores(x_clr, pca_rob, i)
  )
})

insolation <- read.table("./Spatio-Temporal-Variations/Linear Analysis/data all/Insolation.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

insolation <- insolation %>%
  rename(age = anu) %>% 
  filter(age <= 547000) %>% 
  mutate(
    age = age * 0.001,
    ins_chb_s = rescale(ins_chb_s, to = c(-1, 1)),
    ins_chb_a = rescale(ins_chb_a, to = c(-1, 1)),
    ins_tibet_s = rescale(ins_tibet_s, to = c(-1, 1))
  ) %>% 
  select(age, ins_tibet_s)

data_long <- melt(insolation, id.vars = "age", 
                  variable.name = "variable", 
                  value.name = "value")

gam_fit1 <- gam(score ~ s(age, k=520, bs="cr"), data=pc_time_data_list[[1]])

pred_data <- data.frame(age = seq(min(pc_time_data_list[[1]]$age), 
                                 max(pc_time_data_list[[1]]$age), 
                                 length.out=200))
pred_data$fit <- predict(gam_fit1, pred_data)

gam_fit2 <- gam(score ~ s(age, k=520, bs="cr"), data=pc_time_data_list[[2]])

pred_data_2 <- data.frame(age = seq(min(pc_time_data_list[[2]]), 
                                 max(pc_time_data_list[[2]]$age), 
                                 length.out=200))
pred_data_2$fit <- predict(gam_fit2, pred_data_2)

pred_se <- predict(gam_fit1, pred_data, se.fit=TRUE)
pred_data$upper <- pred_data$fit + 1.96 * pred_se$se.fit
pred_data$lower <- pred_data$fit - 1.96 * pred_se$se.fit

pred_se_2 <- predict(gam_fit2, pred_data_2, se.fit=TRUE)
pred_data_2$upper <- pred_data_2$fit + 1.96 * pred_se_2$se.fit
pred_data_2$lower <- pred_data_2$fit - 1.96 * pred_se_2$se.fit

plot1 <- ggplot(pc_time_data_list[[1]], aes(x = age, y = score)) +
    geom_point(alpha=0.5) +
    geom_line(data=pred_data, aes(y=fit), color="#8856a7", linewidth=0.6) +
    geom_ribbon(data=pred_data, aes(y=fit, ymin=lower, ymax=upper), 
                alpha=0.2) +
    geom_line(data = data_long, aes(x = age, y = value), color = "#9ebcda", linewidth = 0.8, alpha = 0.8) +
    labs(x = "", 
         y = "PC1 Score") +
    theme_minimal() +
    ylim(-1.5,1.5) +
    theme(legend.position = "none") 

plot2 <- ggplot(pc_time_data_list[[2]], aes(x = age, y = score)) +
    geom_point(alpha=0.5) +
    geom_line(data=pred_data_2, aes(y=fit), color="#9ebcda", linewidth=0.6) +
    geom_ribbon(data=pred_data_2, aes(y=fit, ymin=lower, ymax=upper), 
                alpha=0.2) +
    geom_line(data = data_long, aes(x = age, y = value), color = "#8856a7", linewidth = 0.8, , alpha = 0.8) +
    labs(x = "Age (ka)", 
         y = "PC2 Score") +
    theme_minimal() +
    ylim(-1.5,1.5) +
    theme(legend.position = "none") 

png("./scripts/figures/figure_9.png", width = 12, height = 8, units = "in", res = 300)
grid.arrange(plot1, plot2, ncol = 1, top = textGrob("Temporal trend of principal component scores and precession circles", gp = gpar(fontsize = 16)))
dev.off()