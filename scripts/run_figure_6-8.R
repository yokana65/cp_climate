required_packages <- c("ggplot2", "gridExtra", "grid", "compositions",
                        "robCompositions")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

source("scripts/read_data_KL15_XRF.R")
source("scripts/helper_functions.R")

set_1 <- get_color_palette()

dir <- paste0(getwd(), "/data/Africa_NE_200/data/")
data_kl15_xrf <- read.table(paste0(dir,'data_KL15_XRF.txt'), header = TRUE, sep = "\t")
data_kl15_xrf <- rename(data_kl15_xrf, depth = User_ID)
data_kl15_agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
data_kl15_agem <- rename(data_kl15_agem, age = best)
results <- read_data_kl15_xrf(data_kl15_xrf, data_kl15_agem)
data_kl15 <- results$data_kl15

#*******reproduction figure 6***********#
data_sel <- data[4:ncol(data_kl15)-1]
colnames(data_sel) <- gsub("_cts", "", colnames(data_sel))
x <- acomp(data_sel)
x_clr <- clr(data_sel)
# computing the mean of an acomp object returns the compositional mean (cf. Boogart et al. 2013, p.74)
x_m <- mean(x)
x_m_clr <- mean(clr(x))

df_row <- data.frame(name = names(x_m), value = as.vector(x_m), value_clr = as.vector(x_m_clr))

plot1 <- ggplot(df_row, aes(x = name, y = value)) +
  geom_bar(stat = "identity", fill = "#9ebcda") +
  labs(x = "Element", y = "Percentage") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
  axis.text.x = element_text(angle = 45, hjust = 1, , color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"))

plot3 <- ggplot(df_row, aes(x = name, y = value_clr)) +
  geom_bar(stat = "identity", fill = "#9ebcda") +
  labs(x = "Element", y = "clr-coefficients") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
  axis.text.x = element_text(angle = 45, hjust = 1, , color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"))


png("./scripts/figures/figure_6.png", width = 10, height = 5, units = "in", res = 300)
grid.arrange(plot1, plot3, ncol = 2, top = textGrob("Compositional mean on the simplex and in clr coefficients", gp = gpar(fontsize = 16)))
dev.off()

#*******reproduction figure 7***********#
pca_rob <- pcaCoDa(data_sel)

png("./scripts/figures/figure_7.png", width = 8, height = 8, units = "in", res = 300)
plot_pca_rotations(pca_rob$loadings, components = c(1,2), main = "PC1 vs. PC2", fixed=TRUE, pos_vector = c(4,4,4,2,4,4,2,2,4,2,4,2,2))
dev.off()

#*******reproduction biplots for appendix***********#
png("./scripts/figures/figure_D.png", width = 15, height = 10, units = "in", res = 300)

par(mfrow=c(2,3))

plot_pca_rotations(pca_rob$loadings, components = c(1,2), main = "PC1 vs. PC2", fixed=TRUE, pos_vector = c(4,4,4,2,4,4,2,2,4,2,4,2,2))
plot_pca_rotations(pca_rob$loadings, components = c(1,3), main = "PC1 vs. PC3", fixed=TRUE, pos_vector = c(4,4,4,2,4,4,2,2,4,2,4,2,2))
plot_pca_rotations(pca_rob$loadings, components = c(1,4), main = "PC1 vs. PC4", fixed=FALSE, pos_vector = c(4,4,4,2,4,4,2,2,4,2,4,2,2))
plot_pca_rotations(pca_rob$loadings, components = c(1,5), main = "PC1 vs. PC5", fixed=TRUE, pos_vector = c(4,4,4,2,4,4,2,2,4,2,4,2,2))
plot_pca_rotations(pca_rob$loadings, components = c(1,6), main = "PC1 vs. PC6", fixed=TRUE, pos_vector = c(4,4,4,2,4,4,2,2,4,2,4,2,2))
plot_pca_rotations(pca_rob$loadings, components = c(2,3), main = "PC2 vs. PC3", fixed=TRUE, pos_vector = c(4,4,4,2,4,4,2,2,4,2,4,2,2))

dev.off()


#*******reproduction figure 8***********#
pca_rotation_clr <- pca_rob$loadings

mean_comp_list <- lapply(1:6, function(i) {
  data.frame(
    component = 1:13,
    mean_value = 0,
    pc1_plus = pca_rotation_clr[,i] * sqrt(pca_rob$eigenvalues[i]),
    pc1_minus = -pca_rotation_clr[,i] * sqrt(pca_rob$eigenvalues[i]),
    names = colnames(data_sel)
  )
})

plot1 <- ggplot(mean_comp_list[[1]]) +
    geom_bar(aes(x = names, y = pc1_plus, color = "Factor 0.35"), 
             stat = "identity", fill = NA) +
    geom_bar(aes(x = names, y = mean_value), 
             stat = "identity", fill = NA, color = "black") +
    scale_color_manual(values = c("Factor 0.35" = "#9ebcda", "Factor -0.35" = "blue"),
                      name = "PC1",
                      breaks = c("Factor 0.35", "Factor -0.35")) +
    labs(x = "", y = "",
         title = "PC1 variations: 43%") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    ylim(-0.15,0.25)

plot2 <- ggplot(mean_comp_list[[2]]) +
    geom_bar(aes(x = names, y = pc1_plus, color = "Factor 0.27"), 
             stat = "identity", fill = NA) +
    geom_bar(aes(x = names, y = mean_value), 
             stat = "identity", fill = NA, color = "black") +
    scale_color_manual(values = c("Factor 0.27" = "#9ebcda", "Factor -0.23" = "blue"),
                      name = "PC2",
                      breaks = c("Factor 0.27", "Factor -0.23")) +
    labs(x = "", y = "",
         title = "PC2 variations: 26%") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")+
    ylim(-0.15,0.25)


plot3 <- ggplot(mean_comp_list[[3]]) +
    geom_bar(aes(x = names, y = pc1_plus, color = "Factor 0.18"), 
             stat = "identity", fill = NA) +
    geom_bar(aes(x = names, y = mean_value), 
             stat = "identity", fill = NA, color = "black") +
    scale_color_manual(values = c("Factor 0.18" = "#9ebcda", "Factor -0.11" = "blue"),
                      name = "PC3",
                      breaks = c("Factor 0.18", "Factor -0.11")) +
    labs(x = "", y = "",
         title = "PC3 variations: 12%") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    ylim(-0.15, 0.25)

plot4 <- ggplot(mean_comp_list[[4]]) +
    geom_bar(aes(x = names, y = pc1_plus, color = "Factor 0.15"), 
             stat = "identity", fill = NA) +
    geom_bar(aes(x = names, y = mean_value), 
             stat = "identity", fill = NA, color = "black") +
    scale_color_manual(values = c("Factor 0.15" = "#9ebcda", "Factor -0.09" = "blue"),
                      name = "PC4",
                      breaks = c("Factor 0.15", "Factor -0.09")) +
    labs(x = "", y = "",
         title = "PC4 variations: 8%") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") 

plot5 <- ggplot(mean_comp_list[[5]]) +
    geom_bar(aes(x = names, y = pc1_plus, color = "Factor 0.12"), 
             stat = "identity", fill = NA) +
    geom_bar(aes(x = names, y = mean_value), 
             stat = "identity", fill = NA, color = "black") +
    scale_color_manual(values = c("Factor 0.12" = "#9ebcda", "Factor -0.09" = "blue"),
                      name = "PC4",
                      breaks = c("Factor 0.12", "Factor -0.09")) +
    labs(x = "", y = "",
         title = "PC5 variations: 5%") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") 

plot6 <- ggplot(mean_comp_list[[6]]) +
    geom_bar(aes(x = names, y = pc1_plus, color = "Factor 0.09"), 
             stat = "identity", fill = NA) +
    geom_bar(aes(x = names, y = mean_value), 
             stat = "identity", fill = NA, color = "black") +
    scale_color_manual(values = c("Factor 0.09" = "#9ebcda", "Factor -0.09" = "blue"),
                      name = "PC4",
                      breaks = c("Factor 0.09", "Factor -0.09")) +
    labs(x = "", y = "",
         title = "PC6 variations: 3%") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") 

png("./scripts/figures/figure_8.png", width = 10, height = 5, units = "in", res = 300)
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 3)
dev.off()

