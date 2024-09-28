### rough example to demonstrate a balancing effect
# with the orthonormal basis being the principal components

# define number of components
n_components <- 13
x_grid <- seq(1, n_components)
# make the neutral element the center for a start
raw_data <- data.frame("x" = x_grid, "y" = rep(0, n_components))

center_function_comp <- function(comp_data){
  mean <- mean(comp_data[,2])
  comp_data[,2] <- comp_data[,2] - mean
  comp_data
}

clr_mean <- center_function_comp(raw_data)

# zwei Hauptkomponenten im clr-Raum
# PC1 is a balancing effect between the the first five elements and the sixths element
# alternative: use Egozcue et al. (2003) orthogonal basis
pc_1 <- data.frame("x" = x_grid, "y" = c(rep(sqrt(5/(5+1)*5^(-1)), 5),-1, rep(0, 7)))
pc_1 <- center_function_comp(pc_1)
# centering resolves in the sum of all elements being zero
# Normalisierung
pc_1[,2] <- pc_1[,2]/norm(pc_1[,2],type="2")


pc_2 <- data.frame("x" = x_grid, "y" = c(rep(sqrt(12/(12+1)*12^(-1)), 12),0))
pc_2 <- center_function_comp(pc_2)
pc_2[,2] <- pc_2[,2]/norm(pc_2[,2], type="2")

lambda_1 <- 0.5
lambda_2 <- 0.2

# Schritt 2: Simulieren Sie die Scores und konstruieren Sie die Dichten
n_data <- 30
true_observed_clr_comp <- sapply(1:n_data, function(i){
  clr_mean[,2] + rnorm(1, 0, lambda_1)*pc_1[,2] + rnorm(1, 0, lambda_2)*pc_2[,2]
})

# error checking: Sum to zero
check_columns_sum_to(true_observed_clr_densities, 0)

true_observed_comp <- lapply(1:n_data, function(i){
  # TODO: remove x_grid
  clr_density <- data.frame(x_grid, true_observed_clr_comp[,i])
  inverse_clr_trafo(clr_density)
})

# Schritt 3: Datenpunkte aus den Dichten ziehen
# Define number of counts when dealing with multinomial distributions
n_counts <- 2000
n_samples <- 40
# change the structure to a list of compositions
# x_data <- lapply(1:n_data, function(i){
#   probs <- true_observed_densities[[i]][,2]
#   t(rmultinom(n_samples, n_counts, probs))
# })
x_data <- unlist(lapply(1:n_data, function(i) {
  probs <- true_observed_comp[[i]][,2]
  samples <- t(rmultinom(n_samples, n_counts, probs))
  lapply(1:nrow(samples), function(j) samples[j,])
}), recursive = FALSE)

# turn x_data into a matrix
# TODO: Why do I need a list at all? 
x_data_matrix <- do.call(rbind, x_data)

library(compositions)

x_acomp <- acomp(x_data_matrix)
x_clr <- clr(x_acomp)
# dies ist equivalent zur Summe über k bis D-1 con clr(x_acomp[,k]*e_k))
# TODO

# # Schritt 4: Dichte-PCA schätzen

pcx <- princomp(x_acomp)
plot(pcx)
biplot(pcx)

# compare loadings
pc_1_inv <- inverse_clr_trafo(pc_1)
                                                                                                                                                                                                                                                            
# density_pca <- fit_density_pca(x_data, max_iter = 10)

# # Schritt 5: Ergebnisse plotten
# plot_pca(density_pca$pca, x_grid = density_pca$x_grid)

# # Schritt 6: Synthetische Kompositionen konstruieren
# predicted <- predict_latent_densities(density_pca)
# predicted_clr_densities <- predicted$clr_densities
# predicted_densities <- lapply(1:ncol(predicted_clr_densities), function(i){
#   inverse_clr_trafo(data.frame(x_grid, predicted_clr_densities[,i]))
# })

# # Plotten Sie einige der synthetischen Kompositionen
# par(mfrow = c(2, 2))
# for (i in 1:4) {
#   barplot(predicted_densities[[i]][,2], names.arg = x_grid, main = paste("Synthetic Composition", i))
# }

check_columns_sum_to <- function(data, integer) {
  n_cols <- ncol(data)
  
  columns_sum_to_zero <- logical(n_cols)
  
  for (i in 1:n_cols) {
    column_sum <- sum(data[, i])
    columns_sum_to_zero[i] <- (column_sum == integer)
  }
  
  return(columns_sum_to_zero)
}

inverse_clr_trafo <- function(clr_density){
  f_integral <- sum(exp(clr_density[,2]))
  data.frame("x" = clr_density[,1], "y" = exp(clr_density[,2])/f_integral)
}

scores <- rnorm(n_components, 0, 1)
