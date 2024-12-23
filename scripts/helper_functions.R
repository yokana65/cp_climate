stabilize_weights <- function(log_weights) {
  max_log_weight <- max(log_weights)

  stable_weights <- exp(log_weights - max_log_weight)

  stable_weights / sum(stable_weights)
}

monitor_weights <- function(weights, iteration) {
  for (i in seq_along(weights)) {
    ess <- 1/sum(weights[[i]]^2)
    max_weight <- max(weights[[i]])
    cat(sprintf("Iteration %d:\n
     ESS: %.2f\n 
     Max weight: %.2e\n 
     Observation: %d\n",
                iteration,
                ess,
                max_weight,
                i))
  }
}

monitor_global_ess <- function(all_weights, k) {
  mean_ess <- mean(sapply(all_weights, function(w) 1/sum(w^2)))
  cat(sprintf("Iteration %d: Mean ESS = %.2f\n", k, mean_ess))
}

generate_orthonormal_basis <- function(k, D) {
  scaling_factor <- -sqrt(k / (k + 1))

  basis_vector <- c(rep(1 / k, k), -1, rep(0, D - k - 1))

  basis_vector <- scaling_factor * basis_vector

  return(basis_vector)
}

plot_pca_rotation <- function(rotation, scale = 1, main = "PCA - clr") {
  plot(rotation[, 1], rotation[, 2],
       xlab = "PC1",
       ylab = "PC2",
       main = main)

  arrows(0, 0,
         rotation[, 1] * scale,
         rotation[, 2] * scale,
         length = 0.1, col = "blue")

  labels <- if (!is.null(rownames(rotation))) rownames(rotation)
  else seq_len(nrow(rotation))

  text(rotation[, 1], rotation[, 2],
       labels = labels,
       pos = 4,
       cex = 0.8)

  abline(h = 0, v = 0, lty = 2, col = "gray")
}

plot_pca_rotations <- function(rotation, components = c(1,2), scale = 1, main = "PCA - clr", fixed = FALSE) {
    # Set up the plot with optional fixed axes
    if (fixed) {
        plot(rotation[, components[1]], rotation[, components[2]],
             xlab = paste0("PC", components[1]),
             ylab = paste0("PC", components[2]),
             main = main,
             xlim = c(-0.6, 0.6),
             ylim = c(-0.6, 0.6))
    } else {
        plot(rotation[, components[1]], rotation[, components[2]],
             xlab = paste0("PC", components[1]),
             ylab = paste0("PC", components[2]),
             main = main)
    }

    arrows(0, 0,
           rotation[, components[1]] * scale,
           rotation[, components[2]] * scale,
           length = 0.1, col = "blue")
    
    labels <- if (!is.null(rownames(rotation))) rownames(rotation)
              else seq_len(nrow(rotation))
    
    text(rotation[, components[1]], rotation[, components[2]],
         labels = labels,
         pos = 4,
         cex = 0.8)
    
    abline(h = 0, v = 0, lty = 2, col = "gray")
}

plot_marginal_scores <- function(proposal_scores, weights, iteration) {
  # Extract dimensions
  n_obs <- length(proposal_scores)
  n_components <- nrow(proposal_scores[[1]])
  
  # Set up plotting layout
  par(mfrow = c(2, ceiling(n_components/2)))
  
  # For each component
  for(j in 1:n_components) {
    # Combine scores and weights from all observations
    all_scores <- unlist(lapply(1:n_obs, function(i) proposal_scores[[i]][j,]))
    all_weights <- unlist(lapply(1:n_obs, function(i) weights[[i]]))
    
    # Create weighted density plot
    density_est <- density(all_scores, weights = all_weights)
    plot(density_est, main = paste("Component", j, "- Iteration", iteration),
         xlab = "Score", ylab = "Density")
    
    # Add rug plot of individual scores
    rug(all_scores)
  }
  par(mfrow = c(1,1))
}

clrInv_long <- function(clr_coords) {
    # Exponentiate the clr coordinates
    exp_coords <- exp(clr_coords)
    
    # Calculate the geometric mean normalization constant
    norm_const <- sum(exp_coords)
    
    # Return normalized compositions
    exp_coords / norm_const
}

sample_from_density <- function(n, density_estimate) {
    sample(density_estimate$x, size = n, prob = density_estimate$y, replace = TRUE)
}

seq_ilr_transform <- function(x_data) {
  # sequential calculation of ilr Coordinates (Cp. robCompositions::pcaCoDa)
  x_ilr = matrix(NA, nrow = nrow(x_data), ncol = ncol(x_data) - 1)
  for (i in seq_len(x_ilr)) {
    x_ilr[, i] <- sqrt((i) / (i + 1)) *
      log(((apply(as.matrix(x_data[, 1:i]), 1, prod))^(1/i)) / (x_data[, i + 1]))
  }
  return(x_ilr)
}
# TODO: check if this results in same coorddinates as basis_matrix multiplication

get_helmert <- function(x) {    
  V <- matrix(0, nrow = ncol(x), ncol = ncol(x) - 1)
    for (i in 1:ncol(V)) {
      V[1:i, i] <- 1/i
      V[i + 1, i] <- (-1)
      V[, i] <- V[, i] * sqrt(i/(i + 1))
    }
  return(V)
}
# TODO: make connection between Helmert Matrix and Bases in theoretical part
fix_sign <- function(A) {
    mysign <- function(x) ifelse(x < 0, -1L, 1L)
    A[] <- apply(A, 2L, function(x) x * mysign(x[1L]))
    A
}