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

get_helmert <- function(D) {    
  V <- matrix(0, nrow = D, ncol = D - 1)
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

clr_transform <- function(x, replace_zeros = "fraction") {
  replace_zeros <- match.arg(replace_zeros, c("neutral", "fraction"))
  if (!is.vector(x) && !is.matrix(x)) {
    stop("Input must be a vector or a matrix")
  }
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  if (any(is.na(x))) {
    stop("Input contains NA values")
  }
  if (any(!is.finite(x))) {
    stop("Input contains infinite values")
  }
  if (any(x < 0)) {
    stop("All values must be positive for log-ratio transformation")
  }
  
  # Behandlung von Nullwerten
  if (replace_zeros == "fraction") {
    x[x == 0] <- 0.0333
  }
  log_x <- log(x)
  # Berechne CLR für jede Zeile
  clr_values <- t(apply(log_x, 1, function(row) {
    if (replace_zeros == "neutral") {
      row - mean(row[row != -Inf])
    } else {
      row - mean(row)
    }
  }))
  if (replace_zeros == "neutral") {
    clr_values[x == 0] <- 0
  }
  # Rückgabe als Vektor wenn Eingabe ein Vektor war
  if (nrow(x) == 1 && is.vector(x)) {
    clr_values <- as.vector(clr_values)
  }
  
  return(clr_values)
}

clrTransVec <- function(x) {
  if (!is.vector(x)) {
    stop("Input must be a vector")
  }  
  if (any(is.na(x))) {
    stop("Input contains NA values")
  }
  if (any(!is.finite(x))) {
    stop("Input contains infinite values")
  }  
  if (any(x <= 0)) {
    stop("All values must be positive for log-ratio transformation")
  }
  log_x <- log(x)
  
  mean_log_x <- mean(log_x)

  clr_values <- log_x - mean_log_x
  
  return(clr_values)
}

inv_clr <- function(y, check_sum = FALSE) {
  if (!is.vector(x) && !is.matrix(x)) {
    stop("Input must be a vector or a matrix")
  }
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  if (check_sum) {
    row_sums <- rowSums(y)
    if (any(abs(row_sums) > 1e-10)) {
      stop("The sum of CLR coordinates must be close to zero for each row.")
    }
  }
  exp_values <- exp(y)
  
  # Normalisiere die Werte, damit sie sich zu 1 summieren
  normalized_values <- exp_values / rowSums(exp_values)
  
  return(normalized_values)
}

ilr_transform <- function(x) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  if (any(x < 0)) {
    stop("Alle Werte müssen positiv sein, um die ILR-Transformation anzuwenden.")
  }
  V <- get_helmert(ncol(x))

  clr_x <- clr_transform(x)

  # Anwendung der ILR-Transformation
  ilr_x <- clr_x %*% V
  
  return(ilr_x)
}

inv_ilr <- function(ilr_x, use_transpose = FALSE) {
  # Überprüfen, ob ilr_x eine Matrix oder eine Liste ist
  if (!is.matrix(ilr_x) && !is.list(ilr_x)) {
    stop("Input must be a matrix or a list")
  }
  
  # Anzahl der Komponenten bestimmen
  if (is.matrix(ilr_x)) {
    D <- ncol(ilr_x) + 1
  } else { # ilr_x ist eine Liste
    D <- length(ilr_x[[1]]) + 1
    ilr_x <- do.call(rbind, ilr_x) # Konvertiere Liste in Matrix
  }
  # Erstellung der Helmert-Matrix
  V <- get_helmert(D)
  
  # Rücktransformation in den CLR-Raum
  if (use_transpose) {
    clr_x <- ilr_x %*% t(V)
  } else {
    # Berechnung der inversen Matrix von V nur wenn benötigt
    V_inv <- MASS::ginv(V)
    clr_x <- ilr_x %*% V_inv
  }
  # Rücktransformation in den kompositionellen Raum
  x <- exp(clr_x)
  x <- x / rowSums(x)
  
  return(x)
}

replace_zeros <- function(x, repl = 0.0333) {
  # Überprüfen, ob x ein Vektor, eine Matrix oder ein Data Frame ist
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  } else if (!is.vector(x) && !is.matrix(x)) {
    stop("Die Eingabe muss ein Vektor, eine Matrix oder ein Data Frame sein.")
  }
  
  # Ersetzen der Nullwerte durch 0.0333
  x[x == 0] <- repl
  
  return(x)
}

egv_clr2ilr <- function(egv_clr, E) {
  # Get number of relevant eigenvectors (D-1)
  n_eigenvectors <- ncol(E)
  
  # Transform each eigenvector separately
  transformed_vectors <- lapply(1:n_eigenvectors, function(i) {
    t(egv_clr[,i] %*% E)
  })
  
  # Combine transformed vectors
  result <- do.call(cbind, transformed_vectors)
  
  return(result)
}