validate_dimensions <- function(x_data) {
  lengths <- unique(sapply(x_data, length))
  if (length(lengths) != 1) {
    stop("All observations must have the same number of components")
  }

  D <- lengths[1]

  return(D)
}

prepare_data <- function(x_data) {
  if (!is.list(x_data) && !is.matrix(x_data)) {
    stop("Input x_data must be a list or a matrix")
  }

  if (is.data.frame(x_data) || is.matrix(x_data)) {
    x_data <- apply(x_data, 1, function(x) x, simplify = FALSE)
  }

  D <- validate_dimensions(x_data)    
  H <- get_helmert(D)

  list(
    x_data = x_data,
    H = H,
    D = D
  )
}

finalize_pca <- function(pca, prepared_data, x_data, weights, start_time) {

  pca$rotation_ilr <- pca$rotation
  pca$rotation <- prepared_data$H %*% pca$rotation_ilr
  pca$center <- as.vector(prepared_data$H %*% pca$center)
  pca$eigenvalues <- pca$sdev^2
  pca$sdev <- pca$sdev / sum(pca$sdev)

  rownames(pca$rotation) <- names(x_data[[1L]])
  cn <- paste0("Comp.", seq_len(ncol(pca$rotation)))
  colnames(pca$rotation) <- cn

  ess <- sapply(weights, function(w) 1 / sum(w^2))
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time

  return(list(
    pca          = pca,
    ess          = ess,
    elapsed_time = elapsed_time
  ))
}


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


plot_pca_rotations <- function(rotation, components = c(1,2), scale = 1, main = "PCA - clr", fixed = FALSE, pos_vector = NULL) {
    if (fixed) {
        plot(rotation[, components[1]], rotation[, components[2]],
             xlab = paste0("PC", components[1]),
             ylab = paste0("PC", components[2]),
             main = main,
             xlim = c(-0.7, 0.7),
             ylim = c(-0.7, 0.7))
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
    
    default_pos_vector <- ifelse(rotation[, components[1]] < 0, 2, 4)

    if (is.null(pos_vector)) {
      pos_vector <- default_pos_vector
    } else {
      if (length(pos_vector) != length(labels)) {
        warning("Length of 'pos_vector' does not match the number of labels. Using default positions for all labels.")
        pos_vector <- default_pos_vector
      }
    }

    text(rotation[, components[1]], rotation[, components[2]],
         labels = labels,
         pos = pos_vector,
         cex = 0.8)
    
    abline(h = 0, v = 0, lty = 2, col = "gray")
}

calculate_pc_scores <- function(data, pca, pc_number) {
  centered_data <- scale(data, center = TRUE, scale = FALSE)

  scores <- centered_data %*% pca$loadings[, pc_number]

  return(scores)
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

clrInverse <- function(clr_coords) {
    exp(clr_coords) / sum(exp(clr_coords))
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
      V[, i] <- V[, i] * sqrt(i/(i + 1)) * (-1)
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
  
  if (replace_zeros == "fraction") {
    x[x == 0] <- 0.65 * 0.5/1
  }
  log_x <- log(x)

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

summary_stats <- function(x_data) {
  data.frame(
    mean = colMeans(x_data),
    sd = apply(x_data, 2, sd),
    min = apply(x_data, 2, min),
    max = apply(x_data, 2, max),
    n_zeros_components = apply(x_data, 2, function(x) sum(x == 0)),
    n_zeros_composition = sum(apply(x_data, 1, function(x) any(x == 0)))
  )
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

predict_latent_compositions_clr <- function(comp_pca, basis_matrix, sc_factor){
  predicted_scores <- sapply(seq_along(comp_pca$x_data), function(i){
    optim_result <- optim(rep(0, length = length(comp_pca$pca$sdev)), conditional_scores_log_ilr , gr = gradient_cslc_vs1,
                          x_data_i = comp_pca$x_data[[i]], pca = comp_pca$pca,
                          basis_matrix = basis_matrix,
                          sc_factor = sc_factor,
                          control = list(fnscale = -1), method = "BFGS")
    as.vector(optim_result$par)
  })
  
  compositions_list <- vector("list", ncol(predicted_scores))

  for(i in 1:ncol(predicted_scores)) {
      # Extract current score vector
      current_scores <- predicted_scores[,i]

      # Transform to CLR space
      ilr_center <- as.vector(comp_pca$pca$center)
      clr_center <- ilr_center %*% basis_matrix
      clr_rotation <- t(basis_matrix) %*% comp_pca$pca$rotation
      clr_composition <- clr_center + as.vector(clr_rotation %*% current_scores)

      # Transform to composition
      compositions_list[[i]] <- clrInv(clr_composition)
  }
  return(list("predicted_scores" = predicted_scores, "compositions" = compositions_list))
}

predict_coordinates <- function(comp_pca, x_df, basis_matrix){
  predicted_scores <- sapply(seq_along(comp_pca$x_data), function(i){
    optim_result <- optim(rep(0, length = length(comp_pca$pca$sdev)), conditional_scores_log_ilr_vs3b , gr = gradient_cslc_vs1c,
                          x_data_i = x_df[i, ], pca = comp_pca$pca,
                          basis_matrix = basis_matrix,
                          control = list(fnscale = -1), method = "BFGS")
    as.vector(optim_result$par)
  })
  
  coordinates_list <- vector("list", ncol(predicted_scores))

  for(i in 1:ncol(predicted_scores)) {
      current_scores <- predicted_scores[,i]

      coordinates_list[[i]] <- comp_pca$pca$center + as.vector(comp_pca$pca$rotation %*% current_scores)
  }
  return(list("predicted_scores" = predicted_scores, "coordinates" = coordinates_list))
}

load_required_packages <- function() {
  required_packages <- c("ggplot2", "gridExtra", "compositions",
                        "robCompositions", "zCompositions", "targets", "dplyr")
  
  new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  invisible(lapply(required_packages, library, character.only = TRUE))
}

get_color_palette <- function() {
  c("#e0ecf4", "#9ebcda", "#8856a7")
}

process_simulation_data <- function(simulation_data) {
  pca_results <- list(
    clr_rob = list(),
    clr_std = list()
  )
  
  for(i in 1:length(simulation_data)) {
    x_data <- do.call(rbind, simulation_data[[i]]$x_data)
    x_data_repl <- if(any(x_data == 0)) {
      cmultRepl(x_data, method = "CZM", suppress.print=TRUE)
    } else {
      x_data
    }
    
    pca_results$clr_std[[i]] <- prcomp(clr(x_data_repl))
    pca_results$clr_rob[[i]] <- pcaCoDa(x_data_repl, method = "robust")
  }
  
  return(pca_results)
}

normalize <- function(v) {
    v / sqrt(sum(v^2))
}