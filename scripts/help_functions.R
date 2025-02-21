#*** MCEM functions ***#
stabilize_weights <- function(log_weights) {
  max_log_weight <- max(log_weights)

  stable_weights <- exp(log_weights - max_log_weight)

  stable_weights / sum(stable_weights)
}

finalize_pca <- function(pca, prepared_data, x_data, weights, start_time) {

  pca$scores <- sapply(seq_along(prepared_data$x_data), function(i){
    optim_result <- optim(rep(0, length = length(pca$sdev)), log_conditional_scores ,
                          x_data_i = prepared_data$x_data[[i]], pca = pca,
                          basis_matrix =  prepared_data$H,
                          control = list(fnscale = -1), method = "Nelder-Mead")
    as.vector(optim_result$par)
  })
  pca$rotation_ilr <- pca$rotation
  pca$rotation <- prepared_data$H %*% pca$rotation_ilr
  pca$center <- as.vector(prepared_data$H %*% pca$center)
  pca$eigenvalues <- pca$sdev^2
  pca$sdev <- pca$sdev / sum(pca$sdev)

  rownames(pca$rotation) <- names(x_data)
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

validate_dimensions <- function(x_data) {
  lengths <- unique(sapply(x_data, length))
  if (length(lengths) != 1) {
    stop("All observations must have the same number of components")
  }

  D <- lengths[1]

  return(D)
}

get_helmert <- function(D) {    
  V <- matrix(0, nrow = D, ncol = D - 1)
    for (i in 1:ncol(V)) {
      V[1:i, i] <- 1/i
      V[i + 1, i] <- (-1)
      V[, i] <- V[, i] * sqrt(i/(i + 1)) * (-1)
    }
  return(V)
}

#***** functions for reproduction of thesis figures*****#
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

clrInverse <- function(clr_coords) {
    exp(clr_coords) / sum(exp(clr_coords))
}

sample_from_density <- function(n, density_estimate) {
    sample(density_estimate$x, size = n, prob = density_estimate$y, replace = TRUE)
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

calculate_pc_scores <- function(data, pca, pc_number) {
  centered_data <- scale(data, center = TRUE, scale = FALSE)

  scores <- centered_data %*% pca$loadings[, pc_number]

  return(scores)
}