#*** MCEM functions ***#
stabilize_weights <- function(log_weights) {
  max_log_weight <- max(log_weights)

  stable_weights <- exp(log_weights - max_log_weight)

  stable_weights / sum(stable_weights)
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