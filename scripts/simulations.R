build_setting_2comp_5parts <- function(n_observations = 100,
                            eigenvalues = c(0.6, 0.3, 0.05, 0.05),
                            mean = c(0, 2, 0.5, -2, -0.5),
                            n_counts = 500) {
  # set.seed(123)
  v1 <- c(0, 1/sqrt(2), -1/sqrt(2), 0, 0)  # Contrast between parts 2 and 3
  v2 <- c(-1/3, -1/3, -1/3, 1, 0)  # Focus on part 4
  v3 <- c(1/sqrt(2), 0, 0, -1/sqrt(2), 0)  # Additional contrast
  v4 <- c(1/2, 1/2, 0, 0, -1)  # Additional contrast

  V <- cbind(v1, v2, v3, v4)
  Sigma <- V %*% diag(eigenvalues) %*% t(V)
  
  clr_coords <- rmvnorm(n_observations, mean = mean, sigma = Sigma)

  compositions <- clrInv(clr_coords)
  composition_list <- apply(compositions, 1, function(x) x, simplify = FALSE)

  x_data <- lapply(1:n_observations, function(i) {
    probs <- composition_list[[i]]
    rmultinom(1, n_counts, probs)[, 1]
  })
  x_data_matrix <- do.call(rbind, x_data)

  return(list("x_data" = x_data, "Sigma" = Sigma,"x_data_matrix" = x_data_matrix))
}