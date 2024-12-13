build_setting_4comp_13parts_vs2 <- function(data,
    eigenvalues = c(0.5, 0.2, 0.12, 0.08),
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532),
    scale = 0.01, 
    n_observations = nrow(data)) {
  # simulate the sample size
  data <- data * scale

  density_estimate <- density(data$aggregate)

  n_counts <- sample_from_density(n_observations, density_estimate)

  pc_1 <- c(0.7, 0, 0, -0.2, 0, 0, -0.2, -0.2, 0.4, -0.1, 0, -0.2, -0.2) 
  pc_2 <- c(0.5, 0, -0.4, 0.4, 0, 0, 0, 0, 0, 0, -0.5, 0, 0)
  pc_3 <- c(0, 0, 0.3, 0, 0.4, 0, 0, 0, -0.7, 0, 0, 0, 0)
  pc_4 <- c(0.3, -0.2, 0, -0.3, 0.6, 0.3, 0, 0, -0.4, 0, 0, 0, -0.3)

  normalize <- function(v) {
      v / sqrt(sum(v^2))
  }
  
  # Normalize each PC vector
  pc_1_norm <- normalize(pc_1)
  pc_2_norm <- normalize(pc_2)
  pc_3_norm <- normalize(pc_3)
  pc_4_norm <- normalize(pc_4)

  lambda_1 <- eigenvalues[1]
  lambda_2 <- eigenvalues[2]
  lambda_3 <- eigenvalues[3]
  lambda_4 <- eigenvalues[4]
  
  clr_coords <- lapply(1:n_observations, function(i){
    mean + rnorm(1, 0, sqrt(lambda_1)) * pc_1_norm +
      rnorm(1, 0, sqrt(lambda_2)) * pc_2_norm +
      rnorm(1, 0, sqrt(lambda_3)) * pc_3_norm +
      rnorm(1, 0, sqrt(lambda_4)) * pc_4_norm
  })

  composition_list <- lapply(clr_coords, clrInv)

  x_data <- lapply(1:n_observations, function(i) {
    probs <- composition_list[[i]]
    rmultinom(1, n_counts[i] , probs)[, 1]
  })
  x_data_matrix <- do.call(rbind, x_data)

  return(list("x_data" = x_data,"x_data_matrix" = x_data_matrix))
}
