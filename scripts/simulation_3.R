build_setting_4comp_5parts_vs2 <- function(n_observations = 100,
                            eigenvalues = c(0.5, 0.3, 0.1, 0.1),
                            mean = c(0, 0.9, 0.3, -0.8, -0.2),
                            n_counts = 500) {
  # set.seed(123)
  # follow a more complex structure with contrary influences -> a more realistic example
  # with five dimensions the last dimension is a constant 
  pc_1 <- c(0, 1/sqrt(2), -1/sqrt(2), 0, 0)  # Contrast between parts 2 and 3
  pc_2 <-  c(-1/3, -1/3, -1/3, 1, 0)  # Contrast between parts 1 and 4
  pc_3 <- c(1/2, 1/2, 0, 0, -1)
  pc_4 <- c(1/sqrt(2), 0, 0, -1/sqrt(2), 0)

  lambda_1 <- eigenvalues[1]
  lambda_2 <- eigenvalues[2]
  lambda_3 <- eigenvalues[3]
  lambda_4 <- eigenvalues[4]
  
  clr_coords <- lapply(1:n_observations, function(i){
    mean + rnorm(1, 0, sqrt(lambda_1)) * pc_1 +
      rnorm(1, 0, sqrt(lambda_2)) * pc_2 +
      rnorm(1, 0, sqrt(lambda_3)) * pc_3 +
      rnorm(1, 0, sqrt(lambda_4)) * pc_4
  })

  composition_list <- lapply(clr_coords, clrInv)

  x_data <- lapply(1:n_observations, function(i) {
    probs <- composition_list[[i]]
    rmultinom(1, n_counts, probs)[, 1]
  })
  x_data_matrix <- do.call(rbind, x_data)

  return(list("x_data" = x_data,"x_data_matrix" = x_data_matrix))
}
