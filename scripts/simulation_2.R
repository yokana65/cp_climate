build_setting_2comp_5parts_vs2 <- function(n_observations = 100,
                            eigenvalues = c(0.6, 0.3, 0.05, 0.05),
                            mean = c(0, 1, 0.5, -1, -0.5),
                            n_counts = 500) {
  # set.seed(123)
  # follow procedure outlined by Steyer and Greven
  # with five dimensions the last dimension is a constant 
  pc_1 <- c(0, 1, -1, 0, 0)  # Contrast between parts 2 and 3 in clr space
  pc_2 <- c(-1, 0, 0, 1, 0)  # Contrast between parts 1 and 4 in clr space
  # the problem is that those clr-coordinates are not calculated with the same basis

  lambda_1 <- eigenvalues[1]
  lambda_2 <- eigenvalues[2]

  
  clr_coords <- lapply(1:n_observations, function(i){
    mean + rnorm(1, 0, sqrt(lambda_1))*pc_1 + rnorm(1, 0, sqrt(lambda_2))*pc_2
  })

  composition_list <- lapply(clr_coords, clrInv)

  x_data <- lapply(1:n_observations, function(i) {
    probs <- composition_list[[i]]
    rmultinom(1, n_counts, probs)[, 1]
  })
  x_data_matrix <- do.call(rbind, x_data)

  return(list("x_data" = x_data,"x_data_matrix" = x_data_matrix))
}
