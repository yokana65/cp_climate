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


complex_setting <- function(data,
  eigenvalues = c(0.43, 0.26, 0.12, 0.08),
  mean = c(-1.21, -1.76, 1.70, -0.6,
           1.04, -3.62, -1.26, 1.04, -0.83,
           0.50, 4.0, -0.49, 1.51),
  n_observations = nrow(data)) {

  mean_data <- mean(data$aggregate)
  sd_data <- sd(data$aggregate)

  n_counts <- round(rnorm(n_observations, mean = mean_data, sd = sd_data))

  pc_1 <- c(0.4, 0, 0.3, -0.3, 0, 0, -0.2, -0.2, 0.3, -0.1, 0.3, -0.3, -0.2)
  pc_2 <- c(0.7, 0, -0.3, 0.2, 0, -0.2, 0, 0, 0, 0, -0.4, 0, 0)
  pc_3 <- c(0.2, 0, -0.3, 0, -0.4, 0.6, -0.3, 0.2, 0, 0, 0, 0, 0)
  pc_4 <- c(0.3, 0, 0, 0, 0, 0.2, 0, 0, -0.8, 0, 0.3, 0, 0)

  normalize <- function(v) {
    v / sqrt(sum(v^2))
  }

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

complex_setting_autoc <- function(data,
                          eigenvalues = c(0.43, 0.26, 0.12, 0.08),
                          mean = c(-1.21, -1.76, 1.70, -0.6,
                                 1.04, -3.62, -1.26, 1.04, -0.83,
                                 0.50, 4.0, -0.49, 1.51),
                          n_observations = nrow(data),
                          rho = 0.7,
                          lag = 10) {
pc_matrix <- matrix(c(
    c(0.4, 0, 0.3, -0.3, 0, 0, -0.2, -0.2, 0.3, -0.1, 0.3, -0.3, -0.2),
    c(0.7, 0, -0.3, 0.2, 0, -0.2, 0, 0, 0, 0, -0.4, 0, 0),
    c(0.2, 0, -0.3, 0, -0.4, 0.6, -0.3, 0.2, 0, 0, 0, 0, 0),
    c(0.3, 0, 0, 0, 0, 0.2, 0, 0, -0.8, 0, 0.3, 0, 0)
), ncol = 4)

pc_norm <- apply(pc_matrix, 2, function(v) v / sqrt(sum(v^2)))

library(stats)

ar_coefs <- c(0.3, 0.1, -0.05, 0.1, 0.05, 0.04, -0.02, 0.02, 0.01, -0.01)

coefs <- lapply(eigenvalues, function(ev) {
    arima.sim(
      model = list(ar = ar_coefs),
      n = n_observations,
      sd = sqrt(eigenvalues)
    )
})

clr_coords <- lapply(1:n_observations, function(i) {
    mean + Reduce('+', Map(function(coef, pc) coef[i] * pc,
                          coefs, split(pc_norm, col(pc_norm))))
})

composition_list <- lapply(clr_coords, clrInv)

mean_data <- mean(data$aggregate)
sd_data <- sd(data$aggregate)
n_counts <- round(rnorm(n_observations, mean = mean_data, sd = sd_data))

x_data <- lapply(1:n_observations, function(i) {
    probs <- composition_list[[i]]
    rmultinom(1, n_counts[i], probs)[, 1]
})

x_data_matrix <- do.call(rbind, x_data)

return(list("x_data" = x_data, "x_data_matrix" = x_data_matrix))
}

generate_ar<- function(n, rho, sigma, lag_max = 10) {
    eps <- rnorm(n, 0, sigma)
    x <- numeric(n)
    x[1:lag_max] <- eps[1:lag_max]
    
    for(t in (lag_max+1):n) {
        x[t] <- sum(rho^(1:lag_max) * rev(x[(t-lag_max):(t-1)])) + eps[t]
    }
    return(x)
}

complex_setting_auto <- function(data,
                          eigenvalues = c(0.43, 0.26, 0.12, 0.08),
                          mean = c(-1.21, -1.76, 1.70, -0.6,
                                 1.04, -3.62, -1.26, 1.04, -0.83,
                                 0.50, -4.0, -0.49, 1.51),
                          n_observations = nrow(data),
                          n_periods = 10) {
  
  mean_data <- mean(data$aggregate)
  sd_data <- sd(data$aggregate)
  n_counts <- round(rnorm(n_observations, mean = mean_data, sd = sd_data))
  
  pc_1 <- c(0.4, 0, 0.3, -0.3, 0, 0, -0.2, -0.2, 0.3, -0.1, 0.3, -0.3, -0.2)
  pc_2 <- c(0.7, 0, -0.3, 0.2, 0, -0.2, 0, 0, 0, 0, -0.4, 0, 0)
  pc_3 <- c(0.2, 0, -0.3, 0, -0.4, 0.6, -0.3, 0.2, 0, 0, 0, 0, 0)
  pc_4 <- c(0.3, 0, 0, 0, 0, 0.2, 0, 0, -0.8, 0, 0.3, 0, 0)

  
  pc_1_norm <- normalize(pc_1)
  pc_2_norm <- normalize(pc_2)
  pc_3_norm <- normalize(pc_3)
  pc_4_norm <- normalize(pc_4)



  # coef1 <- generate_sin_coef(n_observations, n_periods, sqrt(eigenvalues[1]))
  # coef2 <- generate_cos_coef(n_observations, n_periods, sqrt(eigenvalues[2]))
  # coef3 <- generate_sin_coef(n_observations, n_periods, sqrt(eigenvalues[3]))
  # coef4 <- generate_cos_coef(n_observations, n_periods, sqrt(eigenvalues[4]))

  coef1 <- generate_temporal_scores(n_observations, n_periods, 
                                    eigenvalues[1], phase = 0)
  coef2 <- generate_temporal_scores(n_observations, n_periods, 
                                    eigenvalues[2], phase = pi / 4)
  coef3 <- generate_temporal_scores(n_observations, n_periods, 
                                    eigenvalues[3], phase = pi / 2)
  coef4 <- generate_temporal_scores(n_observations, n_periods, 
                                    eigenvalues[4], phase = pi / 3)

  clr_coords <- lapply(1:n_observations, function(i){
    mean + coef1[i] * pc_1_norm + coef2[i]*pc_2_norm +
      coef3[i] * pc_3_norm + coef4[i] * pc_4_norm
  })

  composition_list <- lapply(clr_coords, clrInv)

  x_data <- lapply(1:n_observations, function(i) {
      probs <- composition_list[[i]]
      rmultinom(1, n_counts[i], probs)[, 1]
  })

  x_data_matrix <- do.call(rbind, x_data)

  return(list("x_data" = x_data, "x_data_matrix" = x_data_matrix))
}

generate_sin_coef <- function(n, periods, sigma) {
    t <- seq(0, 2*pi*periods, length.out = n)
    sin_wave <- sin(t)
    return(sin_wave * sigma + rnorm(n, 0, sigma))
}

generate_cos_coef <- function(n, periods, sigma) {
    t <- seq(0, 2*pi*periods, length.out = n)
    sin_wave <- cos(t)
    return(sin_wave * sigma + rnorm(n, 0, sigma))
}

generate_temporal_scores <- function(n_observations, n_periods, eigenvalue, phase = 0) {
  base_scores <- rnorm(n_observations, mean = 0, sd = sqrt(eigenvalue))
  t <- seq(0, 2*pi*n_periods, length.out = n_observations)
  temporal_component <- sin(t + phase)
  scores <- base_scores + temporal_component * sqrt(eigenvalue)
  return(scores)
}

# plot(coef1, type = "l", col = "blue")
# plot(coef1, type = "l", col = "blue", ylim = range(c(coef1, coef2, coef3, coef4)))
# lines(coef2, col = "red")
# lines(coef3, col = "green")
# lines(coef4, col = "purple")
# legend("topright", legend = c("PC1", "PC2", "PC3", "PC4"), 
#        col = c("blue", "red", "green", "purple"), lty = 1)