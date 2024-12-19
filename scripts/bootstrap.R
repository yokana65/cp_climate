bbootstrap_pca_ts <- function(data,
                              block_length = 50,
                              replicates = 1000,
                              reference_pca = NULL,
                              scale = 1,
                              eps = 0.01) {
  n <- nrow(data)
  results <- list()

  D <- ncol(data)
  basis_vectors <- lapply(1:(D - 1), generate_orthonormal_basis, D)
  basis_matrix <- do.call(rbind, basis_vectors)

  data <- data * scale

  if (is.null(reference_pca)) {
    reference_pca <- prcomp(clr(data), scale. = FALSE)
    # TODO: remove clr dependency
  }

  for(i in 1:replicates) {
    cat("Starts Replicate:", i, "\n")
    blocks <- ceiling(n/block_length)
    start_points <- sample(1:(n - block_length + 1), blocks, replace = TRUE)

    boot_indices <- unlist(lapply(start_points, function(x) x:(x + block_length - 1)))
    boot_indices <- boot_indices[boot_indices <= n][1:n]

    boot_data <- data[boot_indices, ]
    boot_results <- fit_pca_ilr_vs_4(boot_data,
                       max_iter = 50,
                       r = 10,
                       lambda = 1,
                       eps = eps,
                       sc_factor = 1,
                       sum_exp = TRUE)
    
    boot_pca <- boot_results$pca

    boot_pca$rotation <- t(basis_matrix) %*% boot_pca$rotation %*% basis_matrix

    for(j in seq_len(ncol(boot_pca$rotation))) {
      cor_sign <- sign(sum(boot_pca$rotation[, j] * reference_pca$rotation[, j]))
      if (cor_sign < 0) {
        boot_pca$rotation[, j] <- -boot_pca$rotation[, j]
        # boot_pca$x[, j] <- -boot_pca$x[, j]
      }
    }

    results[[i]] <- list(
      rotation = basis_matrix %*% boot_pca$rotation %*% t(basis_matrix),
      sdev = boot_pca$sdev
    )
  }

  return(results)
}

bbootstrap_pca_ts_par <- function(data,
                             block_length = 50,
                             replicates = 1000,
                             reference_pca = NULL,
                             scale = 1,
                             eps = 0.01,
                             workers = 6) {
    start_time <- Sys.time()
    n <- nrow(data)
    data <- data * scale

    D <- ncol(data)
    basis_vectors <- lapply(1:(D - 1), generate_orthonormal_basis, D)
    basis_matrix <- do.call(rbind, basis_vectors)

    if (is.null(reference_pca)) {
        reference_pca <- prcomp(clr(data), scale. = FALSE)
    }
    
    plan(multisession, workers = workers)

    results <- future_lapply(1:replicates, function(i) {
      cat("Starts Replicate:", i, "\n")
      blocks <- ceiling(n/block_length)
      start_points <- sample(1:(n - block_length + 1), blocks, replace = TRUE)
      
      boot_indices <- unlist(lapply(start_points, function(x) x:(x + block_length - 1)))
      boot_indices <- boot_indices[boot_indices <= n][1:n]
      
      boot_data <- data[boot_indices, ]
      boot_results <- fit_pca_ilr_vs_4(boot_data,
                                      max_iter = 50,
                                      r = 10,
                                      lambda = 1,
                                      eps = eps,
                                      sc_factor = 1,
                                      sum_exp = TRUE)
      boot_pca <- boot_results$pca
      
      for(j in seq_len(ncol(boot_pca$rotation))) {
          cor_sign <- sign(sum(boot_pca$rotation[, j] * reference_pca$rotation[, j]))
          if (cor_sign < 0) {
              boot_pca$rotation[, j] <- -boot_pca$rotation[, j]
              # boot_pca$x[, j] <- -boot_pca$x[, j]
          }
      }
      cat("Ends Replication:", i, "\n")      
      list(
          rotation = boot_pca$rotation,
          sdev = boot_pca$sdev
      )
    }, future.seed = TRUE)
    
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    print(paste("The bootstrap finished after:", elapsed_time, attr(elapsed_time, "units")))
    return(results)
}

bbootstrap_pca_ts_par_vs2 <- function(data,
                             block_length = 50,
                             replicates = 1000,
                             reference_pca = NULL,
                             scale = 1,
                             eps = 0.01,
                             workers = 6) {
    start_time <- Sys.time()
    n <- nrow(data)
    data <- data * scale

    D <- ncol(data)
    basis_vectors <- lapply(1:(D - 1), generate_orthonormal_basis, D)
    basis_matrix <- do.call(rbind, basis_vectors)

    if (is.null(reference_pca)) {
        reference_pca <- prcomp(clr(data), scale. = FALSE)
    }
    
    plan(multisession, workers = workers)

    results <- future_lapply(1:replicates, function(i) {
      cat("Starts Replicate:", i, "\n")
      blocks <- ceiling(n/block_length)
      start_points <- sample(1:(n - block_length + 1), blocks, replace = TRUE)
      
      boot_indices <- unlist(lapply(start_points, function(x) x:(x + block_length - 1)))
      boot_indices <- boot_indices[boot_indices <= n][1:n]
      
      boot_data <- data[boot_indices, ]
      boot_results <- fit_pca_ilr_vs_4(boot_data,
                                      max_iter = 50,
                                      r = 10,
                                      lambda = 1,
                                      eps = eps,
                                      sc_factor = 1,
                                      sum_exp = TRUE)
      boot_pca <- boot_results$pca
      
      for(j in seq_len(ncol(boot_pca$rotation))) {
        max_idx <- which.max(abs(boot_pca$rotation[, j]))
        cor_sign <- sign(boot_pca$rotation[max_idx, j] * reference_pca$rotation[max_idx, j])
          if (cor_sign < 0) {
              boot_pca$rotation[, j] <- -boot_pca$rotation[, j]
          }
      }
      cat("Ends Replication:", i, "\n")      
      list(
          rotation = boot_pca$rotation,
          sdev = boot_pca$sdev
      )
    }, future.seed = TRUE)
    
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    print(paste("The bootstrap finished after:", elapsed_time, attr(elapsed_time, "units")))
    return(results)
}


# test #############################
# x_data <- tar_read(data_kl15_comp)
# x_data <- tar_read(sim_comp_2_smi_nSim_80)
# x_matrix <- x_data[[1]]$x_data_matrix
# results <-
#   bbootstrap_pca_ts(x_matrix,
#                         block_length = 50,
#                         replicates = 10,
#                         scale = 1,
#                         eps = 0.04)
# length(results)
# results[[1]]$rotation
# # needs around 200 seconds (without using any workers)
# results_2 <-
#   bbootstrap_pca_ts_par(x_matrix,
#                         block_length = 50,
#                         replicates = 10,
#                         scale = 1,
#                         eps = 0.04)

# length(results_2)
# results_2[[1]]$rotation
# needs only 1/3 third (using 6 workers)