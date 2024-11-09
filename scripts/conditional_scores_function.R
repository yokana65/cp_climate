conditional_scores_log_ilr <- function(scores,
                                          x_data_i,
                                          pca,
                                          basis_matrix,
                                          sc_factor) {
  scaling_factor <- sc_factor
  ilr_comp <- as.vector(pca$center + pca$rotation %*% scores)
  clr_comp <- ilr2clr(ilr_comp)
  norm_constant <- sum(exp(clr_comp))

  # Compute scaled log likelihood
  log_likelihood <- sum(x_data_i * clr_comp) - sum(x_data_i) * log(norm_constant)
  # apply scaling when sc_factor ungleich 1
  if (scaling_factor != 1) {
    log_likelihood <- log_likelihood - scaling_factor
  }

  # Prior remains unchanged as it's already well-scaled
  log_prior <- - sum(0.5 * scores^2 / (pca$sdev^2))

  return(log_likelihood + log_prior)
}