
log_conditional_scores <- function(scores,
                                   x_data_i,
                                   pca,
                                   basis_matrix) {
  clr_coef <- as.vector(
                        basis_matrix %*% pca$center +
                          basis_matrix %*% pca$rotation %*% scores)
  norm_constant <- sum(exp(clr_coef))

  log_likelihood <- sum(
                        x_data_i * clr_coef) -
                          sum(x_data_i) * log(norm_constant)

  log_prior <- - sum(0.5 * scores^2 / (pca$sdev^2))

  return(log_likelihood + log_prior)
}
