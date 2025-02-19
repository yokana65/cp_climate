gradient_lcs <- function(scores,
                              x_data_i,
                              pca,
                              basis_matrix) {
  m_i <- sum(x_data_i)
  clr_comp <- as.vector(
                        basis_matrix %*% pca$center +
                          basis_matrix %*%  pca$rotation %*% scores)


  grad_vecs <- sapply(seq_along(scores), function(k) {
    e_k <- basis_matrix[, k]
    v_k <- pca$rotation[, k]
    term1 <- sum(x_data_i * e_k)
    term2 <- m_i * sum(clr_comp * e_k)

    grad_k <- v_k * (term1 - term2)

    return(grad_k)
  })
  grad <- colSums(grad_vecs)

  grad <- grad - scores / (pca$sdev^2)

  return(grad)
}
