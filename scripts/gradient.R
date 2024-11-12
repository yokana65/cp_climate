gradient_cslc_vs1 <- function(scores,
                                 x_data_i,
                                 pca,
                                 basis_matrix) {
  m_i <- sum(x_data_i)
  ilr_comp <- as.vector(pca$center + pca$rotation %*% scores)
  clr_comp <- ilr2clr(ilr_comp)
  composition <- clrInv_long(clr_comp)


  grad_vecs <- sapply(seq_along(scores), function(k) {
    e_k <- basis_matrix[k, ] 
    v_k <- pca$rotation[, k]
    term1 <- sum(x_data_i * e_k)
    term2 <- m_i * sum(composition * e_k)

    grad_k <- v_k * (term1 - term2) 

    return(grad_k)
  })
  grad <- rowSums(grad_vecs)

  grad <- grad - scores / (pca$sdev^2)

  return(grad)
}

gradient_cslc_vs2 <- function(scores,
                                 x_data_i,
                                 pca,
                                 basis_matrix,
                                 sc_factor) {
  scaling_factor <- sc_factor
  m_i <- sum(x_data_i)
  ilr_comp <- as.vector(pca$center + pca$rotation %*% scores)
  clr_comp <- ilr2clr(ilr_comp)
  composition <- clrInv_long(clr_comp)


  grad_vecs <- sapply(seq_along(scores), function(k) {
    e_k <- basis_matrix[k, ] 
    v_k <- pca$rotation[, k]
    term1 <- sum(x_data_i * e_k)
    term2 <- m_i * sum(composition * e_k)

    grad_k <- scaling_factor * v_k * (term1 - term2) - scores[k] / (pca$sdev[k]^2)

    return(grad_k)
  })
  grad <- rowSums(grad_vecs)

  return(grad)
}