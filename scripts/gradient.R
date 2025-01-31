gradient_cslc_vs1 <- function(scores,
                              x_data_i,
                              pca,
                              basis_matrix,
                              sc_factor) {
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
  grad <- rowSums(grad_vecs) # TODO: shouldnt be colSums?

  grad <- grad - scores / (pca$sdev^2)

  return(grad)
}

gradient_cslc_vs1b <- function(scores,
                              x_data_i,
                              pca,
                              basis_matrix) {
  m_i <- sum(x_data_i)
  clr_comp <- as.vector(t(basis_matrix) %*% pca$center + t(basis_matrix) %*%  pca$rotation %*% scores)
  composition <- clrInverse(clr_comp)


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

gradient_cslc_vs1c <- function(scores,
                              x_data_i,
                              pca,
                              basis_matrix) {
  m_i <- sum(x_data_i)
  clr_comp <- as.vector(t(basis_matrix) %*% pca$center + t(basis_matrix) %*%  pca$rotation %*% scores)
  composition <- clrInverse(clr_comp)


  grad_vecs <- sapply(seq_along(scores), function(k) {
    e_k <- basis_matrix[k, ]
    v_k <- pca$rotation[, k]
    term1 <- sum(x_data_i * e_k)
    term2 <- m_i * sum(composition * e_k)

    grad_k <- v_k * (term1 - term2) 

    return(grad_k)
  })
  grad <- colSums(grad_vecs) 

  grad <- grad - scores / (pca$sdev^2)

  return(grad)
}

gradient_cslc_vs1d <- function(scores,
                              x_data_i,
                              pca,
                              basis_matrix) {
  m_i <- sum(x_data_i)
  clr_comp <- as.vector(t(basis_matrix) %*% pca$center + t(basis_matrix) %*%  pca$rotation %*% scores)
  composition <- clrInverse(clr_comp)


  grad_vecs <- sapply(seq_along(scores), function(k) {
    e_k <- basis_matrix[k, ]
    v_k <- pca$rotation[, k]
    term1 <- sum(x_data_i * e_k)
    term2 <- m_i * sum(composition * e_k)

    grad_k <- v_k * (term1 - term2) 

    return(grad_k)
  })
  grad <- colSums(grad_vecs) 

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

gradient_cslc_vs3 <- function(scores,
                              x_data_i,
                              pca,
                              basis_matrix,
                              sc_factor) {
  m_i <- sum(x_data_i)
  ilr_comp <- as.vector(pca$center + pca$rotation %*% scores)
  composition <- inv_ilr(ilr_comp, use_transpose = TRUE)

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

gradient_cslc_vs4 <- function(scores,
                              x_data_i,
                              pca,
                              basis_matrix) {
  m_i <- sum(x_data_i)
  ilr_comp <- as.vector(pca$center + pca$loading %*% scores)
  clr_comp <- ilr2clr(ilr_comp)
  composition <- clrInv_long(clr_comp)


  grad_vecs <- sapply(seq_along(scores), function(k) {
    e_k <- basis_matrix[k, ]
    v_k <- pca$loading[, k]
    term1 <- sum(x_data_i * e_k)
    term2 <- m_i * sum(composition * e_k)

    grad_k <- v_k * (term1 - term2) 

    return(grad_k)
  })
  grad <- rowSums(grad_vecs)

  grad <- grad - scores / (pca$sdev^2)

  return(grad)
}

gradient_lcs <- function(scores,
                              x_data_i,
                              pca,
                              basis_matrix) {
  m_i <- sum(x_data_i)
  clr_comp <- as.vector(
                        basis_matrix %*% pca$center +
                          basis_matrix %*%  pca$rotation %*% scores)
  composition <- clrInverse(clr_comp)


  grad_vecs <- sapply(seq_along(scores), function(k) {
    e_k <- basis_matrix[, k]
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

gradient_lcs_vs2 <- function(scores,
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
  grad <- rowSums(grad_vecs)

  grad <- grad - scores / (pca$sdev^2)

  return(grad)
}


