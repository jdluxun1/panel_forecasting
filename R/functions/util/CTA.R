CTA <- function(Y, X, N, K, A_, sqrt_ht, iV, iVb_prior, PAI) {
  y_til <- Y %*% t(A_)
  for (j in seq_len(N)) {
    PAI[, j] <- 0
    j_to_N <- j:N
    lambda <- as.vector(sqrt_ht[, j_to_N, drop = FALSE])
    Y_j <- as.vector(y_til[, j_to_N, drop = FALSE] - X %*% PAI %*% t(A_[j_to_N, , drop = FALSE])) / lambda
    X_j <- kronecker(A_[j_to_N, j, drop = FALSE], X) / lambda
    index <- (K * (j - 1) + 1):(K * j)
    V_post <- solve(iV[index, index, drop = FALSE] + crossprod(X_j))
    b_post <- V_post %*% (iVb_prior[index] + crossprod(X_j, Y_j))
    PAI[, j] <- as.vector(b_post + t(chol(V_post)) %*% rnorm(K))
  }
  PAI
}
