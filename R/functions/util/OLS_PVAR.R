OLS_PVAR <- function(Yraw, N, G, p) {
  Yraw <- as.matrix(Yraw)
  Traw <- nrow(Yraw)
  NG <- ncol(Yraw)
  if (NG != N * G) stop("wrong specification of N and G")

  Ylag <- mlag2(Yraw, p)
  k <- p * NG
  Tobs <- Traw - p
  X <- Ylag[(p + 1):Traw, , drop = FALSE]
  x <- kronecker(diag(NG), X)
  Y <- Yraw[(p + 1):Traw, , drop = FALSE]
  y <- as.vector(Y)

  alpha_OLS_vec <- solve(crossprod(x), crossprod(x, y))
  alpha_OLS_mat <- solve(crossprod(X), crossprod(X, Y))
  err <- Y - X %*% alpha_OLS_mat
  SSE <- crossprod(err)
  sigma_OLS <- SSE / (Tobs - (k - 1))

  list(alpha_OLS_vec = as.vector(alpha_OLS_vec), sigma_OLS = sigma_OLS, err = err)
}
