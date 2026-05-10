nearestSPD <- function(A) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("A must be a square matrix.")
  if (length(A) == 1 && A <= 0) return(matrix(.Machine$double.eps, 1, 1))

  B <- (A + t(A)) / 2
  sv <- svd(B)
  H <- sv$v %*% diag(sv$d, nrow = length(sv$d)) %*% t(sv$v)
  Ahat <- (B + H) / 2
  Ahat <- (Ahat + t(Ahat)) / 2

  k <- 0
  while (inherits(try(chol(Ahat), silent = TRUE), "try-error")) {
    k <- k + 1
    mineig <- min(eigen(Ahat, symmetric = TRUE, only.values = TRUE)$values)
    Ahat <- Ahat + (-mineig * k^2 + .Machine$double.eps) * diag(nrow(A))
  }
  Ahat
}
