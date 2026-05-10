mlag2 <- function(X, p) {
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0, Traw, N * p)
  for (ii in seq_len(p)) {
    Xlag[(p + 1):Traw, (N * (ii - 1) + 1):(N * ii)] <- X[(p + 1 - ii):(Traw - ii), , drop = FALSE]
  }
  Xlag
}
