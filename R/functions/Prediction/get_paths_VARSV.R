get_paths_VARSV <- function(PAI, h, L, N, PHI, inv_A, hlast, y) {
  hhat <- matrix(0, h + L, N)
  hhat[L, ] <- hlast[nrow(hlast), ]

  yhat <- matrix(0, h + L, N)
  yhat[seq_len(L), ] <- y[(nrow(y) - L + 1):nrow(y), , drop = FALSE]
  Cchol <- chol(PHI)

  for (m in (L + 1):(h + L)) {
    hhat[m, ] <- hhat[m - 1, ] + rnorm(N) %*% Cchol
    xhat <- numeric(0)
    for (i in seq_len(L)) xhat <- c(xhat, yhat[m - i, ])
    xhat <- c(1, xhat)
    smat <- inv_A %*% diag(exp(hhat[m, ]), N) %*% t(inv_A)
    yhat[m, ] <- xhat %*% matrix(PAI, nrow = N * L + 1, ncol = N) + rnorm(N) %*% chol_new(smat)
  }

  t(yhat[(L + 1):(h + L), , drop = FALSE])
}
