KSC <- function(v_tilda, h, Qt, states_pmean, states_pvar) {
  v_tilda <- as.matrix(v_tilda)
  h <- as.matrix(h)
  Tobs <- nrow(v_tilda)
  N <- ncol(v_tilda)

  pi_mix <- c(0.00609, .04775, .13057, .20674, .22715, .18842, .12047, .05591, .01575, .00115)
  mi <- c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
  si <- c(0.11265, 0.177788, 0.26768, 0.40611, 0.62699, 0.98583, 1.57469, 2.54498, 4.16591, 7.33342)

  S <- matrix(0L, Tobs, N)
  for (i in seq_len(N)) {
    q <- matrix(pi_mix, Tobs, 10, byrow = TRUE) *
      dnorm(matrix(v_tilda[, i], Tobs, 10),
            mean = matrix(h[, i], Tobs, 10) + matrix(mi, Tobs, 10, byrow = TRUE),
            sd = matrix(sqrt(si), Tobs, 10, byrow = TRUE))
    q <- q / rowSums(q)
    S[, i] <- 10 - rowSums(matrix(runif(Tobs), Tobs, 10) < t(apply(q, 1, cumsum))) + 1
  }

  CK(v_tilda - matrix(mi[S], Tobs, N), matrix(si[S], Tobs, N), Qt, N, Tobs, states_pmean, states_pvar)
}

CK <- function(y, Ht, Qt, N, Tobs, S0, P0) {
  y <- t(y)
  St_collect <- matrix(0, N, Tobs)
  Pt_collect <- array(0, dim = c(N, N, Tobs))
  St_draw <- matrix(0, N, Tobs)

  St <- as.matrix(S0)
  Pt <- as.matrix(P0)
  for (tidx in seq_len(Tobs)) {
    St_1 <- St
    Pt_1 <- Pt + Qt
    vt <- y[, tidx, drop = FALSE] - St_1
    Varvt <- Pt_1 + diag(Ht[tidx, ], N)
    Kt <- Pt_1 %*% solve(Varvt)
    St <- St_1 + Kt %*% vt
    Pt <- Pt_1 - Kt %*% Pt_1
    St_collect[, tidx] <- St
    Pt_collect[, , tidx] <- Pt
  }

  St_draw[, Tobs] <- St_collect[, Tobs] + t(chol(Pt_collect[, , Tobs])) %*% rnorm(N)
  if (Tobs > 1) {
    for (tidx in seq(Tobs - 1, 1)) {
      Kt <- Pt_collect[, , tidx] %*% solve(Pt_collect[, , tidx] + Qt)
      Smean <- St_collect[, tidx, drop = FALSE] + Kt %*% (St_draw[, tidx + 1, drop = FALSE] - St_collect[, tidx, drop = FALSE])
      Svar <- Pt_collect[, , tidx] - Kt %*% Pt_collect[, , tidx]
      St_draw[, tidx] <- Smean + t(chol(Svar)) %*% rnorm(N)
    }
  }
  t(St_draw)
}
