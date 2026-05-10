get_grhyp <- function(parameters, index_dom, index_di, N, G, L) {
  parameters <- as.matrix(parameters)
  NG <- N * G
  k2 <- G * L * (N - 1)
  int <- parameters[1, ]
  olag <- numeric(0)
  cvlag <- numeric(0)
  ditemp <- numeric(0)

  for (ii in seq_len(N)) {
    indx_v <- ((ii - 1) * G + 1):(ii * G)
    indexc <- index_dom[ii, ]
    for (jj in seq_len(G)) {
      inda <- seq_len(G * L)
      ind <- seq(jj, G * L, by = G)
      inda <- inda[-ind]
      olag <- c(olag, parameters[indexc[ind], indx_v[jj]])
      cvlag <- c(cvlag, parameters[indexc[inda], indx_v[jj]])
    }
  }

  for (ii in seq_len(NG)) {
    inc <- ceiling(ii / G)
    ind_di <- index_di[inc, ]
    ditemp <- c(ditemp, parameters[ind_di, ii])
  }

  di <- matrix(ditemp, nrow = k2 * G, ncol = N)
  list(int = int, olag = olag, cvlag = cvlag, di = di)
}
