get_hypall <- function(int, olag, cvlag, di, index_dom, index_for, N, G, L) {
  NG <- N * G
  K <- 1 + NG * L
  olag_mat <- matrix(olag, nrow = L, ncol = NG)
  cvlag_mat <- matrix(cvlag, nrow = (G - 1) * L, ncol = NG)
  k_d <- G * L * (N - 1)
  parameters <- matrix(0, K, NG)
  parameters[1, ] <- int

  for (ii in seq_len(N)) {
    indx_v <- ((ii - 1) * G + 1):(ii * G)
    olag_temp <- olag_mat[, indx_v, drop = FALSE]
    cvlag_temp <- cvlag_mat[, indx_v, drop = FALSE]
    indexc <- index_dom[ii, ]
    for (jj in seq_len(G)) {
      inda <- seq_len(G * L)
      ind <- seq(jj, G * L, by = G)
      inda <- inda[-ind]
      parameters[indexc[ind], indx_v[jj]] <- olag_temp[, jj]
      parameters[indexc[inda], indx_v[jj]] <- cvlag_temp[, jj]
    }
  }

  for (ii in seq_len(N)) {
    indx_v <- ((ii - 1) * G + 1):(ii * G)
    di_temp <- matrix(di[, ii], nrow = k_d, ncol = G)
    parameters[index_for[ii, ], indx_v] <- di_temp
  }
  parameters
}
