get_priorV2 <- function(omega, index_dom, index_di, N, G, L) {
  local_V <- as.matrix(omega)
  K <- nrow(local_V)
  NG <- ncol(local_V)
  sigmab <- matrix(0, K, NG)
  sigmab[1, ] <- local_V[1, ]

  for (ii in seq_len(N)) {
    indx_v <- ((ii - 1) * G + 1):(ii * G)
    indexc <- index_dom[ii, ]
    for (jj in seq_len(G)) {
      inda <- seq_len(G * L)
      ind <- seq(jj, G * L, by = G)
      inda <- inda[-ind]
      sigmab[indexc[ind], indx_v[jj]] <- local_V[indexc[ind], indx_v[jj]]
      sigmab[indexc[inda], indx_v[jj]] <- local_V[indexc[inda], indx_v[jj]]
    }
  }

  for (ii in seq_len(NG)) {
    inc <- ceiling(ii / G)
    ind_di <- index_di[inc, ]
    sigmab[ind_di, ii] <- local_V[ind_di, ii]
  }
  sigmab
}
