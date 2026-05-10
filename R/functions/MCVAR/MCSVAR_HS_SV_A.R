MCSVAR_HS_SV_A <- function(data, reps, burnin, N, G, L, h) {
  data <- as.matrix(data)
  Traw <- nrow(data)
  NG <- ncol(data)
  Tobs <- Traw - L
  K_eq <- NG * L + 1

  nm1 <- NG * L
  nm2 <- (G - 1) * L * NG
  nd <- G * L * (N - 1)
  m <- NG * (NG - 1) / 2
  alp <- numeric(m)
  fsize <- reps - burnin
  Y_pred <- array(0, dim = c(fsize, NG, h))

  Y <- data[(L + 1):Traw, , drop = FALSE]
  Y_lag <- mlag2(data, L)
  X <- cbind(1, Y_lag[(L + 1):Traw, , drop = FALSE])

  index_dom <- matrix(0, N, G * L + 1)
  index_for <- matrix(0, N, (NG - G) * L)
  for (i in seq_len(N)) {
    for (j in seq_len(L)) {
      index <- ((j - 1) * NG + 2):(j * NG + 1)
      index_dom[i, ((j - 1) * G + 2):(j * G + 1)] <- ((i - 1) * G + (j - 1) * NG + 2):(i * G + (j - 1) * NG + 1)
      index <- index[-(((i - 1) * G + 1):(i * G))]
      index_for[i, ((j - 1) * (NG - G) + 2):(j * (NG - G) + 1)] <- index
    }
  }
  index_dom <- index_dom[, -1, drop = FALSE]
  index_for <- index_for[, -1, drop = FALSE]

  Omega <- matrix(0.0001, K_eq, NG)
  Gamma <- matrix(1, K_eq, NG)
  gamma_blocks <- get_grhyp(Gamma, index_dom, index_for, N, G, L)
  Gamma_int <- as.vector(t(gamma_blocks$int))
  Gamma_olag <- gamma_blocks$olag
  Gamma_cvlag <- gamma_blocks$cvlag
  Gamma_di <- gamma_blocks$di
  Omega_a <- rep(1, m)
  Gamma_a <- rep(1, m)
  lambda <- rep(1, 2)
  c_hyp <- rep(1, 2)
  Sigma_b <- get_priorV2(Omega, index_dom, index_for, N, G, L)
  iV <- diag(1 / as.vector(Sigma_b), K_eq * NG)
  Sigma_a <- Omega_a

  Vol_0mean <- matrix(0, NG, 1)
  Vol_0var <- 100 * diag(NG)

  ARresid <- matrix(0, Traw - 1, NG)
  for (i in seq_len(NG)) {
    ARresid[, i] <- OLS_PVAR(data[, i, drop = FALSE], 1, 1, 1)$err
  }
  htemp <- rowMeans(t(ARresid)^2)
  sqrt_ht <- matrix(rep(sqrt(htemp), each = Tobs), Tobs, NG)
  Vol_states <- 2 * log(sqrt_ht)
  PHI_ <- 0.0001 * diag(NG)

  d_PHI <- NG + 2
  s_PHI <- 0.01 * diag(NG)
  A_ <- diag(NG)

  comp <- if (L > 1) cbind(diag(NG * (L - 1)), matrix(0, NG * (L - 1), NG)) else matrix(0, 0, NG)

  for (irep in seq_len(reps)) {
    if (irep == 1) PAI <- solve(crossprod(X), crossprod(X, Y))

    if (irep <= burnin) {
      PAI <- CTA(Y, X, NG, K_eq, A_, sqrt_ht, iV, rep(0, K_eq * NG), PAI)
    } else {
      stationary <- FALSE
      while (!stationary) {
        PAI <- CTA(Y, X, NG, K_eq, A_, sqrt_ht, iV, rep(0, K_eq * NG), PAI)
        companion <- rbind(t(PAI[2:K_eq, , drop = FALSE]), comp)
        if (max(abs(eigen(companion, only.values = TRUE)$values)) < 1) stationary <- TRUE
      }
    }
    RESID <- Y - X %*% PAI

    count <- 0
    ivaprior <- diag(1 / Sigma_a, m)
    for (ii in 2:NG) {
      y_spread_adj <- RESID[, ii] / sqrt_ht[, ii]
      X_spread_adj <- matrix(0, Tobs, ii - 1)
      for (vv in seq_len(ii - 1)) X_spread_adj[, vv] <- RESID[, vv] / sqrt_ht[, ii]
      ZZ <- crossprod(X_spread_adj)
      Zz <- crossprod(X_spread_adj, y_spread_adj)
      idx <- (count + 1):(count + ii - 1)
      Valpha_post <- solve(ZZ + ivaprior[idx, idx, drop = FALSE])
      alpha_post <- Valpha_post %*% Zz
      alphadraw <- alpha_post + t(chol(Valpha_post)) %*% rnorm(ii - 1)
      a1 <- -as.vector(alphadraw)
      A_[ii, 1:(ii - 1)] <- a1
      alp[idx] <- a1
      count <- count + ii - 1
    }
    invA_ <- solve(A_)

    omega_blocks <- get_grhyp(Omega, index_dom, index_for, N, G, L)
    Omega_int <- omega_blocks$int
    Omega_olag <- omega_blocks$olag
    Omega_cvlag <- omega_blocks$cvlag
    Omega_di <- omega_blocks$di
    beta_blocks <- get_grhyp(PAI, index_dom, index_for, N, G, L)
    Beta_int <- beta_blocks$int
    Beta_olag <- beta_blocks$olag
    Beta_cvlag <- beta_blocks$cvlag
    Beta_di <- beta_blocks$di

    for (ie in seq_len(NG)) Omega_int[ie] <- gigrnd(0, 2 * Gamma_int[ie], Beta_int[ie]^2, 1)
    for (ie in seq_len(nm1)) Omega_olag[ie] <- gigrnd(0, 2 * Gamma_olag[ie], Beta_olag[ie]^2, 1)
    for (ie in seq_len(nm2)) Omega_cvlag[ie] <- gigrnd(0, 2 * Gamma_cvlag[ie], Beta_cvlag[ie]^2, 1)
    for (jt in seq_len(N)) {
      for (ie in seq_len(nd * G)) Omega_di[ie, jt] <- gigrnd(0, 2 * Gamma_di[ie, jt], Beta_di[ie, jt]^2, 1)
    }
    for (ie in seq_len(m)) Omega_a[ie] <- gigrnd(0, 2 * Gamma_a[ie], alp[ie]^2, 1)

    for (ie in seq_len(NG)) Gamma_int[ie] <- rgamma(1, shape = 1, scale = 1 / (lambda[1] + Omega_int[ie]))
    for (ie in seq_len(nm1)) Gamma_olag[ie] <- rgamma(1, shape = 1, scale = 1 / (lambda[1] + Omega_olag[ie]))
    for (ie in seq_len(nm2)) Gamma_cvlag[ie] <- rgamma(1, shape = 1, scale = 1 / (lambda[1] + Omega_cvlag[ie]))
    for (jt in seq_len(N)) {
      for (ie in seq_len(nd * G)) Gamma_di[ie, jt] <- rgamma(1, shape = 1, scale = 1 / (lambda[1] + Omega_di[ie, jt]))
    }
    for (ie in seq_len(m)) Gamma_a[ie] <- rgamma(1, shape = 1, scale = 1 / (lambda[2] + Omega_a[ie]))

    lambda[1] <- rgamma(1, shape = 1 + (K_eq * NG - 1) / 2, scale = 1 / (c_hyp[1] + sum(c(Gamma_int, Gamma_olag, Gamma_cvlag, as.vector(Gamma_di)))))
    c_hyp[1] <- rgamma(1, shape = 1, scale = 1 / (lambda[1] + 1))
    lambda[2] <- rgamma(1, shape = 1 + (m - 1) / 2, scale = 1 / (c_hyp[2] + sum(Gamma_a)))
    c_hyp[2] <- rgamma(1, shape = 1, scale = 1 / (lambda[2] + 1))

    Omega <- get_hypall(Omega_int, Omega_olag, Omega_cvlag, Omega_di, index_dom, index_for, N, G, L)
    Sigma_b <- get_priorV2(Omega, index_dom, index_for, N, G, L)
    Sigma_a <- Omega_a
    Sigma_b[Sigma_b < 1e-10] <- 1e-10
    iV <- diag(1 / as.vector(Sigma_b), K_eq * NG)
    Sigma_a[Sigma_a < 1e-10] <- 1e-10

    Vol_states <- KSC(log((RESID %*% t(A_))^2 + 1e-6), Vol_states, PHI_, Vol_0mean, Vol_0var)
    sqrt_ht <- exp(Vol_states / 2)

    eta <- Vol_states[2:nrow(Vol_states), , drop = FALSE] - Vol_states[1:(nrow(Vol_states) - 1), , drop = FALSE]
    temp <- t(chol(solve(s_PHI + crossprod(eta)))) %*% matrix(rnorm(NG * (Tobs + d_PHI)), NG, Tobs + d_PHI)
    PHI_ <- solve(temp %*% t(temp))

    if (irep > burnin) {
      Y_pred[irep - burnin, , ] <- get_paths_VARSV(as.vector(PAI), h, L, NG, PHI_, invA_, Vol_states, Y)
    }
  }
  Y_pred
}

MCVAR_HS_SV_A <- MCSVAR_HS_SV_A
