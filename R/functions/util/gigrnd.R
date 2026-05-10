gigrnd <- function(P, a, b, sampleSize = 1L) {
  if (a <= 0) stop("a must be positive")
  if (b <= 0) b <- .Machine$double.eps

  lambda <- P
  omega <- sqrt(a * b)
  swap <- FALSE
  if (lambda < 0) {
    lambda <- -lambda
    swap <- TRUE
  }
  alpha <- sqrt(omega^2 + lambda^2) - lambda

  psi <- function(x, alpha, lambda) -alpha * (cosh(x) - 1) - lambda * (exp(x) - x - 1)
  dpsi <- function(x, alpha, lambda) -alpha * sinh(x) - lambda * (exp(x) - 1)
  gfun <- function(x, sd, td, f1, f2) {
    if (x >= -sd && x <= td) return(1)
    if (x > td) return(f1)
    f2
  }

  x <- -psi(1, alpha, lambda)
  if (x >= 0.5 && x <= 2) {
    tt <- 1
  } else if (x > 2) {
    tt <- sqrt(2 / (alpha + lambda))
  } else {
    tt <- log(4 / (alpha + 2 * lambda))
  }

  x <- -psi(-1, alpha, lambda)
  if (x >= 0.5 && x <= 2) {
    s <- 1
  } else if (x > 2) {
    s <- sqrt(4 / (alpha * cosh(1) + lambda))
  } else {
    s <- min(if (lambda == 0) Inf else 1 / lambda,
             log(1 + 1 / alpha + sqrt(1 / alpha^2 + 2 / alpha)))
  }

  eta <- -psi(tt, alpha, lambda)
  zeta <- -dpsi(tt, alpha, lambda)
  theta <- -psi(-s, alpha, lambda)
  xi <- dpsi(-s, alpha, lambda)
  p <- 1 / xi
  r <- 1 / zeta
  td <- tt - r * eta
  sd <- s - p * theta
  q <- td + sd
  X <- numeric(sampleSize)

  for (sample in seq_len(sampleSize)) {
    done <- FALSE
    while (!done) {
      U <- runif(1)
      V <- runif(1)
      W <- runif(1)
      if (U < q / (p + q + r)) {
        X[sample] <- -sd + q * V
      } else if (U < (q + r) / (p + q + r)) {
        X[sample] <- td - r * log(V)
      } else {
        X[sample] <- -sd + p * log(V)
      }

      f1 <- exp(-eta - zeta * (X[sample] - tt))
      f2 <- exp(-theta + xi * (X[sample] + s))
      if (W * gfun(X[sample], sd, td, f1, f2) <= exp(psi(X[sample], alpha, lambda))) {
        done <- TRUE
      }
    }
  }

  X <- exp(X) * (lambda / omega + sqrt(1 + (lambda / omega)^2))
  if (swap) X <- 1 / X
  X / sqrt(a / b)
}
