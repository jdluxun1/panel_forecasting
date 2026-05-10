chol_new <- function(A) {
  R <- try(chol(A), silent = TRUE)
  if (!inherits(R, "try-error")) return(R)
  chol(nearestSPD(A))
}
