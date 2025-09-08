# Helpers for tests -----------------------------------------------

gen_lagged_pair <- function(N, lag, noise_sd = 0.05, seed = 1) {
  set.seed(seed)
  t <- seq_len(N)
  x <- sin(2 * pi * 2 * t / 100) + rnorm(N, 0, noise_sd)
  if (lag >= 0) {
    y <- c(rep(NA_real_, lag), x[1:(N - lag)]) + rnorm(N, 0, noise_sd)
  } else {
    L <- -lag
    y <- c(x[(L + 1):N], rep(NA_real_, L)) + rnorm(N, 0, noise_sd)
  }
  list(x = x, y = y)
}
