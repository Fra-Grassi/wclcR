#' Windowed Cross-Lagged Correlation (Boker-style, split-lag, NA-safe)
#'
#' Computes windowed cross-lagged Pearson correlations between two numeric
#' vectors using “split-lag” windows (positive lags shift \code{y}, negative lags
#' shift \code{x}). Within each window, lags are scanned in 1-sample steps.
#'
#' @param x,y Numeric vectors of equal length.
#' @param win_size Integer > 0. Window length in samples.
#' @param win_inc  Integer > 0. Window increment (hop) in samples.
#' @param max_lag  Integer >= 0. Maximum absolute lag in samples.
#' @param output One of c("r", "z", "both"). If "z" or "both", Fisher z-transform is applied via atanh() to (optionally clamped) r.
#' @param z_eps Small epsilon to use to clamp e to (-1+z_eps, 1-z_eps) before atanh to avoid infinite values. Default to 1e-12; set to NA to disable clamping.
#'
#' @return An object of class \code{"wclc"}: a list with
#' \itemize{
#'   \item \code{R}     — matrix with dimensions (num_windows x (2*max_lag+1)) of Pearson \eqn{r} (present if output %in% c("r", "both"))
#'   \item \code{Z}     — matrix with dimensions (num_windows x (2*max_lag+1)) of Fisher z = atanh(\eqn{r}) (present if output %in% c("z", "both"))
#'   \item \code{lags}  — integer vector of lags (-max_lag … +max_lag)
#'   \item \code{starts} — 1-based window start indices in the original series
#' }
#' @references Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological Methods, 7(3), 338–355. https://doi.org/10.1037/1082-989X.7.3.338

#' @examples
#' #' set.seed(1)
#' N <- 1000
#' x <- sin(2*pi*2*(0:(N-1))/100) + rnorm(N, 0, 0.05)
#' y <- c(rep(NA_real_, 5), x[1:(N-5)]) + rnorm(N, 0, 0.05)
#' res <- wclc(x, y, win_size = 200, win_inc = 50, max_lag = 20, output = "both", z_eps = 1e-12)
#' @export
wclc <- function(x, y, win_size, win_inc, max_lag, output = c("r", "z", "both"), z_eps = 1e-12)  {
  output <- match.arg(output)
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors.", call. = FALSE)
  }
  if (length(x) != length(y)) {
    stop("x and y must have the same length.", call. = FALSE)
  }
  N <- length(x)
  if (win_size <= 0L || win_inc <= 0L || max_lag < 0L) {
    stop("win_size, win_inc must be > 0, and max_lag >= 0.", call. = FALSE)
  }
  if (win_size > N) {
    stop("win_size must be <= length(x).", call. = FALSE)
  }
  if (2L * max_lag + win_size > N) {
    stop("Series too short for given win_size and max_lag.", call. = FALSE)
  }

  lags <- -max_lag:max_lag
  nL   <- length(lags)

  # Valid starts so all lagged windows stay in-bounds
  first_start <- max_lag + 1L
  last_start  <- N - max_lag - win_size + 1L
  if (last_start < first_start) stop("No valid window fits with these parameters.", call. = FALSE)
  starts <- seq.int(first_start, last_start, by = win_inc)
  nW <- length(starts)

  # Precompute base window indices: matrix [win_size x nW]
  # idx_mat[r, c] = original index for row r of window starting at starts[c]
  idx_mat <- outer(0:(win_size - 1L), starts, `+`)          # rows: 1..win_size, cols: windows

  # Column-wise, NA-safe Pearson correlation
  col_cor_na <- function(A, B) {
    # A, B: matrices with same dims [win_size x nW]
    vapply(seq_len(ncol(A)), function(j) {
      aj <- A[, j]; bj <- B[, j]
      ok <- is.finite(aj) & is.finite(bj)
      nv <- sum(ok)
      if (nv < 2L) return(NA_real_)
      aj <- aj[ok]; bj <- bj[ok]
      sx <- stats::sd(aj); sy <- stats::sd(bj)
      if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(NA_real_)
      ajc <- aj - mean(aj); bjc <- bj - mean(bj)
      sum(ajc * bjc) / ((nv - 1) * sx * sy)
    }, numeric(1))
  }

  # Allocate result: rows = windows (time), cols = lags
  R <- matrix(NA_real_, nrow = nW, ncol = nL,
              dimnames = list(start = starts, lag = as.character(lags)))

  # Iterate over lags (cheap), build lagged index mats, pull slices, correlate by columns
  for (i in seq_along(lags)) {
    L <- lags[i]
    if (L >= 0L) {
      Wx <- matrix(x[idx_mat], nrow = win_size)
      Wy <- matrix(y[idx_mat + L], nrow = win_size)
    } else {
      Wx <- matrix(x[idx_mat - L], nrow = win_size)  # note: -L is positive shift
      Wy <- matrix(y[idx_mat], nrow = win_size)
    }
    # Column-wise correlation -> a length-nW vector; store as column i
    R[, i] <- col_cor_na(Wx, Wy)
  }

  # Define return list based on output selection
  out <- list(lags = lags, starts = starts)
  if (output %in% c("r", "both")) {  # add computed correlations
    out$R <- R
  }
  if (output %in% c("z", "both")) {  # add z-transformed correlations
    if (is.na(z_eps)) {
      Z <- atanh(R)  # no clamping, may yield ±Inf if |r|==1
    } else {
      Z <- atanh(pmax(pmin(R, 1 - z_eps), - 1 + z_eps))  # clamping
    }
    out$Z <- Z
  }

  class(out) <- "wclc"

  return(out)
}
