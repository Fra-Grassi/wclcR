test_that("wclc returns correct shape and metadata", {
  N <- 1000
  xy <- gen_lagged_pair(N, lag = 5, noise_sd = 0.02, seed = 42)
  res <- wclc(xy$x, xy$y, win_size = 200, win_inc = 50, max_lag = 20, output = "both")

  expect_s3_class(res, "wclc")
  expect_type(res$R, "double")
  expect_type(res$Z, "double")
  expect_equal(nrow(res$R), length(res$starts))
  expect_equal(ncol(res$R), length(res$lags))
  expect_equal(nrow(res$Z), length(res$starts))
  expect_equal(ncol(res$Z), length(res$lags))
  expect_equal(res$lags, -20:20)
  expect_true(all(diff(res$starts) %in% 50))
  # Windows should be in-bounds given the constraints
  expect_true(min(res$starts) >= 20 + 1)
  expect_true(max(res$starts) + 200 - 1 <= N - 20)
})

test_that("wclc recovers a known positive lag (within tolerance)", {
  skip_on_cran()
  N <- 1200
  true_lag <- 7L
  xy <- gen_lagged_pair(N, lag = true_lag, noise_sd = 0.03, seed = 123)
  res <- wclc(xy$x, xy$y, win_size = 250, win_inc = 50, max_lag = 20, output = "r")

  # Per-window peak lag (ties -> first)
  peak_lags <- res$lags[max.col(res$R, ties.method = "first")]
  # Be robust to local noise: median error small
  med_err <- median(peak_lags - true_lag, na.rm = TRUE)
  expect_lte(abs(med_err), 2)
  # And most windows should be near the true lag
  frac_near <- mean(abs(peak_lags - true_lag) <= 2, na.rm = TRUE)
  expect_gt(frac_near, 0.6)
})

test_that("wclc recovers a known negative lag (y leads x)", {
  skip_on_cran()
  N <- 1200
  true_lag <- -6L
  xy <- gen_lagged_pair(N, lag = true_lag, noise_sd = 0.03, seed = 321)
  res <- wclc(xy$x, xy$y, win_size = 250, win_inc = 50, max_lag = 20, output = "r")

  peak_lags <- res$lags[max.col(res$R, ties.method = "first")]
  med_err <- median(peak_lags - true_lag, na.rm = TRUE)
  expect_lte(abs(med_err), 2)
})

test_that("NA handling and zero-variance windows yield NA correlations", {
  N <- 600
  xy <- gen_lagged_pair(N, lag = 4, noise_sd = 0.02, seed = 99)
  x <- xy$x; y <- xy$y

  # Inject NA segment that will contaminate some windows
  x[101:180] <- NA_real_  # 80 points == win_size
  # Make a long constant segment to force zero variance in several windows
  y[300:380] <- 1  # 81 points > win_size

  res <- wclc(x, y, win_size = 80, win_inc = 20, max_lag = 15, output = "r")

  # Some entries should be NA (due to NA or zero variance)
  expect_true(any(is.na(res$R)))

  # If x is entirely constant, everything must be NA
  res_all_const <- wclc(rep(0, N), y, win_size = 80, win_inc = 20, max_lag = 15)
  expect_true(all(is.na(res_all_const$R)))
})

test_that("Fisher z output matches atanh(r) with clampiing", {
  set.seed(1)
  x <- rnorm(1000); y <- x + rnorm(1000, sd = 0.2)
  a <- wclc(x, y, win_size = 100, win_inc = 50, max_lag = 10, output = "both", z_eps = 1e-9)
  # Construct expected Z by clamping then atanh
  R <- a$R
  Rcl <- pmax(pmin(R, 1 - 1e-9), -1 + 1e-9)
  expect_equal(a$Z, atanh(Rcl), tolerance = 1e-12)
})

test_that("when z_eps = NA, Z equals atanh(R) exactly (within FP tolerance)", {
  set.seed(1)
  x <- rnorm(600); y <- x + rnorm(600, sd = 0.2)
  res <- wclc(x, y, win_size = 150, win_inc = 50, max_lag = 6, output = "both", z_eps = NA)
  expect_equal(res$Z, atanh(res$R), tolerance = 1e-12)
})

test_that("z_eps = NA disables clamping (z grows or stays same, equals atanh(R))", {
  # Construct strongly correlated identical signals
  x <- rep(1:50, 10)
  y <- x

  # No clamping
  b_nc <- wclc(x, y, win_size = 100, win_inc = 50, max_lag = 0,
               output = "both", z_eps = NA)
  # With clamping
  b_cl <- wclc(x, y, win_size = 100, win_inc = 50, max_lag = 0,
               output = "z", z_eps = 1e-9)

  # 1) No-clamp Z must equal atanh(R) exactly (modulo floating error)
  expect_equal(b_nc$Z, atanh(b_nc$R), tolerance = 1e-12)

  # 2) Removing clamping cannot reduce |Z|
  expect_gte(max(abs(b_nc$Z), na.rm = TRUE), max(abs(b_cl$Z), na.rm = TRUE))
})

test_that("output='r' returns only R, not Z", {
  set.seed(1)
  x <- rnorm(500); y <- x + rnorm(500, sd = 0.3)
  res <- wclc(x, y, win_size = 100, win_inc = 50, max_lag = 5, output = "r")
  expect_s3_class(res, "wclc")
  expect_true(!is.null(res$R))
  expect_null(res$Z)
  expect_type(res$R, "double")
})

test_that("output='z' returns only Z, not R", {
  set.seed(1)
  x <- rnorm(500); y <- x + rnorm(500, sd = 0.3)
  res <- wclc(x, y, win_size = 100, win_inc = 50, max_lag = 5, output = "z")
  expect_s3_class(res, "wclc")
  expect_true(!is.null(res$Z))
  expect_null(res$R)
  expect_type(res$Z, "double")
})

test_that("output='both' returns both R and Z and Z == atanh(clamped R)", {
  set.seed(1)
  x <- rnorm(800); y <- x + rnorm(800, sd = 0.25)
  z_eps <- 1e-9
  res <- wclc(x, y, win_size = 120, win_inc = 40, max_lag = 8, output = "both", z_eps = z_eps)
  expect_s3_class(res, "wclc")
  expect_true(!is.null(res$R))
  expect_true(!is.null(res$Z))
  Rcl <- pmax(pmin(res$R, 1 - z_eps), -1 + z_eps)
  expect_equal(res$Z, atanh(Rcl), tolerance = 1e-12)
})

test_that("input validation errors are informative", {
  x <- 1:100; y <- 1:100
  expect_error(wclc(x, y[-1], 20, 10, 5, output = "r"), "same length", fixed = FALSE)
  expect_error(wclc(x, y, 0, 10, 5, output = "r"), "must be > 0", fixed = TRUE)
  expect_error(wclc(x, y, 20, -1, 5, output = "r"), "must be > 0", fixed = TRUE)
  expect_error(wclc(x, y, 200, 10, 5, output = "r"), "must be <= length", fixed = FALSE)
  expect_error(wclc(x, y, 60, 10, 30, output = "r"), "too short", fixed = FALSE)
  expect_error(wclc(x, y, 20, 10, 5, output = "nonsense"), "should be one of", fixed = FALSE)
})

test_that("edges='ragged' yields NA tails per-lag (no padding, full starts)", {
  set.seed(1)
  N <- 100; win_size <- 20; win_inc <- 5; max_lag <- 10
  x <- rnorm(N); y <- x + rnorm(N, 0.2)

  b <- wclc(x, y, win_size, win_inc, max_lag, edges = "ragged", output = "r")

  # 1) Full start coverage: from 1 to N - win_size + 1
  expect_equal(min(b$starts), 1L)
  expect_equal(max(b$starts), N - win_size + 1L)

  # Helper to compute expected number of invalid (NA) starts for a lag
  expected_na_count <- function(L) {
    bound <- if (L >= 0L) (N - win_size + 1L - L) else (N - win_size + 1L + L)
    sum(b$starts > bound)
  }

  # Column indices
  col0  <- which(b$lags == 0)
  colP  <- which(b$lags ==  max_lag)
  colN  <- which(b$lags == -max_lag)

  # 2) Lag 0: no NAs anywhere
  expect_false(any(is.na(b$R[, col0])))

  # 3) Positive max lag: only late starts should be NA; count must match formula
  na_pos <- is.na(b$R[, colP])
  expect_false(any(na_pos[1:3]))  # no early NAs
  expect_true(any(na_pos))        # there are some at the tail
  expect_equal(sum(na_pos), expected_na_count(+max_lag))

  # 4) Negative max lag: only late starts should be NA; count matches
  na_neg <- is.na(b$R[, colN])
  expect_false(any(na_neg[1:3]))  # no early NAs
  expect_true(any(na_neg))        # tail NAs exist
  expect_equal(sum(na_neg), expected_na_count(-max_lag))
})
