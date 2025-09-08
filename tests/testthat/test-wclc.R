test_that("wclc returns correct shape and metadata", {
  N <- 1000
  xy <- gen_lagged_pair(N, lag = 5, noise_sd = 0.02, seed = 42)
  res <- wclc(xy$x, xy$y, win_size = 200, win_inc = 50, max_lag = 20)

  expect_s3_class(res, "wclc")
  expect_type(res$R, "double")
  expect_equal(nrow(res$R), length(res$starts))
  expect_equal(ncol(res$R), length(res$lags))
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
  res <- wclc(xy$x, xy$y, win_size = 250, win_inc = 50, max_lag = 20)

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
  res <- wclc(xy$x, xy$y, win_size = 250, win_inc = 50, max_lag = 20)

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

  res <- wclc(x, y, win_size = 80, win_inc = 20, max_lag = 15)

  # Some entries should be NA (due to NA or zero variance)
  expect_true(any(is.na(res$R)))

  # If x is entirely constant, everything must be NA
  res_all_const <- wclc(rep(0, N), y, win_size = 80, win_inc = 20, max_lag = 15)
  expect_true(all(is.na(res_all_const$R)))
})

test_that("input validation errors are informative", {
  x <- 1:100; y <- 1:100
  expect_error(wclc(x, y[-1], 20, 10, 5), "same length", fixed = FALSE)
  expect_error(wclc(x, y, 0, 10, 5), "must be > 0", fixed = TRUE)
  expect_error(wclc(x, y, 20, -1, 5), "must be > 0", fixed = TRUE)
  expect_error(wclc(x, y, 200, 10, 5), "must be <= length", fixed = FALSE)
  expect_error(wclc(x, y, 60, 10, 30), "too short", fixed = FALSE)
})
