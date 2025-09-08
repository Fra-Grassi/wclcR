
# wclcR

<!-- badges: start -->

<!-- badges: end -->

`wclcR` provides tools for computing windowed cross-lagged correlations
(Boker et al., 2002) in R.  
It implements a fast, NA-safe split-lag algorithm and will grow to
include peak extraction, plotting, and EEG/EMG use-cases.

## Installation

You can install the development version of wclcR like so:

``` r
# install.packages("remotes")
remotes::install_github("Fra-Grassi/wclcR")
```

## Example

``` r
set.seed(1)
N <- 1000
x <- sin(2*pi*2*(0:(N-1))/100) + rnorm(N, 0, 0.05)
y <- c(rep(NA_real_, 5), x[1:(N-5)]) + rnorm(N, 0, 0.05)

library(wclcR)
res <- wclc(x, y, win_size = 200, win_inc = 50, max_lag = 20, output = "r")
```
