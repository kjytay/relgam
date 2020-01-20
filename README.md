<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
[![CRAN
status](https://www.r-pkg.org/badges/version/relgam)](https://cran.r-project.org/package=relgam)
<!-- badges: end -->

relgam
======

`relgam` is a package implementing **Reluctant Generalized Additive
Modeling (RGAM)**, a new technique for fitting sparse generalized
additive models (GAMs). RGAMs are computationally scalable and work with
continuous, binary, count and survival data. For more details on RGAM,
please see [the preprint](https://arxiv.org/abs/1912.01808). At a high
level, RGAM is fit in the following way: letting `y` denote the response
variable and `X` denote the design matrix,

1.  Fit the lasso of `y` on `X` using the `glmnet` package. Compute the
    residuals `r` at the `lambda` hyperparameter selected by
    cross-validation.
2.  For each of the original features, fit a `d`-degree smoothing spline
    of `r` on the feature to get a new, non-linear feature. Let `F`
    denote the matrix of new features.
3.  Fit the lasso of `y` on `X` and `F` using the `glmnet` package.

An example
----------

Here is a simple example to illustrate how to use this package. First,
let’s generate some data:

``` r
set.seed(1)
n <- 100; p <- 12
x = matrix(rnorm((n) * p), ncol = p)
f4 = 2 * x[,4]^2 + 4 * x[,4] - 2
f5 = -2 * x[, 5]^2 + 2
f6 = 0.5 * x[, 6]^3
mu = rowSums(x[, 1:3]) + f4 + f5 + f6
y = mu + sqrt(var(mu) / 4) * rnorm(n)
```

We can fit an RGAM using the `rgam()` function:

``` r
library(relgam)
fit <- rgam(x, y, verbose = FALSE)
#  init_nz not specified: setting to default (all features)
#  using default value of gamma for RGAM: 0.6
```

(If `verbose = TRUE` (default), model-fitting is tracked with a progress
bar in the console.)
