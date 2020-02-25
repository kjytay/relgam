context("Match v1.0 rgam()")
library(relgam)

# set up data
set.seed(1)
n <- 100; p <- 20
x <- matrix(rnorm(n * p), n, p)
beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
y <- x %*% beta + rnorm(n)
bin_y <- ifelse(y > 0, 1, 0)
poi_y <- rpois(n, exp(x %*% beta))
offset <- rnorm(n)

# internal CV folds
nfolds <- 5
foldid <- sample(rep(seq(nfolds), length = n))

# default fit
fit1 <- rgam(x, y, foldid = foldid)
oldfit1 <- readRDS("saved_results/default.rds")
expect_equal(fit1, oldfit1)

# gaussian fit, non-linear features for only those from Step 1
fit2 <- rgam(x, y, foldid = foldid, init_nz = c())
oldfit2 <- readRDS("saved_results/gaussian_hierarchical.rds")
expect_equal(fit2, oldfit2)

# gaussian fit, different gamma and df parameters
fit3 <- rgam(x, y, foldid = foldid, gamma = 1, df = 6)
oldfit3 <- readRDS("saved_results/gaussian_diff_params.rds")
expect_equal(fit3, oldfit3)

# binomial family
fit4 <- rgam(x, bin_y, foldid = foldid, family = "binomial")
oldfit4 <- readRDS("saved_results/binomial.rds")
expect_equal(fit4, oldfit4)

# Poisson family
fit5 <- rgam(x, poi_y, foldid = foldid, family = "poisson")
oldfit5 <- readRDS("saved_results/poisson.rds")
expect_equal(fit5, oldfit5)

# Poisson with offset
fit6 <- rgam(x, poi_y, foldid = foldid, family = "poisson", offset = offset)
oldfit6 <- readRDS("saved_results/poisson_offset.rds")
expect_equal(fit6, oldfit6)