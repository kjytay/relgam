context("Match v1.0 rgam()")
library(relgam)

compare_rgam <- function(oldfit, newfit, test_name) {
    for (key in names(oldfit)) {
        if (!(key %in% c("call", "spline_fit", "lin_comp_fit", "removeLin"))) {
            expect_equal(oldfit[[key]], newfit[[key]],
                                                label = test_name)
        }
    }
}

compare_predictions <- function(oldpred, newfit, x, test_name, ...) {
    newpred <- predict(newfit, x, ...)
    expect_equal(newpred, oldpred, label = paste(test_name, "prediction"))
}

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
fit1 <- rgam(x, y, foldid = foldid, verbose = FALSE)
oldfit1 <- readRDS("saved_results/default.rds")
oldpred1 <- readRDS("saved_results/pred_default.rds")
compare_rgam(oldfit1, fit1, "default")
compare_predictions(oldpred1, fit1, x, "default")

# gaussian fit, non-linear features for only those from Step 1
fit2 <- rgam(x, y, foldid = foldid, init_nz = c(), verbose = FALSE)
oldfit2 <- readRDS("saved_results/gaussian_hierarchical.rds")
oldpred2 <- readRDS("saved_results/pred_gaussian_hierarchical.rds")
compare_rgam(oldfit2, fit2, "gaussian_hierarchical")
compare_predictions(oldpred2, fit2, x, "gaussian_hierarchical")

# gaussian fit, different gamma and df parameters
fit3 <- rgam(x, y, foldid = foldid, gamma = 1, df = 6, verbose = FALSE)
oldfit3 <- readRDS("saved_results/gaussian_diff_params.rds")
oldpred3 <- readRDS("saved_results/pred_gaussian_diff_params.rds")
compare_rgam(oldfit3, fit3, "gaussian_diff_params")
compare_predictions(oldpred3, fit3, x, "gaussian_diff_params")

# binomial family
fit4 <- rgam(x, bin_y, foldid = foldid, family = "binomial", verbose = FALSE)
oldfit4 <- readRDS("saved_results/binomial.rds")
oldpred4 <- readRDS("saved_results/pred_binomial.rds")
compare_rgam(oldfit4, fit4, "binomial")
compare_predictions(oldpred4, fit4, x, "binomial")

# Poisson family
fit5 <- rgam(x, poi_y, foldid = foldid, family = "poisson", verbose = FALSE)
oldfit5 <- readRDS("saved_results/poisson.rds")
oldpred5 <- readRDS("saved_results/pred_poisson.rds")
compare_rgam(oldfit5, fit5, "poisson")
compare_predictions(oldpred5, fit5, x, "poisson")

# Poisson with offset
fit6 <- rgam(x, poi_y, foldid = foldid, family = "poisson", offset = offset, 
             verbose = FALSE)
oldfit6 <- readRDS("saved_results/poisson_offset.rds")
oldpred6 <- readRDS("saved_results/pred_poisson_offset.rds")
compare_rgam(oldfit6, fit6, "poisson_offset")
compare_predictions(oldpred6, fit6, x, "poisson_offset", newoffset = offset)