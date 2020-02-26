context("Match v1.0 rgam()")
library(relgam)

compare_cvrgam <- function(oldfit, newfit) {
    for (key in names(oldfit)) {
        if (!(key %in% c("call", "glmfit"))) expect_equal(oldfit[[key]], newfit[[key]])
    }

    for (key in names(oldfit$glmfit)) {
        if (!(key %in% c("call", "spline_fit", "lin_comp_fit", "removeLin"))) {
            expect_equal(oldfit$glmfit[[key]], newfit$glmfit[[key]])
        }
    }
}

# set up data
set.seed(1)
n <- 100; p <- 20
x <- matrix(rnorm(n * p), n, p)
beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
y <- x %*% beta + rnorm(n)

# CV folds
nfolds <- 10
foldid <- sample(rep(seq(nfolds), length = n))

# specify number of folds
cvfit <- cv.rgam(x, y, foldid = foldid, keep = TRUE, verbose = FALSE)
oldfit <- readRDS("saved_results/cvfit.rds")
compare_cvrgam(oldfit, cvfit)

