#' Fit reluctant generalized additive model
#'
#' Fits a reluctant generalized additive model (RGAM) for an entire regularization
#' path indexed by the parameter \code{lambda}. Fits linear, logistic, Poisson
#' and Cox regression models. RGAM is a three-step algorithm: Step 1 fits the
#' lasso and computes residuals, Step 2 constructs the non-linear features, and
#' Step 3 fits a lasso of the response on both the linear and non-linear features.
#'
#' If there are variables which the user definitely wants to compute non-linear
#' versions for in Step 2 of the algorithm, they can be specified as a vector for
#' the \code{init_nz} option. The algorithm will compute non-linear versions for
#' these features as well as the features suggested by Step 1 of the algorithm.
#'
#' If \code{standardize = TRUE}, the standard deviation of the linear and
#' non-linear features would be 1 and gamma respectively. If
#' \code{standardize = FALSE}, linear features will remain on their original
#' scale while non-linear features would have standard deviation gamma times
#' the mean standard deviation of the linear features.
#'
#' For \code{family="gaussian"}, \code{rgam} standardizes \code{y} to have unit
#' variance (using \code{1/n} rather than \code{1/(n-1)} formula).
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is
#' an observation vector.
#' @param y Response variable. Quantitative for \code{family = "gaussian"} or
#' \code{family = "poisson"} (non-negative counts). For \code{family="binomial"},
#' should be a numeric vector consisting of 0s and 1s. For \code{family="cox"},
#' y should be a two-column matrix with columns named 'time' and 'status'.
#' The latter is a binary variable, with '1' indicating death, and '0'
#' indicating right-censored.
#' @param lambda A user-supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence; supplying a value of
#' lambda overrides this.
#' @param lambda.min.ratio Smallest value for lambda as a fraction of the
#' largest lambda value. The default depends on the sample size nobs relative to
#' the number of variables nvars. If nobs > nvars, the default is 0.0001,
#' close to zero. If nobs < nvars, the default is 0.01.
#' @param standardize If \code{TRUE} (default), the columns of the input matrix
#' are standardized before the algorithm is run. See details section for more
#' information.
#' @param nl_maker This is a function that takes in one feature \code{x} and one 
#' response \code{r} and fits a model of \code{r} on \code{x}. It returns a list
#' of two items: (i) \code{f} which is a vector of fitted values, and (ii) a 
#' function \code{nl_predictor} which returns the fit for a given \code{newx} 
#' value. Additional arguments for \code{nl_maker} can be passed via \code{...}. 
#' The default is \code{relgam::makef}, which fits a smoothing spline with a 
#' user-specified degrees of freedom.
#' @param family Response type. Either \code{"gaussian"} (default) for linear
#' regression, \code{"binomial"} for logistic regression, \code{"poisson"} for
#' Poisson regression or \code{"cox"} for Cox regression.
#' @param offset A vector of length nobs. Useful for the "poisson" family (e.g.
#' log of exposure time), or for refining a model by starting at a current fit.
#' Default is NULL. If supplied, then values must also be supplied to the predict
#' function.
#' @param init_nz A vector specifying which features we must include when
#' computing the non-linear features. Default is to construct non-linear features
#' for all given features.
#' @param nfolds Number of folds for CV in Step 1 (default is 5). Although
#' \code{nfolds} can be as large as the sample size (leave-one-out CV), it is
#' not recommended for large datasets. Smallest value allowable is \code{nfolds = 3}.
#' @param foldid An optional vector of values between 1 and \code{nfolds}
#' identifying what fold each observation is in. If supplied, \code{nfolds} can
#' be missing.
#' @param gamma Scale factor for non-linear features (vs. original features), to
#' be between 0 and 1. Default is 0.8 if \code{init_nz = c()}, 0.6 otherwise.
#' @param parallel If TRUE, the \code{cv.glmnet()} call in Step 1 is
#' parallelized. Must register parallel before hand, such as doMC or others.
#' Default is FALSE.
#' @param verbose If \code{TRUE} (default), model-fitting is tracked with a
#' progress bar.
#' @param ... Any additional arguments to be the non-linear fitter in Step 2.
#'
#' @return An object of class \code{"rgam"}.
#' \item{full_glmfit}{The glmnet object resulting from Step 3: fitting a \code{glmnet}
#' model for the response against the linear & non-linear features.}
#' \item{nl_predictor}{List of functions used to get the non-linear features for 
#' new data. For internal use only.}
#' \item{init_nz}{Column indices for the features which we allow to have
#' non-linear relationship with the response.}
#' \item{step1_nz}{Indices of features which CV in Step 1 chose.}
#' \item{mxf}{Means of the features (both linear and non-linear).}
#' \item{sxf}{Scale factors of the features (both linear and non-linear).}
#' \item{feat}{Column indices of the non-zero features for each value of
#' \code{lambda}.}
#' \item{linfeat}{Column indices of the non-zero linear components for each value of
#' \code{lambda}.}
#' \item{nonlinfeat}{Column indices of the non-zero non-linear components for each value
#' of \code{lambda}.}
#' \item{nzero_feat}{The number of non-zero features for each value of
#' \code{lambda}.}
#' \item{nzero_lin}{The number of non-zero linear components for each value of
#' \code{lambda}.}
#' \item{nzero_nonlin}{The number of non-zero non-linear components for each value
#' of \code{lambda}.}
#' \item{lambda}{The actual sequence of \code{lambda} values used.}
#' \item{p}{The number of features in the original data matrix.}
#' \item{family}{Response type.}
#' \item{call}{The call that produced this object.}
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#'
#' fit <- rgam(x, y)
#'
#' # construct non-linear features for only those selected by Step 1
#' fit <- rgam(x, y, init_nz = c())
#'
#' # specify scale factor gamma and degrees of freedom
#' fit <- rgam(x, y, gamma = 1, df = 6)
#'
#' # binomial family
#' bin_y <- ifelse(y > 0, 1, 0)
#' fit2 <- rgam(x, bin_y, family = "binomial")
#'
#' # Poisson family
#' poi_y <- rpois(n, exp(x %*% beta))
#' fit3 <- rgam(x, poi_y, family = "poisson")
#' # Poisson with offset
#' offset <- rnorm(n)
#' fit3 <- rgam(x, poi_y, family = "poisson", offset = offset)
#'
#' @import utils
#' @export
rgam <- function(x, y, lambda = NULL, lambda.min.ratio = ifelse(nrow(x) < ncol(x),
                0.01, 1e-04), standardize = TRUE, nl_maker = relgam:::makef,
                family = c("gaussian","binomial", "poisson", "cox"), offset = NULL,
                init_nz, nfolds = 5, foldid = NULL,
                gamma, parallel = FALSE, verbose = TRUE, ...) {
    this.call <- match.call()

    n <- nrow(x); p <- ncol(x)
    if (is.null(offset)) {
        offset <- y * 0
        is.offset = FALSE
    } else {
        is.offset = TRUE
    }

    family <- match.arg(family)
    if (family == "binomial" && any(!(unique(y) %in% c(0, 1)))) {
        stop("If family is binomial, y can only contain 0s and 1s")
    }
    if (family == "poisson" && any(y < 0)) {
        stop("If family is Poisson, y values cannot be negative")
    }

    # specify init_nz and gamma
    if (missing(init_nz)) {
        message("init_nz not specified: setting to default (all features)")
        init_nz <- 1:p
    }
    if (missing(gamma)) {
        if (length(init_nz) == 0) {
            message("using default value of gamma for RGAM_SEL: 0.8")
            gamma <- 0.8
        } else {
            message("using default value of gamma for RGAM: 0.6")
            gamma <- 0.6
        }
    }

    # center (and scale) columns of x and y
    mx <- colMeans(x)
    if (standardize) {
        sx <- apply(x, 2, sd) * sqrt((n-1) / n)
    } else {
        sx <- rep(1, ncol(x))
    }
    x <- scale(x, mx, sx)
    if (family == "gaussian") {
        my <- mean(y); sy <- sd(y) * sqrt((n-1) / n)
        y <- scale(y, my, sy)
    }

    # remove features with sx == 0 from consideration
    exclude <- which(sx == 0)
    x[, exclude] <- 0

    # if fold IDs not given, randomly create them
    # if given, count the number of folds
    if (is.null(foldid)) {
        foldid <- sample(rep(seq(nfolds), length = n))
    } else {
        nfolds <- length(unique(foldid))
    }
    if (nfolds < 3) {
        stop("nfolds must be 3 or larger; nfolds = 5 recommended")
    }

    # start progress bar
    if (verbose)  pb <- utils::txtProgressBar(min = 0, max = 3, style = 3)

    # Step 1: Get residuals after "best linear fit"
    cv <- glmnet::cv.glmnet(x, y, standardize = FALSE, foldid = foldid,
                            family = family, parallel = parallel, exclude = exclude,
                            lambda.min.ratio = lambda.min.ratio, offset = offset)

    if (family %in% c("gaussian", "binomial", "poisson")) {
        yhat0 <- predict(cv, x, s = cv$lambda.min, type = "response", newoffset = offset)
        r <- y - yhat0
    } else if (family == "cox") {
        yhat0 <- predict(cv, x, s = cv$lambda.min, type = "link", newoffset = offset)
        r <- coxgrad(yhat0,y[,1],y[,2])
    }

    # get index of features which have non-zero coefficients for step 2
    step1_nz <- predict(cv, s = cv$lambda.min, type="nonzero")
    # to handle special case of nz_feat being empty
    if (is.null(dim(step1_nz))) {
        step1_nz <- numeric()
    } else {
        step1_nz <- step1_nz[, 1]
    }

    init_nz <- sort(unique(c(step1_nz, init_nz)))

    # update progress bar
    if (verbose) utils::setTxtProgressBar(pb, 1)

    # Step 2: make the non-linear features x
    nl_predictor <- list()
    f <- matrix(NA, n, length(init_nz))
    for (j in 1:length(init_nz)) {
        out <- nl_maker(x[, init_nz[j]], r, ...)
        f[, j] <- out$f
        nl_predictor[[j]] <- out$nl_predictor
    }

    # standardize non-linear features
    mf <- colMeans(f)
    if (standardize) {
        sf <- apply(f, 2, sd) * sqrt((n-1) / n) / gamma
    } else {
        sf <- (apply(f, 2, sd) / mean(apply(x, 2, sd))) / gamma
    }
    f <- scale(f, mf, sf)

    # update the mean and sd vectors, and exclude
    mxf <- c(mx, mf)
    sxf <- c(sx, sf)
    exclude <- c(exclude, p + which(sf == 0))
    f[, which(sf == 0)] <- 0

    # update progress bar
    if (verbose) utils::setTxtProgressBar(pb, 2)

    # Step 3: fit lasso with both linear and non-linear features
    full_glmfit <- glmnet::glmnet(cbind(x, f), y, lambda = lambda, standardize = FALSE,
                                  family = family, exclude = exclude,
                                  lambda.min.ratio = lambda.min.ratio, offset = offset)

    # rescale intercept values and beta values
    beta <- matrix(full_glmfit$beta, nrow = p + length(init_nz))
    if (family == "gaussian") {
        full_glmfit$a0 <- full_glmfit$a0 * sy - colSums(beta * mxf / sxf) * sy + my
        full_glmfit$beta <- t(scale(t(beta), F, sxf)) * sy
    } else {
        full_glmfit$a0 <- full_glmfit$a0 - colSums(beta * mxf / sxf)
        full_glmfit$beta <- t(scale(t(beta), F, sxf))
    }

    # get no. of nonzeros for linear and non-linear components
    nz <- predict(full_glmfit, type = "nonzero")
    linfeat <- lapply(nz, function(x) x[x <= p])
    nonlinfeat <- lapply(nz, function(x) init_nz[x[x > p] - p])
    nzero_lin <- unlist(lapply(linfeat, length))
    nzero_nonlin <- unlist(lapply(nonlinfeat, length))

    # figure out no. of nonzero features (whether linear or non-linear)
    feat <- lapply(1:length(nz),
                   function(i) unique(c(linfeat[[i]], nonlinfeat[[i]])))
    nzero_feat <- unlist(lapply(feat, length))

    full_glmfit$offset <- is.offset

    # update progress bar
    if (verbose) utils::setTxtProgressBar(pb, 3)

    out <- list(full_glmfit = full_glmfit, nl_predictor = nl_predictor,
                init_nz = init_nz, step1_nz = step1_nz, mxf = mxf, sxf = sxf,
                feat = feat, linfeat = linfeat, nonlinfeat = nonlinfeat,
                nzero_feat = nzero_feat, nzero_lin = nzero_lin,
                nzero_nonlin = nzero_nonlin, lambda = full_glmfit$lambda,
                p = p, family = family, call = this.call)
    class(out) <- "rgam"
    return(out)
}
