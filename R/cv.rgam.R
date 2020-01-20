#' Cross-validation for reluctant generalized additive model (rgam)
#'
#' Does \code{k}-fold cross-validation for \code{rgam}.
#'
#' The function runs \code{rgam} nfolds+1 times; the first to get the lambda
#' sequence, and then the remainder to compute the fit with each of the folds
#' omitted. The error is accumulated, and the average error and standard
#' deviation over the folds is computed.
#'
#' Note that \code{cv.rgam} only does cross-validation for lambda but not for
#' the degrees of freedom hyperparameter.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is
#' an observation vector.
#' @param y Response \code{y} as in \code{rgam}.
#' @param lambda A user-supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence; supplying a value of
#' lambda overrides this.
#' @param family Response type. Either \code{"gaussian"} (default) for linear
#' regression, \code{"binomial"} for logistic regression, \code{"poisson"} for
#' Poisson regression or \code{"cox"} for Cox regression.
#' @param offset Offset vector as in \code{rgam}.
#' @param init_nz A vector specifying which features we must include when
#' computing the non-linear features. Default is to construct non-linear
#' features for all given features.
#' @param gamma Scale factor for non-linear features (vs. original features),
#' to be between 0 and 1. Default is 0.8 if \code{init_nz = c()}, 0.6 otherwise.
#' @param nfolds Number of folds for CV (default is 10). Although \code{nfolds}
#' can be as large as the sample size (leave-one-out CV), it is not recommended
#' for large datasets. Smallest value allowable is \code{nfolds = 4}.
#' @param foldid An optional vector of values between 1 and \code{nfolds}
#' identifying what fold each observation is in. If supplied, \code{nfolds} can
#' be missing.
#' @param keep If \code{keep = TRUE}, a prevalidated array is returned
#' containing fitted values for each observation at each value of lambda. This
#' means these fits are computed with this observation and the rest of its fold
#' omitted. Default is \code{FALSE}.
#' @param parallel If TRUE, use parallel foreach to fit each fold. Must
#' register parallel before hand, such as doMC or others. Note that this also
#' passes \code{parallel = TRUE} to the \code{rgam()} call within each fold.
#' Default is FALSE.
#' @param verbose Print information as model is being fit? Default is
#' \code{TRUE}.
#' @param ... Other arguments that can be passed to \code{rgam}.
#'
#' @return An object of class \code{"cv.rgam"}.
#' \item{glmfit}{A fitted \code{rgam} object for the full data.}
#' \item{lambda}{The values of \code{lambda} used in the fits.}
#' \item{nzero_feat}{The number of non-zero features for the model \code{glmfit}.}
#' \item{nzero_lin}{The number of non-zero linear components for the model
#' \code{glmfit}.}
#' \item{nzero_nonlin}{The number of non-zero non-linear components for the
#' model \code{glmfit}.}
#' \item{fit.preval}{If \code{keep=TRUE}, this is the array of prevalidated
#' fits.}
#' \item{cvm}{The mean cross-validated error: a vector of length
#' \code{length(lambda)}.}
#' \item{cvse}{Estimate of standard error of \code{cvm}.}
#' \item{cvlo}{Lower curve = \code{cvm - cvsd}.}
#' \item{cvup}{Upper curve = \code{cvm + cvsd}.}
#' \item{lambda.min}{The value of \code{lambda} that gives minimum
#'   \code{cvm}.}
#' \item{lambda.1se}{The largest value of \code{lambda} such that the CV
#'   error is within one standard error of the minimum.}
#' \item{foldid}{If \code{keep=TRUE}, the fold assignments used.}
#' \item{name}{Name of error measurement used for CV.}
#' \item{call}{The call that produced this object.}
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#'
#' cvfit <- cv.rgam(x, y)
#'
#' # specify number of folds
#' cvfit <- cv.rgam(x, y, nfolds = 5)
#'
#' @import foreach
#' @export
cv.rgam <- function(x, y, lambda = NULL, family = c("gaussian","binomial",
                   "poisson", "cox"), offset = NULL, init_nz, gamma, nfolds = 10,
                   foldid = NULL, keep = FALSE, parallel = FALSE,
                   verbose = TRUE, ...) {
    this.call <- match.call()
    family <- match.arg(family)
    n <- nrow(x); p <- ncol(x)
    if (is.null(offset)) {
        offset <- y * 0
        is.offset = FALSE
    } else {
        is.offset = TRUE
    }

    # if fold IDs not given, randomly create them
    # if given, count the number of folds
    if (is.null(foldid)) {
        foldid <- sample(rep(seq(nfolds), length = n))
    } else {
        nfolds <- length(unique(foldid))
    }
    if (nfolds <= 3) {
        stop("nfolds must be bigger than 3; nfolds = 10 recommended")
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

    # get rgam fit for all of the data
    if (verbose) cat("Get initial fit", fill = TRUE)
    fit0 <- rgam(x, y, lambda = lambda, family = family, foldid = foldid,
                offset = offset, init_nz = init_nz, gamma = gamma, verbose = verbose, ...)
    fit0$full_glmfit$offset <- is.offset
    if (verbose) cat("\nInitial fit done", fill = TRUE)

    # fit rgam on folds
    # note: we run each rgam with K-1 folds, reusing the foldid vector
    fits <- vector("list", nfolds)
    if (parallel) {
        fits = foreach(ii = 1:nfolds, .packages = c("rgam")) %dopar% {
            oo <- foldid == ii
            xc <- x[!oo, , drop = FALSE]
            if (family == "cox") {
                yy <- y[!oo, , drop = FALSE]
            } else {
                yy <- y[!oo]
            }
            temp_foldid <- foldid[!oo]
            if (ii != nfolds) {
                temp_foldid[temp_foldid == nfolds] <- ii
            }
            rgam(xc, yy, lambda = fit0$lambda, family = family,
                offset = offset[!oo], foldid = temp_foldid,
                init_nz = init_nz, gamma = gamma, ...)
        }
    } else {
        for (ii in 1:nfolds) {
            if (verbose) {
                cat(c("\nFold = ", ii), fill = TRUE)
            }
            oo <- foldid == ii
            xc <- x[!oo, , drop = FALSE]
            if (family == "cox") {
                yy <- y[!oo, , drop = FALSE]
            } else {
                yy <- y[!oo]
            }
            temp_foldid <- foldid[!oo]
            if (ii != nfolds) {
                temp_foldid[temp_foldid == nfolds] <- ii
            }
            fits[[ii]] <- rgam(xc, yy, lambda = fit0$lambda, family = family,
                              offset = offset[!oo], foldid = temp_foldid,
                              init_nz = init_nz, gamma = gamma, verbose = verbose, ...)
        }
    }

    # get predictions
    yhat <- matrix(NA, n, length(fit0$lambda))
    for (ii in 1:nfolds) {
        oo <- foldid == ii
        out <- predict(fits[[ii]], x[oo, , drop = F], newoffset = offset[oo])
        yhat[oo, 1:ncol(out)] <- out
    }

    # if keep = TRUE, keep the pre-validated fits
    yhat.preval <- NULL
    foldid_copy <- NULL
    if (keep) {
        yhat.preval <- yhat
        foldid_copy <- foldid
    }

    # compute CV error
    if (family == "cox") {
        name <- "Partial Likelihood Deviance"
        err <- matrix(0, nrow = nfolds, ncol = length(fit0$lambda))
        for (ii in 1:nfolds) {
            oo <- foldid == ii
            yhat <- predict(fits[[ii]], x, newoffset = offset)
            err[ii, ] = glmnet::coxnet.deviance(pred = yhat, y = y) -
                glmnet::coxnet.deviance(pred = yhat[!oo, ], y = y[!oo,])
        }
        weights <- as.vector(tapply(y[, "status"], foldid, sum))
        err <- err / weights
        cvm <- apply(err, 2, weighted.mean, w = weights, na.rm = TRUE)
        cvse <- sqrt(apply(scale(err, cvm, FALSE)^2, 2, weighted.mean,
                   w = weights, na.rm = TRUE) / (nfolds - 1))
    } else {
        # assumption: family is gaussian or binomial
        name <- switch(family,
                       "gaussian" = "Mean-Squared Error",
                       "binomial" = "Deviance",
                       "poisson" = "Poisson Deviance")
        errfun <- switch(family,
                         "gaussian" = msefun,
                         "binomial" = binfun,
                         "poisson" = poifun)
        if (family == "binomial") {
            yhat <- 1 / (1 + exp(-yhat))
        }
        ym <- array(y, dim(yhat))
        err <- errfun(yhat, ym)
        cvm <- apply(err, 2, mean, na.rm = T)
        err2 <- matrix(NA, nfolds, length(fit0$lambda))
        for (ii in 1:nfolds) {
            err2[ii, ] <- colMeans(err[foldid == ii, ])
        }
        cvse <- sqrt(apply(err2, 2, var, na.rm = T) / nfolds)
    }
    cvlo <- cvm - cvse
    cvup <- cvm + cvse

    # get best and 1se lambda
    imin <- which.min(cvm)
    lambda.min <- fit0$lambda[imin]
    imin.1se <- which(cvm < cvm[imin] + cvse[imin])[1]
    lambda.1se <- fit0$lambda[imin.1se]

    out <- list(glmfit=fit0, lambda=fit0$lambda, nzero_feat = fit0$nzero_feat,
                nzero_lin=fit0$nzero_lin, nzero_nonlin=fit0$nzero_nonlin,
                fit.preval=yhat.preval, cvm=cvm, cvse=cvse, cvlo=cvlo,
                cvup=cvup, lambda.min=lambda.min, lambda.1se=lambda.1se,
                foldid=foldid_copy, name=name, call=this.call)
    class(out) <- "cv.rgam"
    return(out)
}
