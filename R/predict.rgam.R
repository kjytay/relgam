#' Make predictions from a "rgam" object
#'
#' This function returns the predictions from a "\code{rgam}" object
#' for a new data matrix.
#'
#' @param object Fitted "\code{rgam}" object.
#' @param xnew Matrix of new values for \code{x} at which predictions are to
#' be made.
#' @param ... Any other arguments to be passed to \code{predict.glmnet()}.
#'
#' @return Predictions of which the model \code{object} makes at
#' \code{xnew}. The type of predictions depends on whether a \code{type} argument
#' is passed. By default it givs the linear predictors for the regression model.
#'
#' If an offset is used in the fit, then one must be supplied via the
#' \code{newoffset} option.
#'
#' @seealso \code{\link{rgam}}.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#' fit <- rgam(x, y)
#'
#' # predict for full lambda path
#' predict(fit, xnew = x[1:5, ])
#'
#' # predict for specific lambda values
#' predict(fit, xnew = x[1:5, ], s = 0.1)
#'
#' # predictions for binomial family
#' bin_y <- ifelse(y > 0, 1, 0)
#' fit2 <- rgam(x, bin_y, family = "binomial")
#' # linear predictors
#' predict(fit2, xnew = x[1:5, ], s = 0.05)
#' # probabilities
#' predict(fit2, xnew = x[1:5, ], type = "response", s = 0.05)
#'
#' @export
predict.rgam <- function(object, xnew, ...) {
    # make the non-linear features for xnew
    xnew_nonlin = xnew[, object$init_nz, drop = F]
    xnew_nonlin <- scale(xnew_nonlin, object$mxf[object$init_nz],
                         object$sxf[object$init_nz])

    fnew <- matrix(NA, nrow(xnew_nonlin), ncol(xnew_nonlin))
    if (ncol(xnew_nonlin) > 0) {
        for (j in 1:ncol(xnew_nonlin)) {
            temp <- object$spline_fit[[j]]
            fnew[, j] <- predict(temp, xnew_nonlin[, j])$y
            if (object$removeLin) {
                lm_coef <- object$lin_comp_fit[[j]]
                fnew[, j] <- fnew[, j] - lm_coef[1] - lm_coef[2] * xnew_nonlin[, j]
            }
        }
    }

    # make predictions for (xnew, fnew)
    out <- predict(object$full_glmfit, cbind(xnew, fnew), ...)

    return(out)
}
