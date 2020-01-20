#' Get RGAM model component for one feature
#'
#' Returns the additive component of the RGAM model for a given feature at given
#' data points, i.e. f_j(X_j).
#'
#' @param object Fitted \code{rgam} object.
#' @param x Data for which we want the additive component. If \code{x} is a
#' matrix, it assumed that \code{X_j} is the \code{j}th column of this matrix.
#' If \code{x} is a vector, it is assumed to be \code{X_j} itself.
#' @param j The index of the original feature whose additive component we want.
#' @param index Index of lambda value for which plotting is desired. Default is
#' the last lambda value in \code{object$lambda}.
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
#' # get the additive component for the feature 6, x as matrix
#' f6 <- getf(fit, x, 6)  # last value of lambda
#' plot(x[, 6], f6)
#' f6 <- getf(fit, x, 6, index = 20)  # f1 at 20th value of lambda
#' plot(x[, 6], f6)
#'
#' # get the additive component for the feature 6, x as vector
#' new_x6 <- seq(-1, 1, length.out = 30)
#' new_f6 <- getf(fit, new_x6, 6)  # last value of lambda
#' plot(new_x6, new_f6)
#'
#' @export
getf <- function (object, x, j, index) {
    if (missing(index)) {
        index = length(object$lambda)
    }

    if (is.matrix(x)) {
        xval <- x[, j]
    } else {
        xval <- x
    }

    # get the linear part
    beta <- object$full_glmfit$beta[j , index]
    yval <- beta * xval

    # get the non-linear part
    if (j %in% object$init_nz) {
        l <- which(object$init_nz == j)

        temp <- object$spline_fit[[l]]
        fval <- predict(temp, scale(xval, object$mxf[j], object$sxf[j]))$y
        if (object$removeLin) {
            lm_coef <- object$lin_comp_fit[[l]]
            fval <- fval - lm_coef[1] - lm_coef[2] * xval
        }

        beta <- object$full_glmfit$beta[object$p + l, index]
        yval <- yval + beta * fval
    }

    yval
}
