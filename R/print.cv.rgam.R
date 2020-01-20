#' Print a cross-validated rgam object
#'
#' Print a summary of the results of cross-validation for a RGAM model.
#'
#' The call that produced the object x is printed, followed by some information
#' on the performance for \code{lambda.min} and \code{lambda.1se}.
#'
#' @param x Fitted \code{rgam} object.
#' @param digits Significant digits in printout.
#' @param ... Additional print arguments.
#'
#' @seealso \code{\link{cv.rgam}}, \code{\link{print.rgam}}.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#' cvfit <- cv.rgam(x, y)
#' print(cvfit)
#'
#' @export
print.cv.rgam <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")
    cat("Measure:", x$name, "\n\n")

    # get indices for optimal lambda values
    optlams = c(x$lambda.min, x$lambda.1se)
    which = match(optlams, x$lambda)
    out <- cbind(signif(optlams, digits),
                 signif(x$cvm[which], digits), signif(x$cvse[which], digits),
                 x$nzero_feat[which], x$nzero_lin[which], x$nzero_nonlin[which])
    dimnames(out) = list(c("min", "1se"), c("Lambda", "Measure",
                                            "SE", "Nonzero", "Lin", "NonLin"))
    print(out, ...)
}
