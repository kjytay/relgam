#' Print a rgam object
#'
#' Print a summary of the rgam path at each step along the path.
#'
#' The call that produced the object x is printed, followed by a five-column
#' matrix with columns NonZero, Lin, NonLin, %Dev and Lambda. The first three
#' columns say how many nonzero, linear and nonlinear terms there are. %Dev is
#' the percent deviance explained (relative to the null deviance).
#'
#' @param x Fitted \code{rgam} object.
#' @param digits Significant digits in printout.
#' @param ... Additional print arguments.
#'
#' @return The matrix above is silently returned.
#'
#' @seealso \code{\link{rgam}}.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 12
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 3), rep(0, 9)), ncol = 1)
#' y <- x %*% beta + x[, 4]^2 + rnorm(n)
#' fit <- rgam(x, y)
#' print(fit)
#'
#' @export
print.rgam <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")

    out <- cbind(x$nzero_feat, x$nzero_lin, x$nzero_nonlin,
                 signif(x$full_glmfit$dev.ratio, digits),
                 signif(x$full_glmfit$lambda, digits))
    colnames(out) = c("NonZero", "Lin", "NonLin", "%Dev", "Lambda")
    print(out, ...)
    invisible(out)
}
