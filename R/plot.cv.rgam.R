#' Plot the cross-validation curve produced by "cv.rgam" object
#'
#' Plots the cross-validation curve produced by a \code{cv.rgam} object, along
#' with upper and lower standard deviation curves, as a function of the \code{lambda}
#' values used. The plot also shows the number of non-zero features picked for
#' each value of lambda.
#'
#' A plot is produced and nothing is returned.
#'
#' @param x Fitted "\code{cv.rgam}" object.
#' @param sign.lambda Either plot against \code{log(lambda)} (default) or
#' \code{-log(lambda)} (if \code{sign.lambda = -1}).
#' @param ... Other graphical parameters to plot.
#'
#' @seealso \code{\link{rgam}} and \code{\link{cv.rgam}}.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#'
#' cvfit <- cv.rgam(x, y)
#' plot(cvfit)
#'
#' @export
plot.cv.rgam <- function(x, sign.lambda = 1, ...) {
    xlab <- "log(Lambda)"
    if (sign.lambda < 0) {
        xlab = paste("-", xlab, sep = "")
    }
    plot.args <- list(x = sign.lambda * log(x$lambda), y = x$cvm,
                      ylim = range(x$cvup, x$cvlo),
                      xlab = xlab, ylab = x$name, type = "n")
    new.args <- list(...)
    if (length(new.args)) {
        plot.args[names(new.args)] = new.args
    }
    do.call("plot", plot.args)
    error.bars(sign.lambda * log(x$lambda), x$cvup, x$cvlo,
               width = 0.01, col = "darkgrey")
    points(sign.lambda * log(x$lambda), x$cvm, pch = 20, col = "red")

    axis(side = 3, at = sign.lambda * log(x$lambda), labels = paste(x$nzero_feat),
         tick = FALSE, line = 0)

    abline(v = sign.lambda * log(x$lambda.min), lty = 3)
    abline(v = sign.lambda * log(x$lambda.1se), lty = 3)
    invisible()
}
