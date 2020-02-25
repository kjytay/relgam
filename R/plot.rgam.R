#' Make a plot of rgam model fit
#'
#' Produces plots of the estimated functions for specified variables at a given
#' value of lambda.
#'
#' A plot of the specified fitted functions is produced. Nothing is returned.
#'
#' @param x Fitted \code{rgam} object.
#' @param newx Matrix of values of each predictor at which to plot.
#' @param index Index of lambda value for which plotting is desired. Default is
#' the last lambda value in \code{x$lambda}.
#' @param which Which features to plot. Default is the first 4 or \code{nvars}
#' variables, whichever is smaller.
#' @param rugplot If \code{TRUE} (default), adds a rugplot showing the values of x
#' at the bottom of each fitted function plot.
#' @param grid_length The number of points to evaluate the estimated function at.
#' Default is 100.
#' @param names Vector of variable names of features in \code{which}. By default,
#' name of the \code{j}th variable is \code{xj}.
#' @param ... Optional graphical parameters to plot.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 12
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 3), rep(0, 9)), ncol = 1)
#' y <- x %*% beta + x[, 4]^2 + rnorm(n)
#' fit <- rgam(x, y)
#'
#' # default: print functions for first 4 variables
#' opar <- par(mfrow = c(2, 2))
#' plot(fit, newx = x, index = 20)
#' par(opar)
#'
#' # print for variables 5 to 8
#' opar <- par(mfrow = c(2, 2))
#' plot(fit, newx = x, index = 20, which = 5:8)
#' par(opar)
#'
#' @export
plot.rgam <- function(x, newx, index, which = NULL, rugplot = TRUE,
                     grid_length = 100, names, ...) {
    rgam.out = x
    x = newx
    p = ncol(x)
    if (missing(index)) {
        index = length(rgam.out$lambda)
    }

    if (is.null(which)) {
        warning(paste("Plotting first", min(p, 4), "variables by default"))
        which <- 1:(min(p, 4))
    }

    if (!missing(names)) {
        if (length(names) != length(which)) {
            warning("Length of names does not match length of which")
        }
        xlab <- names
        ylab <- paste0("f(", names, ")")
    } else {
        xlab <- paste0("x", which)
        ylab <- paste0("f(x", which, ")")
    }

    for (j in which) {
        # get xrange for the plot and the corresponding y-values
        xval <- seq(min(x[, j]), max(x[, j]), length.out = grid_length)
        yval <- getf(rgam.out, xval, j, index)

        # get color of the line
        if (j %in% rgam.out$nonlinfeat[[index]]) {
            colval <- "red"
        } else if (j %in% rgam.out$linfeat[[index]]) {
            colval <- "green"
        } else {
            colval <- "blue"
        }

        plot(x = xval, y = yval, type = "l", col = colval, lwd = 2,
             xlab = xlab[match(j, which)], ylab = ylab[match(j, which)], ...)
        if (rugplot) {
            rug(x[, j])
        }
    }
}
