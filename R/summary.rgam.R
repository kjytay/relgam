#' rgam summary routine
#'
#' Makes a two-panel plot of the rgam object showing coefficient paths.
#'
#' A two panel plot is produced, that summarizes the linear components and
#' the nonlinear components, as a function of lambda. For the linear components,
#' it is the coefficient for each variable. For the nonlinear components, it is
#' the coefficient of the non-linear variable. Nothing is returned.
#'
#' @param object Fitted \code{rgam} object.
#' @param label If \code{TRUE}, annotate the plot with variable labels. Default
#' is \code{FALSE}.
#' @param index The indices of the lambda hyperparameter which we want the plot
#' for. The default is to plot for the entire lambda path.
#' @param which Which values to plot. Default is all variables.
#' @param ... Additional arguments to summary.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#'
#' fit <- rgam(x, y)
#' opar <- par(mfrow = c(1, 2))
#' summary(fit)
#' par(opar)
#'
#' # with labels, just variables 1 to 5
#' opar <- par(mfrow = c(1, 2))
#' summary(fit, label = TRUE, which = 1:5)
#' par(opar)
#'
#' # as above, but just the first 30 values of lambda
#' opar <- par(mfrow = c(1, 2))
#' summary(fit, label = TRUE, which = 1:5, index = 1:30)
#' par(opar)
#'
#' @export
summary.rgam <- function (object, label = FALSE, index = NULL, which = NULL, ...)
{
    maxvar <- nrow(object$full_glmfit$beta)
    p <- maxvar - length(object$init_nz)
    init_nz <- object$init_nz
    colours <- rainbow(n = p)

    # get the relevant lambda
    if (is.null(index)) {
        lambda <- object$lambda
        index <- 1:length(lambda)
    } else {
        lambda <- object$lambda[index]
    }

    # compute which
    if (is.null(which)) {
        which <- 1:p
    } else {
        which <- sort(which)
    }

    # internal function for drawing line for one feature
    nzlines = function(lambda, alpha, ...) {
        if (any(abs(alpha) > 0)) {
            num_lambda <- length(lambda)
            start = max(1, min(seq(num_lambda)[abs(alpha) > 0]) - 1)
            whichnz = seq(from = start, to = num_lambda)
            if (length(whichnz) > 1) {
                lines(lambda[whichnz], alpha[whichnz], ...)
            }
        }
        invisible()
    }

    # plot for linear coefficients
    lin_beta <- object$full_glmfit$beta[1:p, index, drop = FALSE]
    plot(0, type = "n", xlab = expression(lambda), ylab = "Coefficient",
         xlim = c(max(lambda) * 1.1, min(lambda) * (0.9 - 0.1 * label)),
         ylim = range(lin_beta), main = "Linear Components", log = "x")
    abline(h = 0, lty = 3)
    for (j in 1:p) {
        if (j %in% which) {
            nzlines(lambda, lin_beta[j, ], col = colours[j], lwd = 2,
                    type = "l", pch = "", cex = 0)
        }
    }
    if (label) {
        text(rep(min(lambda) * 0.85, length(which)), lin_beta[which, length(lambda)],
             labels = which, col = colours, cex = 0.6)
    }

    # plot for non-linear coefficients
    nonlin_beta <- object$full_glmfit$beta[(p+1):maxvar, index, drop = FALSE]
    plot(0, type = "n", xlab = expression(lambda), ylab = "Coefficient",
         xlim = c(max(lambda) * 1.1, min(lambda) * (0.9 - 0.1 * label)),
         ylim = range(nonlin_beta), main = "Non-Linear Components", log = "x")
    abline(h = 0, lty = 3)
    j_list <- c()
    label_list <- c()
    if (maxvar > p) {
        for (j in 1:(maxvar - p)) {
            if (init_nz[j] %in% which) {
                j_list <- c(j_list, j)
                label_list <- c(label_list, init_nz[j])
                nzlines(lambda, nonlin_beta[j, ], col = colours[init_nz[j]], lwd = 2,
                        type = "l", pch = "", cex = 0)
            }
        }
    }
    if (label) {
        text(rep(min(lambda) * 0.85, length(j_list)), nonlin_beta[j_list, length(lambda)],
             labels = label_list, col = colours[label_list], cex = 0.6)
    }
}
