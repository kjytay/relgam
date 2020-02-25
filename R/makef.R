#' Make non-linear features
#'
#' Internal function for making non-linear features.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is
#' an observation vector.
#' @param r Vector of residuals.
#' @param df Degrees of freedom for the fit. Default is 4.
#' @param tol A tolerance for same-ness or uniqueness of the x values. To be
#' passed to the \code{smooth.spline()} function. Default is \code{0.01}.
#' @param removeLin If \code{TRUE} (default), removes the linear component from
#' the newly created non-linear features.
#'
#' @return A list:
#' \item{f}{Non-linear features associated with the features in \code{x}.}
#' \item{spline_fit}{A list of the spline fits of the residual against each
#' feature. Useful for creating the non-linear features for new data.}
#' \item{lin_comp_fit}{If \code{removeLin = TRUE}, a list of coefficients for
#' simple linear regression of non-linear feature on original feature. Useful
#' for creating the non-linear features for new data.}
makef <- function(x, r, df = 4, tol = 0.01, removeLin = T) {
    n <- nrow(x); p <- ncol(x)

    f <- matrix(NA, n, p)
    spline_fit <- list()
    lin_comp_fit <- list()
    if (p > 0) {
        for (j in 1:p) {
            temp <- smooth.spline(x[, j], r, df = df, tol = tol)
            f[,j] <- predict(temp, x[, j])$y

            # remove linear component if asked to
            if(removeLin) {
                lm_model <- lsfit(x[,j], f[,j])
                f[,j] <- lm_model$res
                lin_comp_fit[[j]] <- lm_model$coefficients
            }

            spline_fit[[j]] <- temp
        }
    }

    return(list(f = f, spline_fit = spline_fit, lin_comp_fit = lin_comp_fit))
}
