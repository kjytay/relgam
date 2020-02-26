#' Make non-linear features
#'
#' Internal function for making non-linear features.
#'
#' @param x Input vector of length \code{nobs}.
#' @param r Vector of residuals.
#' @param df Degrees of freedom for the fit. Default is 4.
#' @param tol A tolerance for same-ness or uniqueness of the x values. To be
#' passed to the \code{smooth.spline()} function. Default is \code{0.01}.
#' @param removeLin If \code{TRUE} (default), removes the linear component from
#' the newly created non-linear features.
#'
#' @return A list:
#' \item{f}{Non-linear feature associated with \code{x}.}
#' \item{nl_predictor}{A function which, when given new data \code{newx}, returns 
#' the value of the non-linear predictor at those \code{x} values. Needed for 
#' creating the non-linear features for new data.}
makef <- function(x, r, df = 4, tol = 0.01, removeLin = T) {
    n <- nrow(x)

    temp <- smooth.spline(x, r, df = df, tol = tol)
    f <- predict(temp, x)$y

    # remove linear component if asked to
    lin_comp_fit <- NULL
    if (removeLin) {
        lm_model <- lsfit(x, f)
        f <- lm_model$res
        lin_comp_fit <- lm_model$coefficients
    }
    
    # create predictor function to help us get the non-linear function for 
    # future x values
    if (removeLin) {
        nl_predictor <- function(xnew) {
            fnew <- predict(temp, xnew)$y
            fnew - lin_comp_fit[1] - lin_comp_fit[2] * xnew
        }
    } else {
        nl_predictor <- function(xnew) {
            predict(temp, xnew)$y
        }
    }
    
    return(list(f = f, nl_predictor = nl_predictor))
}
