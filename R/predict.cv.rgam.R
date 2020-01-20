#' Make predictions from a "cv.rgam" object
#'
#' This function returns the predictions for a new data matrix from a
#' cross-validated \code{rgam} model by using the stored "\code{glmfit}"
#' object and the optimal value chosen for \code{lambda}.
#'
#' This function makes it easier to use the results of cross-validation to make
#' a prediction.
#'
#' @param object Fitted "\code{cv.rgam}" object.
#' @param xnew Matrix of new values for \code{x} at which predictions are to
#' be made.
#' @param s Value of the penalty parameter \code{lambda} at which predictions are
#' required. Default is the value \code{s="lambda.1se"} stored in the CV
#' \code{fit}. Alternatively, \code{s="lambda.min"} can be used. If \code{s} is
#' numeric, it is taken as the value(s) of lambda to be used.
#' @param ... Other arguments to be passed to \code{predict.rgam())}.
#'
#' @return Predictions which the cross-validated model makes for \code{xnew} at
#' the optimal value of \code{lambda}. Note that the default is the "lambda.1se"
#' for lambda, to make this function consistent with \code{cv.glmnet} in the
#' \code{glmnet} package.
#'
#' The output depends on the \code{...} argument which is passed on to the predict
#' method for \code{rgam} objects.
#'
#' @seealso \code{\link{cv.rgam}} and \code{\link{predict.rgam}}.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)
#' cvfit <- cv.rgam(x, y)
#'
#' # predictions at the lambda.1se value
#' predict(cvfit, xnew = x[1:5, ])
#'
#' # predictions at the lambda.min value
#' predict(cvfit, xnew = x[1:5, ], s = "lambda.min")
#'
#' # predictions at specific lambda value
#' predict(cvfit, xnew = x[1:5, ], s = 0.1)
#'
#' # probability predictions for binomial family
#' bin_y <- ifelse(y > 0, 1, 0)
#' cvfit2 <- cv.rgam(x, bin_y, family = "binomial")
#' predict(cvfit2, xnew = x[1:5, ], type = "response", s = "lambda.min")
#'
#' @export
predict.cv.rgam <- function(object, xnew, s = c("lambda.1se", "lambda.min"),
                           ...) {
    if (is.numeric(s))
        lambda = s
    else if (is.character(s)) {
        s = match.arg(s)
        lambda = object[[s]]
    }
    else stop("Invalid form for s")
    predict(object$glmfit, xnew, s = lambda, ...)
}
