#' Compute ROC and other performance measures for binomial model
#'
#' Given a vector of true outcomes and a vector of predictions, returns a list
#' containing performance measures.
#'
#' We currently evaluate the performance measures at 100 quantiles of the
#' predicted values; this can be adjusted via the \code{N} option.
#'
#' @param ytest True test outcome: vector of 0s and 1s.
#' @param rit Predictions for the true outcome. Should be vector of continuous
#' variables between 0 and 1.
#' @param N Number of breakpoints where we evaluate the performance measures.
#' Default is 100.
#'
#' @return A list of performance measures and intermediate computations.
#' \item{sens}{Vector of sensitivity values.}
#' \item{spec}{Vector of specificity values.}
#' \item{ppv}{Vector of PPV values.}
#' \item{npv}{Vector of NPV values}
#' \item{area}{Area under ROC curve (AUC).}
#' \item{se}{Standard error for AUC.}
#' \item{cutp}{Cut points at which the performance measures were computed.}
#' \item{cutp.max}{Cut point which maximizes (sens + spec) / 2.}
#'
myroc <- function(ytest, rit, N = 100) {
    # compute sensitivity, specificity, NPV and PPV
    spec = sens = ppv = npv = rep(0, N)
    qq <- quantile(rit, probs = seq(0, 1, len = N))
    for(i in 1:N) {
        yhat <- 0 * (rit <= qq[i]) + 1 * (rit > qq[i])
        sens[i] <- sum((ytest == 1) & (yhat == 1)) / sum(ytest == 1)
        spec[i] <- sum((ytest == 0) & (yhat == 0)) / sum(ytest == 0)
        npv[i]  <- sum((ytest == 0) & (yhat == 0)) / sum(yhat == 0)
        ppv[i]  <- sum((ytest == 1) & (yhat == 1)) / sum(yhat == 1)
    }
    sens <- c(1, sens); spec <- c(0, spec)

    # compute AUC (via trapezoids)
    area <- 0
    for (i in 1:(length(sens)-1)) {
        area <- area + 0.5 * (spec[i+1] - spec[i]) * (sens[i+1] + sens[i])
    }

    a = area
    q1 = a / (2 - a)
    q2 = 2 * a * a / (1 + a)
    nn = sum(ytest == 0)
    na = sum(ytest == 1)
    se = sqrt((a*(1-a) + (na-1)*(q1-a*a) + (nn-1)*(q2-a*a)) / (na*nn))

    val = (spec+sens)/2
    cutp.max = qq[which.max(val)]
    return(list(sens=sens, spec=spec, ppv=ppv, npv=npv,
                area=area, se=se, cutp=qq, cutp.max=cutp.max))
}
