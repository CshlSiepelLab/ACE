#' Function get95Int
#' 
#' Returns the 95% confidence interval, based on the fisher information 
#' derived from an optimized parameter's hessian.  
#' Variation is calculated as var = 1/(-Hessian) = 1/(FI), as the inverse
#' of the Fisher information matrix is an estimator of the asymptotic
#' covariance matrix.  Assuming a normal distribution (questionable),
#' @param user_hessian Hessian matrix, returned from optim with hessian = T
#' @param user_par Value of the optimized parameter estimate; used to calculate
#' the confidence interval values.
#' @return Returns 95% confidence interval, with respect to \code{user_par}.
get95Int <- function(user_hessian, user_par) {
  # since we're maxing the log-likelihood, get the "observed info".
  fisher_info<-solve(-user_hessian)
  # hessian results are positive if at local min, negative at local max (intended).
  if (any(fisher_info < 0 )) {
    print('parameter not at local maximum, hessian non-negative.')
    return(NA)
  }
  # get standard error
  if (length(fisher_info) == 1) {
    prop_sigma <- sqrt(fisher_info)
  } else {
    prop_sigma<-sqrt(diag(fisher_info))
    prop_sigma<-diag(prop_sigma)
  }
  # get actual values for 95% confidence interval.
  upper<-user_par+1.96*prop_sigma
  lower<-user_par-1.96*prop_sigma
  interval<-list(upper=upper, lower=lower, se = 1.96*prop_sigma)
  return(interval)
}
# academically, the confidence interval is based on the 'standard error' σ of
# the estimates.  The '95%' CI is ± 1.96σ.
# So for a given value, what's the distribution of observations?
# You can fit a distribution to the probability gradient at the slope, 
# or the second derivative (hessian) around the estimated parameter value.
# 3 options for getting hessian:
# - optimHess(par, fn)
# - optim(..., hessian = T)