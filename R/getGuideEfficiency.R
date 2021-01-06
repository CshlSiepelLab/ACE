#' Function getGuideEfficiency
#' 
#' Converts normalized feature matrix to by-guide percent efficiency.
#' Adds constant term not contained in guide features.
#' @param guide_matrix Data.table of feature values for each guide.
#' @param feature_weight Vector of feature weights for sigmoid transform.
#' @return Returns vector of by-guide percent efficiency.
getGuideEfficiency <- function(guide_matrix, feature_weight) {
  if (!is.matrix(guide_matrix)) stop('must submit numeric matrix of guide features.')
  full_guide_matrix <- cbind(rep(1, nrow(guide_matrix)), guide_matrix)
  if (ncol(full_guide_matrix) != length(feature_weight)) {
    stop('dimensions of feature weights and guide feature data incompatible.')
  }
  fwm <- matrix(feature_weight, ncol=1)
  guide_efficiency <-  exp((guide_matrix %*% fwm))/
    (1+exp((guide_matrix %*% fwm)))
  return(guide_efficiency)
}