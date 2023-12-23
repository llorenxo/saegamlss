#' Gini index Weibull distribution
#' @description Compute the Gini index with Weibull distribution
#'
#' @param sigma The estimated value of sigma
#'
#' @return The Gini index
#'



gini.weibull <- function(sigma){
  gini=1-2^{-1/sigma}
  return(gini)
}
