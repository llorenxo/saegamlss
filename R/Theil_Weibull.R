#' Theil index Weibull distribution
#' @description Compute the Theil index with Weibull distribution
#'
#' @param sigma The estimated value of sigma
#'
#' @return The Theil index


theil.weibull <-function (sigma){
  theil=sigma^{-1}*digamma(1+1/sigma)-log(gamma(1+1/sigma))
  return(theil)
}
