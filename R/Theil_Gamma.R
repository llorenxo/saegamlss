#' Theil index Gamma distribution
#' @description Compute the Theil index with Gamma distribution
#'
#' @param sigma The estimated value of sigma
#'
#' @return The Theil index


theil.gamma <-function (sigma){
  theil=digamma(1+sigma^(-2))-log(sigma^(-2))
  return(theil)
}
