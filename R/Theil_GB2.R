#' Theil index GB2 distribution
#' @description Compute the Theil index with GB2 distribution
#'
#' @param mu The estimated value of mu
#' @param sigma The estimated value of sigma
#' @param nu The estimated value of nu
#' @param tau The estimated value of tau
#'
#' @return The Theil index

theil.gb2 <- function(mu, sigma, nu, tau){
  theil1=((digamma(nu+1/sigma)-digamma(tau-1/sigma))/sigma)
  theil2=log((gamma(nu+1/sigma)*gamma(tau-1/sigma))/(gamma(nu)*gamma(tau)))
  theil=theil1-theil2
  return(theil)
}
