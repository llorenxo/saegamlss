#' Theil index Log-Normal distribution
#' @description Compute the Theil index with Log-Normal distribution
#'
#' @param sigma The estimated value of sigma
#'
#' @return The Theil index
#'

theil.logno <- function(sigma){
  theil=sigma^2/2
  return(theil)
}
