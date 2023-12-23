#' Atkinson index Log-Normal distribution
#' @description Compute the Atkinson index with Log-Normal distribution
#'
#' @param sigma The estimated value of sigma
#' @param epsilon The value for the poverty aversion parameter. Default value is set to 1.
#'
#' @return The Atkinson index
#'


atkinson.logno <- function(sigma, epsilon=NULL){
  if(is.null(epsilon))epsilon=1

  atkinson=1-exp(-0.5*sigma^2*epsilon)
  return(atkinson)

}
