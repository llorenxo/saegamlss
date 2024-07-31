#' Gini index Log-Normal distribution
#' @description Compute the Gini index with Log-Normal distribution
#'
#' @param sigma The estimated value of sigma
#'
#' @return The Gini index
#'


gini.logno <- function(sigma){
  gini <- pracma::erf(sigma/2)
  return(gini)
}


