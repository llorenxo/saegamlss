#' Gini index Gamma distribution
#' @description Compute the Gini index with Gamma distribution
#'
#' @param sigma The estimated value of sigma
#'
#' @return The Gini index
#'



gini.gamma <- function(sigma){
  gini <- gamma((2*sigma^(-2)+1)/2)/(sigma^(-2)*gamma(sigma^(-2)/2)*sqrt(pi))
  return(gini)
}
