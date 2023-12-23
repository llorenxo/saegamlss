#' Atkinson index Pareto distribution
#' @description Compute the Atkinson index with Pareto distribution
#'
#' @param sigma The estimated value of sigma
#' @param epsilon The value for the poverty aversion parameter. Default value is set to 1.
#'
#' @return The Atkinson index
#'


atkinson.pareto <- function(sigma, epsilon=NULL){
  if(is.null(epsilon))epsilon=1

  atkinson=1-(((sigma-1)/sigma)*((sigma)/(sigma+epsilon-1))^(1/(1-epsilon)))
  return(atkinson)
}
