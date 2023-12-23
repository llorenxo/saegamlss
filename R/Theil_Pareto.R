#' Theil index Pareto distribution
#' @description Compute the Theil index with Pareto distribution
#'
#' @param sigma The estimated value of sigma
#'
#' @return The Theil index


theil.pareto <-function (sigma){
  theil=log(1-1/sigma)+(sigma-1)^(-1)
  return(theil)
}
