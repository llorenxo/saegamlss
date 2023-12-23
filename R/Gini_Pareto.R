#' Gini index Pareto distribution
#' @description Compute the Gini index with Pareto distribution
#'
#' @param sigma The estimated value of sigma
#'
#' @return The Gini index
#'



gini.pareto <- function(sigma){
  gini=(2*sigma-1)^(-1)
  return(gini)
}
