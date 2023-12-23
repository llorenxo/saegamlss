#' Theil index exponential distribution
#' @description Compute the Theil index with exponential distribution
#'
#' @param mu The estimated value of mu
#'
#' @return The Theil index
#'

theil.exp <- function (mu){
   theil=1-mu
   return(theil)
  }