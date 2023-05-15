#' Create the non-sample dataset
#' @description Subset the non-sample part of the population
#' @param data population dataset
#' @param sample sample dataset
#' @param id name of id, must be the same in both the datasets
#'
#' @return A dataset of non-sample units
#' @export
#'
#' @examples #sample data with a sample fraction of 0.1
#' @examples library(splitstackshape)
#' @examples data=data.frame("y"=runif(100, 0, 100), "id"=c(1:100), "sa"=rep(1:5, 20))
#' @examples #sample data
#' @examples sample=stratified(data, "sa", size=0.1)
#' @examples #nonsample data
#' @examples nonsample=nonsample(data=data, sample=sample, id=id)
#'
#' @author Lorenzo Mori and Maria Rosaria Ferrante
#' @references Mori, L., & Ferrante, M. R. (2023). Small area estimation under unit-level generalized additive models for location, scale and shape. arXiv e-prints, arXiv-2302.
#'
nonsample= function(data, sample, id){
  '%!in%'<- function(x,y)!('%in%'(x,y))
  nonsample=subset(data, data$id%!in%sample$id)
  return("nonsample"=as.data.frame(nonsample))
}
