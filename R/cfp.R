#' Carbon footprint
#'
#' @description Function to compute the household carbon footprint
#'
#'
#' @param e Vector of emission intensities by industries
#' @param Phi Matrix in which the diagonal blocks represent domestic transaction flows of intermediate goods and services across industries, and the off-diagonal blocks represent the inter-country flows of intermediates via exports and imports
#' @param Xi Matrix having on the principal diagonal the inverse of the total supply
#' @param consumption A dataset where each row corresponds to the vector of household consumption and each column is a specific consumption group
#' @param d Final use matrix
#' @param cc Vector with the expenditure on combustible of the household
#' @param DHE The total volume of direcùt household emissions
#' @param N Number of the units in the population
#' @param eimp Emission intensities derived from countries that play a central role in the import
#' @param Im Matrix of imports
#' @param Ep Matrix having on the principal diagonal the inverse of the total output
#' @param cf Conversion factor matrix
#' @param classic A logical value if TRUE the data are to be rearranged manually before to apply cfp(). Default is FALSE
#'
#' @return The dataset consumption (given in input) with 4 more columns: DIR (Direct emissions), IND (Indirect emissions), EEI (EMbodied emissions in import) and CFP (Carbon footprint)
#' @export
#'
#' @examples
#'
#' # Classic FALSE
#' # 2 groups of classification on consumption
#' # 3 groups of classification on industries sector
#'
#'
#' DHE=10
#' cc=c(2,1)
#' Phi=matrix(data=c(1,2,3,1,4,5,1,0.8, 0.5), nrow=3, ncol=3)
#' Xi=matrix(data=c(0.01,0.02,0.03,0.01,0.04,0.05,0.01,0.08, 0.05), nrow=3,
#'    ncol=3)
#' e=c(0.9,0.4,0.5)
#' consumption=data.frame("Fish"=c(0.1,0.3), "Vegetables"=c(0.2,0.5))
#' d=matrix(data=c(3,2,5,1,4,5,1,0.1, 0.2), nrow=3, ncol=3)
#' N=10
#' eimp=c(0.7,0.3,0.5)
#' Im=matrix(data=c(3,2,5,1,4,5,1,0.1, 0.2), nrow=3, ncol=3)
#' Ep=matrix(data=c(0.001,0.002,0.003,0.001,0.004,0.005,0.01,0.08, 0.05),
#'           nrow=3, ncol=3)
#' cf=matrix(c(1,0,0,1,0.5,0.5) ,nrow=3, ncol=2)
#' cfp(e=e, Phi=Phi, Xi=Xi, consumption=consumption, d=d, cc=cc, DHE=DHE, N=N, eimp=eimp,
#'     Im=Im, Ep=Ep, cf=cf, classic=FALSE)
#'
#' # Classic TRUE
#' # 2 groups of classification on consumption
#' # 3 groups of classification on industries sectors
#' # Consumption are rearranged in such a way that correspond to the three industries sectors
#' # cf matrix is not necessary
#' # Note: here the rearrangement is random. In real application it has to be based on logical reason.
#' rm(cf)
#' consumption=data.frame("Fish"=c(0.1,0.2), "Vegetables"=c(0.1,0.3), "Packaging"=c(0.1, 0.3))
#' cfp(e=e, Phi=Phi, Xi=Xi, consumption=consumption, d=d, cc=cc, DHE=DHE, N=N, eimp=eimp,
#'     Im=Im, Ep=Ep,  classic=TRUE)
#'

cfp <- function(e, Phi, Xi, consumption, d, cc, DHE, N, eimp, Im, Ep, cf, classic=FALSE){

  IA=pracma::inv(diag(ncol=length(e), nrow=length(e))-Phi%*%Ep)
  IAs=pracma::inv(diag(ncol=length(e), nrow=length(e))-Phi%*%Xi)

  D=list()
  IND=array()
  EEI=array()
  DIR=array()
  for (i in 1:nrow(consumption)) {
    if (isFALSE(classic)){
      D[[i]]=(t(cf%*%as.vector(as.numeric(consumption[i,])))%*%d)
    } else {
      D[[i]]=(t(as.vector(as.numeric(consumption[i,])))%*%d)
    }
    IND[[i]]=t(e)%*% IA %*% t(D[[i]])

    EEI[[i]]=t(eimp)%*% Im %*% Xi %*% IAs %*% t(D[[i]])
    DIR[[i]]=cc[i]/sum(cc)*DHE*nrow(consumption)/N

  }


  consumption$EEI <- EEI; consumption$DIR <- DIR
  consumption$IND <- IND; consumption$CFP <- IND+DIR+EEI


  return(consumption)
}

