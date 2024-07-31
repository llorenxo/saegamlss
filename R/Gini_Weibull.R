#Not Exp


gini.weibull <- function(sigma){
  gini=1-2^{-1/sigma}
  return(gini)
}
