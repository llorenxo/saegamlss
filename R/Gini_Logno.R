#Not exp


gini.logno <- function(sigma){
  gini <- pracma::erf(sigma/2)
  return(gini)
}


