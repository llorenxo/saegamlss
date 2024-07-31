#Not Exp


gini.gamma <- function(sigma){
  gini <- gamma((2*sigma^(-2)+1)/2)/(sigma^(-2)*gamma(sigma^(-2)/2)*sqrt(pi))
  return(gini)
}
