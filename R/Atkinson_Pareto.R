#Not Exp


atkinson.pareto <- function(sigma, epsilon=NULL){
  if(is.null(epsilon))epsilon <- 1

  atkinson <- 1-(((sigma-1)/sigma)*((sigma)/(sigma+epsilon-1))^(1/(1-epsilon)))
  return(atkinson)
}
