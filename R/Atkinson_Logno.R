#Not Exp


atkinson.logno <- function(sigma, epsilon=NULL){
  if(is.null(epsilon)) epsilon <- 1

  atkinson <- 1-exp(-0.5*sigma^2*epsilon)
  return(atkinson)

}
