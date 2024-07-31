#Not exp

theil.pareto <-function (sigma){
  theil <- log(1-1/sigma)+(sigma-1)^(-1)
  return(theil)
}
