#Not exp

theil.weibull <-function (sigma){
  theil <- sigma^{-1}*digamma(1+1/sigma)-log(gamma(1+1/sigma))
  return(theil)
}
