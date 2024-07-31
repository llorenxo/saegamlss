#Not exp


theil.gamma <-function (sigma){
  theil <- digamma(1+sigma^(-2))-log(sigma^(-2))
  return(theil)
}
