#Not Exp

theil.gb2 <- function(mu, sigma, nu, tau){
  theil1 <- ((digamma(nu+1/sigma)-digamma(tau-1/sigma))/sigma)
  theil2 <- log((gamma(nu+1/sigma)*gamma(tau-1/sigma))/(gamma(nu)*gamma(tau)))
  theil <- theil1-theil2
  return(theil)
}
