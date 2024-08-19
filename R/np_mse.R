#' Non-parametric MSE for simplified SAE-GAMLSS
#' @description Compute the Non-parametric MSE for simplified SAE-GAMLSS for three possible different
#' indicators (Gini, Theil, Atkinson)
#'
#' @param data A dataset containing sampled household values
#' @param y The dependent variable name
#' @param sa The Small Area domains name
#' @param ncomp The number of components of each household
#' @param R The number of loops to be performed. Default is 200
#' @param sigma.f Logical value if TRUE (default) a random effect is used for sigma
#' @param nu.f Logical value if TRUE (default) a random effect is used for nu
#' @param tau.f Logical value if TRUE (default) a random effect is used for tau
#' @param w Sample weights
#' @param index One index to be estimated ("Gini", "Theil" or "Atkinson"). Default is "all"
#' @param epsilon The value for the poverty aversion parameter. Default value is set to 1
#' @param fdis The assumed distribution. Options are: GB2 (Generalized Beta of 2-type), GA (Gamma), EX (Exponential), LOGNO (Log-Normal), PA (Pareto), WE (Weibull)
#' @param seed The seed. Default is 123
#'
#'
#'
#' @return An object of class "saegamlss_class" containing the values of the MSE for each area and for each index
#'
#' @author Lorenzo Mori and Maria Rosaria Ferrante
#' @export
#'
#' @examples
#'
#' ##################
#' ###Using s_data###
#' ##################
#'
#' set.seed(124)
#'
#'
#' np <- np_mse(data = s_data, y = "y", sa = "sa",
#'                        ncomp = "ncomp", fdis="LOGNO",
#'                        index="Gini", seed = 124,
#'                        R=2)
#'
#' np$Gini.MSE
#'
#'
#'

np_mse <- function(data, y, sa, ncomp,  R = 200, sigma.f = TRUE, nu.f = TRUE,
                               tau.f = TRUE, w = NULL, fdis,
                               index = "all", epsilon = 1, seed = 123){
  set.seed(seed)
  y <- data[[y]]
  sa <- data[[sa]]
  ncomp <- data[[ncomp]]
  f1 <- y ~1 + random(as.factor(sa))
  f2 <- f3 <- f4 <- y ~1

  if(isTRUE(sigma.f))f2 <- f1
  if(isTRUE(nu.f))f3 <- f1
  if(isTRUE(tau.f))f4 <- f1

  if (is.null(w)) data <- data %>% mutate(w = rep(1:nrow(data)))

  s2 <- data %>% dplyr::select("y", "sa", "ncomp", "w")

  s2$id <- c(1:nrow(s2))

  l_sa <- length(data %>% dplyr::distinct(sa) %>% dplyr::pull())


  dif3 <- rep(0, l_sa)
  dif4 <- rep(0, l_sa)
  dif5 <- rep(0, l_sa)


for (t2 in 1:R){

  message("Processing mse loop ", t2, " of ", R)

  s_b <- s2[1,]
  ss <- s2

  for (i in 1:l_sa){

    s_b_sa1 <- subset(ss, sa==i)
    s_b_sa1$id <- c(1: nrow(s_b_sa1))

    if(nrow(s_b_sa1)>1 ){

      s_b_sa <- s_b_sa1[sample(1:nrow(s_b_sa1), nrow(s_b_sa1)-1, replace=TRUE),]

      } else {

      s_b_sa <- s_b_sa1[sample(1:nrow(s_b_sa1), nrow(s_b_sa1), replace=TRUE),]

      }

    for (j in 1:nrow(s_b_sa)){

      s_b_sa [[j,5]] <- s_b_sa[[j,5]] * nrow(s_b_sa1) / nrow(s_b_sa) * rep(as.vector(table(s_b_sa$id)),
                                                                 as.vector(table(s_b_sa$id)))[j]
    }

    s_b_sa$w < -s_b_sa$w * sum(s_b_sa1$w) / sum(s_b_sa$w)

    s_b<-rbind(s_b, s_b_sa)

  }

  s_b <- s_b[-1,]
  s1 <- s_b
  s1 <- as.data.frame(lapply(s1, rep, s1$ncomp))

  cont<-gamlss.control(c.crit = 0.01, n.cyc = 15, mu.step = 1, sigma.step = 1, nu.step = 1,
                      tau.step = 1, gd.tol = Inf, iter = 0, trace = F, autostep = TRUE)


    ga<-gamlss::gamlss(f1, sigma.fo=f2, nu.fo=f3, tau.fo=f4, weights = w,
                      trace = FALSE, family = substitute(fdis), data = as.data.frame(s1),
                      method = RS(100))




  predi <- predictAll(ga, data=s1, newdata=s1)

  if (!is.null(predi$mu)) s1$mu_d <- predi$mu
  if (!is.null(predi$sigma)) s1$sigma_d <- predi$sigma
  if (!is.null(predi$nu)) s1$nu_d <- predi$nu
  if (!is.null(predi$tau)) s1$tau_d <- predi$tau


  if (length(s1$mu_d)<1)   s1$mu_d=s1$mu_d[1:nrow(s1)]

  if (length(s1$nu_d)<1)   s1$nu_d=s1$nu_d[1:nrow(s1)]

  if (length(s1$tau_d)<1)  s1$tau_d=s1$tau_d[1:nrow(s1)]

  if (length(s1$sigma_d)<1) s1$sigma_d=s1$sigma_d[1:nrow(s1)]



  atkinson_gamlss_boot <- 0
  gini_gamlss_boot <- 0
  theil_gamlss_boot <- 0


  for (i in 1:l_sa) {

    a <- subset(data, sa==i)


    if (index == "all" | index == "Gini") {

      gini_gamlss_boot <- rbind(gini_gamlss_boot, p_index(mu=s1$mu_d[1], sigma=s1$sigma_d[1],
                                           nu=s1$nu_d[1], tau=s1$tau_d[1], fdis= fdis, index="Gini", epsilon = epsilon )$index$P_Gini)

     }

    if (index == "all" | index == "Theil"){

       theil_gamlss_boot <- rbind(theil_gamlss_boot, p_index(mu=s1$mu_d[1], sigma=s1$sigma_d[1],
                                             nu=s1$nu_d[1], tau=s1$tau_d[1], fdis=fdis, index="Theil", epsilon = epsilon )$index$P_Theil)

    }

    if (index == "all" | index == "Atkinson"){

    atkinson_gamlss_boot <- rbind(atkinson_gamlss_boot, p_index(mu=s1$mu_d[1], sigma=s1$sigma_d[1],
                                                   nu=s1$nu_d[1], tau=s1$tau_d[1], fdis=fdis, index="Atkinson", epsilon = epsilon )$index$P_Atkinson)
     }

  }

  dif3 <- rbind(dif3, atkinson_gamlss_boot[-1])
  dif4 <- rbind(dif4, gini_gamlss_boot[-1])
  dif5 <- rbind(dif5, theil_gamlss_boot[-1])

  }

  dif3 <- as.data.frame(lapply(as.data.frame(dif3), as.numeric))
  dif4 <- as.data.frame(lapply(as.data.frame(dif4), as.numeric))
  dif5 <- as.data.frame(lapply(as.data.frame(dif5), as.numeric))


  if (index == "all"){

  result <- list("Gini.MSE"=colMeans(dif4[-1,]), "Theil.MSE"=colMeans(dif5[-1,]), "Atkinson.MSE"=colMeans(dif3[-1,]))

  } else if (index=="Gini"){

  result <- list("Gini.MSE"=colMeans(dif4[-1,]))

  } else if (index=="Theil"){

  result <- list("Theil.MSE"=colMeans(dif5[-1,]))

  } else {

  result <- list("Atkinson.MSE"=colMeans(dif3[-1,]))

  }


  result <- lapply(result, function(x) {
    names(x) <- levels(sa)
    return(x)
  })

  attr(result, "class") <- "saegamlss_class"
  return(result)

}
