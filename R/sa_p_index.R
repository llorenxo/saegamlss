#' Estimation of poverty indicators with simplified SAE-GAMLSS
#'
#' @description Estimate the values of poverty indicators in small areas using the simplified SAE-GAMLSS.
#'
#' @param fdis The assumed distribution. Options are: GB2 (Generalized Beta of 2-type), GA (Gamma), EXP (Exponential), LOGNO (Log-Normal), PA (Pareto), WE (Weibull)
#' @param index The index to be estimated ("Gini", "Theil" or "Atkinson"). Default is all
#' @param epsilon The value for the poverty aversion parameter. Default value is set to 1
#' @param sa The name of the variable to be used as "small area"
#' @param data The dataset with sampled data
#' @param y The dependent variable
#' @param sigma.f Logical value if TRUE (default) a random effect is used for sigma
#' @param nu.f Logical value if TRUE (default) a random effect is used for nu
#' @param tau.f Logical value if TRUE (default) a random effect is used for tau
#' @param w Sample weights
#' @param seed The seed. Default is 124
#'
#' @return An object of class "saegamlss_class" with the estimated indicators and the value of
#' the estimated distribution parameters for each area
#' @export
#'
#' @author Lorenzo Mori and Maria Rosaria Ferrante
#' @examples
#'
#'
#' data <- data.frame("y"= rLOGNO(1000, mu=10, sigma=0.8),
#'                    "sa" = as.factor(rep(c(1,2,3,4,5),200)))
#'
#' sa_p_index(data=data, y=data$y, sa = data$sa, fdis="LOGNO", index="Gini")

sa_p_index <- function (data, y, sigma.f=TRUE, nu.f=TRUE, tau.f=TRUE,  sa, w=NULL, fdis,
                        index=NULL, epsilon=NULL, seed = 123){

  mixed <- NULL
  f1 <- y ~1 + random(as.factor(sa))
  f2 <- f3 <- f4 <- y ~1


  if(isTRUE(sigma.f)) f2 <- f1
  if(isTRUE(nu.f)) f3 <- f1
  if(isTRUE(tau.f)) f4 <- f1

  if(is.null(epsilon)) epsilon=1

  set.seed(seed)

  gamlss_reg=gamlss::gamlss(f1, sigma.fo=f2, nu.fo=f3, tau.fo=f4, weights = w,
                              trace = F, family = substitute(fdis), data = as.data.frame(data),
                              method = mixed(10,10))




  data$mu_d <- predictAll(gamlss_reg, data=data, newdata=data)$mu
  data$sigma_d <- predictAll(gamlss_reg, data=data, newdata=data)$sigma
  data$nu_d <- predictAll(gamlss_reg, data=data, newdata=data)$nu
  data$tau_d <- predictAll(gamlss_reg, data=data, newdata=data)$tau



   sa <- data %>% dplyr::pull(sa)
   l_sa <- length(data %>% dplyr::distinct(sa) %>% dplyr::pull())
   gini_gamlss <- array()
   theil_gamlss <- array()
   atkinson_gamlss <- array()
   est_gamlss <- data.frame("sa"=data %>% dplyr::distinct(sa) %>% dplyr::pull(),
                            "mu_est"=rep(NA, l_sa),
                            "sigma_est"=rep(NA, l_sa),
                            "nu_est"=rep(NA, l_sa),
                            "tau_est"=rep(NA, l_sa))

    for (i in 1:l_sa) {
      a <- subset(data, sa==i)
      gini_gamlss <- rbind(gini_gamlss, p_index(mu=a$mu_d[1], sigma=a$sigma_d[1],
                                             nu=a$nu_d[1], tau=a$tau_d[1], fdis= fdis, index="Gini", epsilon = epsilon ))
      theil_gamlss <- rbind(theil_gamlss, p_index(mu=a$mu_d[1], sigma=a$sigma_d[1],
                                               nu=a$nu_d[1], tau=a$tau_d[1], fdis=fdis, index="Theil", epsilon = epsilon ))
      atkinson_gamlss <- rbind(atkinson_gamlss, p_index(mu=a$mu_d[1], sigma=a$sigma_d[1],
                                                     nu=a$nu_d[1], tau=a$tau_d[1], fdis=fdis, index="Atkinson", epsilon = epsilon ))

      if (!is.null(a$mu_d[1]))est_gamlss[i,2]=a$mu_d[1]
      if (!is.null(a$sigma_d[1]))est_gamlss[i,3]=a$sigma_d[1]
      if (!is.null(a$nu_d[1]))est_gamlss[i,4]=a$nu_d[1]
      if (!is.null(a$tau_d[1]))est_gamlss[i,5]=a$tau_d[1]

    }
     est_gamlss <- est_gamlss[, colSums(is.na(est_gamlss)) < nrow(est_gamlss)]

    if (is.null(index)){

      result <- list("Gini"=gini_gamlss[-1], "Theil"=theil_gamlss[-1], "Atkinson"=atkinson_gamlss[-1], "estimates_par"=est_gamlss, "model"=gamlss_reg)

      } else if (index=="Gini"){

         result <- list("Gini"=gini_gamlss[-1], "estimates_par"=est_gamlss, "model"=gamlss_reg)

         } else if (index=="Theil"){

            result <- list("Theil"=theil_gamlss[-1], "estimates_par"=est_gamlss, "model"=gamlss_reg)

            } else {

              result <- list("Atkinson"=atkinson_gamlss[-1], "estimates_par"=est_gamlss, "model"=gamlss_reg)

              }

     attr(result, "class") <- "saegamlss_class"
     return(result)

}








