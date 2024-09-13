#' Estimation of poverty indicators with simplified SAE-GAMLSS
#'
#' @description Estimate the values of poverty indicators in small areas using the simplified SAE-GAMLSS.
#'
#' @param sample A dataset with sampled units
#' @param fdis The assumed distribution. Options are: GB2 (Generalized Beta of 2-type), GAMMA (Gamma), EXP (Exponential), LOGNO (Log-Normal), PARETO (Pareto), WEI (Weibull)
#' @param index The index to be estimated ("Gini", "Theil" or "Atkinson"). Default is all
#' @param epsilon The value for the poverty aversion parameter. Default value is set to 1
#' @param sa The name of the variable to be used as "small area"
#' @param y The dependent variable
#' @param sigma.f Logical value if TRUE (default) a random effect is used for sigma
#' @param nu.f Logical value if TRUE (default) a random effect is used for nu
#' @param tau.f Logical value if TRUE (default) a random effect is used for tau
#' @param w Sample weights. Default is null
#' @param seed The seed. Default is 123
#'
#' @return An object of class "saegamlss" with the estimated indicators and the value of
#' the estimated distribution parameters for each area
#' @export
#'
#' @examples
#'
#' ##################
#' ###Using s_data###
#' ##################
#'
#'
#' index_est <- sa_p_index(sample = s_data, y = "y",
#'                         sa = "sa", fdis = "LOGNO",
#'                         sigma.f = TRUE, index = "all")
#'
#' index_est$estimates_par
#' index_est
#'
#' @author Lorenzo Mori and Maria Rosaria Ferrante

sa_p_index <- function (sample, y, sigma.f = TRUE, nu.f = TRUE, tau.f = TRUE,  sa, w = NULL, fdis,
                        index = "all", epsilon = 1, seed = 123){


  input_var <- list("sample" = sample,
                    "y" = y,
                    "sa" = sa,
                    "sigma.f" = sigma.f,
                    "nu.f" = nu.f,
                    "tau.f" = tau.f,
                    "sa" = sa,
                    "w" = w,
                    "fdis" = fdis,
                    "index" =  index,
                    "epsilon" = epsilon,
                    "seed" = seed
                     )

  sa <- sample[[sa]] %>% as.factor()
  sa_names <- unique(sa)
  y <- sample[[y]]
  f1 <- y ~1 + random(sa)
  f2 <- f3 <- f4 <- y ~1

  if(isTRUE(sigma.f)) f2 <- f1
  if(isTRUE(nu.f)) f3 <- f1
  if(isTRUE(tau.f)) f4 <- f1

  if (is.null(w))  {

    sample <- sample %>% mutate(w = rep(1:nrow(sample)))

  } else {

    sample <- sample %>% mutate(w = (sample %>% dplyr::select(w) %>% dplyr::pull()))
  }


  set.seed(seed)

  gamlss_reg=gamlss::gamlss(f1, sigma.fo=f2, nu.fo=f3, tau.fo=f4, weights = w,
                              trace = F, family = substitute(fdis), data = as.data.frame(sample),
                              method = mixed(100,100))


  sample$mu_d <- predictAll(gamlss_reg, data=sample, newdata=sample)$mu
  sample$sigma_d <- predictAll(gamlss_reg, data=sample, newdata=sample)$sigma
  sample$nu_d <- predictAll(gamlss_reg, data=sample, newdata=sample)$nu
  sample$tau_d <- predictAll(gamlss_reg, data=sample, newdata=sample)$tau



   l_sa <- length(sample %>% dplyr::distinct(sa) %>% dplyr::pull())
   gini_gamlss <- array()
   theil_gamlss <- array()
   atkinson_gamlss <- array()
   est_gamlss <- data.frame("sa"=sample %>% dplyr::distinct(sa) %>% dplyr::pull(),
                            "mu_est"=rep(NA_integer_, l_sa),
                            "sigma_est"=rep(NA_integer_, l_sa),
                            "nu_est"=rep(NA_integer_, l_sa),
                            "tau_est"=rep(NA_integer_, l_sa))
    con = 0

    for (i in sa_names) {

      con = con + 1

      a <- subset(sample, sa==i)

      if (index == "all" | index =="Gini"){


      gini_gamlss <- rbind(gini_gamlss, p_index(mu=a$mu_d[1], sigma=a$sigma_d[1],
                                                nu=a$nu_d[1], tau=a$tau_d[1], fdis= fdis,
                                                index="Gini", epsilon = epsilon )$index$P_Gini)

      }

      if (index == "all" | index =="Theil"){


      theil_gamlss <- rbind(theil_gamlss, p_index(mu=a$mu_d[1], sigma=a$sigma_d[1],
                                               nu=a$nu_d[1], tau=a$tau_d[1], fdis=fdis,
                                               index="Theil", epsilon = epsilon )$index$P_Theil)
      }

      if (index == "all" | index =="Atkinson"){


      atkinson_gamlss <- rbind(atkinson_gamlss, p_index(mu=a$mu_d[1], sigma=a$sigma_d[1],
                                                     nu=a$nu_d[1], tau=a$tau_d[1], fdis=fdis,
                                                     index="Atkinson", epsilon = epsilon )$index$P_Atkinson)

      }

      if (!is.null(a$mu_d[1])) est_gamlss[con,2]=a$mu_d[1]
      if (!is.null(a$sigma_d[1])) est_gamlss[con,3]=a$sigma_d[1]
      if (!is.null(a$nu_d[1])) est_gamlss[con,4]=a$nu_d[1]
      if (!is.null(a$tau_d[1])) est_gamlss[con,5]=a$tau_d[1]

    }

    est_gamlss <- est_gamlss[, colSums(is.na(est_gamlss)) < nrow(est_gamlss)]
    rownames(est_gamlss) <- levels(sa)


    if (index=="all"){

      result <- list("Gini"=gini_gamlss[-1] %>% setNames(levels(sa)) ,
                     "Theil"=theil_gamlss[-1] %>% setNames(levels(sa)),
                     "Atkinson"=atkinson_gamlss[-1] %>% setNames(levels(sa)),
                     "estimates_par"=est_gamlss)

      } else if (index=="Gini"){

         result <- list("Gini"=gini_gamlss[-1]  %>% setNames(levels(sa)),
                        "estimates_par"=est_gamlss)

         } else if (index=="Theil"){

            result <- list("Theil"=theil_gamlss[-1] %>% setNames(levels(sa)),
                           "estimates_par"=est_gamlss)

            } else {

              result <- list("Atkinson"=atkinson_gamlss[-1] %>% setNames(levels(sa)),
                             "estimates_par"=est_gamlss)

            }

    input_var$fit = gamlss_reg
    input_var$sa_names = sa_names
    result$input_var <- input_var

    attr(result, "class") <- "saegamlss"
     return(result)
}









