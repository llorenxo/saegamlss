#' Distribution based Gini, Theil and Atkinson index estimation
#'
#' @description Provide the estimated values for Gini, Theil, and Atkinson indices using the closed-form expression for an assumed distribution.
#' @param mu The estimated value of mu
#' @param sigma The estimated value of sigma
#' @param nu The estimated value of nu
#' @param tau The estimated value of tau
#' @param fdis The assumed distribution. Options are: GB2 (Generalized Beta of 2-type), GAMMA (Gamma), EXP (Exponential), LOGNO (Log-Normal), PARETO (Pareto), WEI (Weibull)
#' @param index The index to be estimated ("Gini", "Theil" or "Atkinson"). Default is all
#' @param epsilon The value for the poverty aversion parameter. Default is  1
#'
#' @return The estimated index (indicators)
#' @export
#' @examples
#'
#'
#' #Using GB2 distribution
#'
#' p_index(mu=5, sigma=2, nu=1, tau=1, fdis="GB2", index="Gini", epsilon=2)
#' p_index(mu=5, sigma=2, nu=1, tau=1, fdis="GB2", index="Theil", epsilon=2)
#' p_index(mu=5, sigma=2, nu=1, tau=1, fdis="GB2", index="Atkinson", epsilon=2)
#' p_index(mu=5, sigma=2, nu=1, tau=1, fdis="GB2", index="all", epsilon=2)
#'
#' #Using Log-Normal distribution
#'
#' p_index(sigma=0.8, fdis="LOGNO", index="Gini", epsilon=2)
#' p_index(sigma=0.8, fdis="LOGNO", index="Theil", epsilon=2)
#' p_index(sigma=0.8, fdis="LOGNO", index="Atkinson", epsilon=2)
#' p_index(sigma=0.8, fdis="LOGNO", index="all", epsilon=2)
#'
#' @author Lorenzo Mori and Maria Rosaria Ferrante

p_index  <- function (mu, sigma, nu, tau, fdis, index="all", epsilon=1) {


  if (index=="all"){

    if (fdis=="LOGNO"){

      gini <- gini.logno(sigma=sigma)
      theil <- theil.logno(sigma=sigma)
      atkinson <- atkinson.logno(sigma=sigma, epsilon=epsilon)

    } else if (fdis=="GB2"){

      gini <- GB2::gini.gb2(shape1=sigma, shape2=nu, shape3=tau)
      theil <- theil.gb2(mu=mu, sigma=sigma, nu=nu, tau=tau)
      atkinson  <- acid::atkinson.GB2(p=nu, a=sigma, b=mu, q=tau,epsilon = epsilon )

    }  else if (fdis=="EXP"){

      gini <- 0.5
      theil <- theil.exp(mu=mu)
      atkinson <- print("The Atkinson index for this distribution is not already implemented")

    } else if (fdis=="WEI"){

      gini <- gini.weibull(sigma=sigma)
      theil <- theil.weibull(sigma=sigma)
      atkinson <- print("The Atkinson index for this distribution is not already implemented")

    } else if (fdis=="GAMMA"){

      gini <- gini.gamma(sigma=sigma)
      theil <- theil.gamma(sigma=sigma)
      atkinson <- print("The Atkinson index for this distribution is not already implemented")

    } else if (fdis=="Pareto"){

      gini <- gini.pareto(sigma=sigma)
      theil <- theil.pareto(sigma=sigma)
      atkinson <- atkinson.pareto(sigma=sigma, epsilon=epsilon)

    }

    index <- list(P_Gini=gini, P_Theil=theil, P_Atkinson=atkinson)

  } else if (index=="Gini"){

    if (fdis=="LOGNO"){

      gini <- gini.logno(sigma=sigma)

    } else if (fdis=="GB2"){

      gini <- GB2::gini.gb2(shape1=sigma, shape2=nu, shape3=tau)

    }else if (fdis=="EXP"){

      gini <- 0.5

    } else if (fdis=="WEI"){

      gini <- gini.weibull(sigma=sigma)

    } else if (fdis=="GAMMA"){

      gini <- gini.gamma(sigma=sigma)

    } else if (fdis=="Pareto"){

      gini <- gini.pareto(sigma=sigma)

    }
    index=list(P_Gini=gini)

  } else if (index=="Theil"){

    if (fdis=="LOGNO"){

      theil <- theil.logno(sigma=sigma)

    } else if (fdis=="GB2"){

      theil <- theil.gb2(mu=mu, sigma=sigma, nu=nu, tau=tau)

    } else if (fdis=="EXP"){

      theil <- theil.exp(mu=mu)

    } else if (fdis=="WEI"){

      theil <- theil.weibull(sigma=sigma)

    } else if (fdis=="GAMMA"){

      theil <- theil.gamma(sigma=sigma)

    } else if (fdis=="Pareto"){

      theil <- theil.pareto(sigma=sigma)

    }

    index <- list(P_Theil=theil)

  } else {
    if (fdis=="LOGNO"){

      atkinson <- atkinson.logno(sigma=sigma, epsilon=epsilon)

    } else if (fdis=="GB2"){

      atkinson <- acid::atkinson.GB2(p=nu, a=sigma, b=mu, q=tau,epsilon = epsilon )

    } else if (fdis=="EXP"){

      atkinson <- print("The Atkinson index for this distribution is not already implemented")

    } else if (fdis=="WEI"){

      atkinson <- print("The Atkinson index for this distribution is not already implemented")

    } else if (fdis=="GAMMA"){

      atkinson <- print("The Atkinson index for this distribution is not already implemented")

    } else if (fdis=="Pareto"){

      atkinson <- atkinson.pareto(sigma=sigma, epsilon=epsilon)
    }

    index <- list(P_Atkinson=atkinson)
  }

  if ("P_Atkinson" %in% names(index)) {

  index = list(index = index, dist = fdis, Atkinson_parameter = epsilon)

  } else {

    index = list(index = index, dist = fdis)

  }

  attr(index, "class") <- "saegamlss"

  return(index)

}
