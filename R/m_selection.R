#' Small area model selection within GAMLSS
#'
#' @description The function select the final SAE-GAMLSS model following three steps as described by Mori and Ferrante (2023)
#'
#' @param sample_data The dataset with sampled units
#' @param y The dependent variable
#' @param kp the penalty for the GAIC (Step 1) with default values kp=2 the standard AIC.
#' @param ndis Number of distribution to be consider at the first step. Default is 3
#' @param R Number of loop to be done within the k-fold cross validation
#' @param k Number of fold in k-fold cross validation. Default is 7
#' @param f_cov A formula containing all the possible covariates and/or additive terms (i.e. x1+x2+x3+random(x4))
#' @param type Type of step 2 to be done. A or B. Default is B (see references)
#' @param fix_dis A distribution to be tested even if is not selected within the ndis at step 1
#' @param supp Suppress the warning in the three step. Default is TRUE
#' @param seed The seed. Default is 123
#'
#' @return An object of class "saegamlss_class", named results, with the results of each steps and the GAIC of all the tested distributions at step 1.
#' @note The summary (Step2) do not reports results on the random effects as usual for GAMLSS. See Step2 select_v for the random effects
#' @export
#'
#' @references Mori, L., & Ferrante, M. R. (2023). Small area estimation under unit-level generalized additive models for location, scale and shape. arXiv e-prints, arXiv-2302.
#'  Stasinopoulus et al. (2019). Flexible regression and smoothing using GAMLSS in R.
#'  Rigby, R. A., & Stasinopoulos, D. M. (2005). Generalized additive models for location, scale and shape. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(3), 507–554.
#' @author Lorenzo Mori and Maria Rosaria Ferrante
#'
#' @examples
#'
#' sample_data=data_gen(
#'   Ni = rep(10, 4), D = 4, M = 1, ty = "no", k = 100, b1 = 4,
#'   x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
#'   x3 = NULL, b4 = NULL, x4 = NULL, xh = NULL,
#'   Dis = rNO, l = c(identity), sigma = 6, sigmah = NULL,
#'   sigmae = 22, costh = NULL
#' )
#'
#' #Adding x2
#' set.seed(1234)
#' sample_data=as.data.frame(sample_data[[1]])
#' sample_data$x2=rnorm(40, 0, 1)
#'
#' m_selection(sample_data=sample_data, y = sample_data$y,
#'             f_cov= ~x1+x2+random(sa), ndis=1,
#'             R=2, k=1, type="A", fix_dis="NO")

m_selection <- function(sample_data, y, f_cov, kp=2,
                        ndis=3, R=200, k=7, type=NULL, fix_dis=NULL,
                        supp=TRUE, seed = 123) {

  set.seed(seed)

  sel2 <- NULL

  if(isTRUE(supp)) options(warn=-1)

  sample_data$y = sample_data %>% dplyr::select(y) %>% dplyr::pull()

#step 1
  res_step1 <- gamlss::fitDist(y, k=kp)
  sel <- rownames(as.data.frame(res_step1$fits[1:ndis]))
  if (!is.null(fix_dis)) sel <- c(sel, fix_dis)
  sel <- unique(sel)

#step 2

  select_v <- list()
  step2 <- list()



  for (j in 1:length(sel)) {

    sel2 <<- sel[[j]]

    if ( type=="B"){

    mod0 <- gamlss::gamlss(y~1, data=as.data.frame(sample_data), family=sel2)

    a <- gamlss::stepGAICAll.B(mod0,
                               scope=list(lower=y~1,upper=f_cov))


    if(!is.null(a$mu.coefficients)){

      "mu_f"= stats::as.formula(paste("y~",
                               paste(names(a$mu.coefficients)
                                     [2:length(a$mu.coefficients)],
                                     collapse = " + ")))
    } else {
      "mu_f"=NULL
    }
    if (!is.null(a$sigma.coefficients)){
      "sigma_f"= stats::as.formula(paste("y~",
                                  paste(names(a$sigma.coefficients)
                                        [2:length(a$sigma.coefficients)],
                                        collapse = " + ")))
    } else {
      "sigma_f"=NULL
    }
    if (!is.null(a$nu.coefficients)){
      "nu_f"=stats::as.formula(paste("y~",
                              paste(names(a$nu.coefficients)
                                    [2:length(a$nu.coefficients)],
                                    collapse = " + ")))
    } else {
      "nu_f"=NULL
    }
    if (!is.null(a$tau.coefficients)){
      "tau_f"=stats::as.formula(paste("y~",
                               paste(names(a$tau.coefficients)
                                     [2:length(a$tau.coefficients)],
                                     collapse = " + ")))
    } else {
      "tau_f"=NULL
    }



    select_v[[j]] <-  list(mu_f, sigma_f, nu_f, tau_f)
    } else {
      mod0 <- gamlss::gamlss(y ~1 , data=sample_data , family=sel2)

      a <- gamlss::stepGAICAll.B(mod0,
                                 scope=list(lower=y~1,upper=f_cov))



      if(!is.null(a$mu.coefficients)){

        "mu_f"= stats::as.formula(paste("y~",
                               paste(names(a$mu.coefficients)
                                     [2:length(a$mu.coefficients)],
                                     collapse = " + ")))
      } else {
        "mu_f"=NULL
      }
      if (!is.null(a$sigma.coefficients)){
        "sigma_f"= stats::as.formula(paste("y~",
                                  paste(names(a$sigma.coefficients)
                                        [2:length(a$sigma.coefficients)],
                                        collapse = " + ")))
      } else {
        "sigma_f"=NULL
      }
      if (!is.null(a$nu.coefficients)){
        "nu_f"=stats::as.formula(paste("y~",
                              paste(names(a$nu.coefficients)
                                    [2:length(a$nu.coefficients)],
                                    collapse = " + ")))
      } else {
        "nu_f"=NULL
      }
      if (!is.null(a$tau.coefficients)){
        "tau_f"=stats::as.formula(paste("y~",
                               paste(names(a$tau.coefficients)
                                     [2:length(a$tau.coefficients)],
                                     collapse = " + ")))
      } else {
        "tau_f"=NULL
      }


      select_v[[j]] <-  list(mu_f, sigma_f, nu_f, tau_f)

    }
  }

#step 3

  a=matrix(ncol=length(sel), nrow=R)

  for (i in 1:R){
        rand_i <- sample (k , nrow(sample_data), replace=TRUE)
      for (j  in 1:length(sel)) {
        f1=stats::as.formula(select_v[[j]][[1]])
        f2=stats::as.formula(select_v[[j]][[2]])
        f3=stats::as.formula(select_v[[j]][[3]])
        f4=stats::as.formula(select_v[[j]][[4]])
      g3 <-tryCatch({
        gamlss::gamlssCV(formula=f1, sigma.fo=f2,
                         nu.fo=f3, tau.fo=f4,
                         data=as.data.frame(sample_data),
                         method = mixed(100, 100),
                         family=substitute(sel[j]),  rand=rand_i)
      }, error = function(e) {
        cat("Error in iteration", j, "\n")
        })

      if (inherits(g3, "gamlssCV", which = FALSE)){
        a[[i, j]] <- c(CV(g3))
      }
    }
  }

  step3=as.data.frame(colMeans(a, na.rm=TRUE))
  colnames(step3) = "Values"
  step3$Dist <- as.vector(sel)
  rm(sel2, envir = .GlobalEnv)
  names(select_v)=sel
  for (i in 1:length(names(select_v))) {
    names(select_v[[i]])=c("mu_f", "sigma_f", "nu_f", "tau_f")

  }
  select_v <- rmNullObs(select_v)
  results=list("step1"=sel, "step2"=select_v, "step3"=step3, "GAIC values"=res_step1,
               "sample_data"=sample_data, y="y")
  attr(results, "class") <- "saegamlss_class"
  return(results)
}

