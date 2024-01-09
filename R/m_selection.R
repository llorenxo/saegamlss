#' Small area model selection within GAMLSS
#'
#' @description The function select the final SAE-GAMLSS model following three steps as described by Mori and Ferrante (2023)
#'
#' @param sample_data The dataset with sampled units
#' @param y The dependent variable
#' @param nRS Number of loop to do with the RS() algorithm
#' @param nCG Number of loop to do with the CG() algorithm
#' @param ndis Number of distribution to be consider at the first step. Default is 3
#' @param R Number of loop to be done within the k-fold cross validation
#' @param k Number of fold in k-fold cross validation. Default is 7
#' @param f_cov A formula containing all the possible covariates and/or additive terms (i.e. x1+x2+x3+random(x4))
#' @param type Type of step 2 to be done. A or B. Default is B (see references)
#' @param fix_dis A distribution to be tested even if is not selected within the ndis at step 1
#' @param seed The seed
#'
#' @return A list with the results of each steps.
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
#'   sigmae = 22, costh = NULL, seed = 124
#' )
#'
#' #Adding x2
#' set.seed(1234)
#' sample_data=as.data.frame(sample_data[[1]])
#' sample_data$x2=rnorm(40, 0, 1)
#'
#' m_selection(sample_data=sample_data, y = sample_data$y,
#'             f_cov= ~x1+x2+random(sa),nRS=20, nCG=20, ndis=2,
#'             R=2, k=2, type="A", fix_dis="NO", seed=123)

m_selection <- function(sample_data, y, f_cov,  nRS = NULL, nCG = NULL,
                        ndis=NULL, R=NULL, k=NULL, type=NULL, fix_dis=NULL, seed=NULL ){
  mixed <- NULL
  sel2 <- NULL
  if (is.null(ndis))  ndis <- 3
  if (is.null(R))  R <- 200
  if (is.null(k))  k <- 7
  if (is.null(seed)) seed <- 123
  options(warn=-1)
  set.seed(seed)
  sample_data$y = sample_data %>% dplyr::select(y) %>% dplyr::pull()

  #usethis

  is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))

  ## Recursively step down into list, removing all such objects
  rmNullObs <- function(x) {
    x <- Filter(Negate(is.NullOb), x)
    lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
  }

#step 1

  sel <- rownames(as.data.frame(gamlss::fitDist(y)$fits[1:ndis]))
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
                         method = mixed(substitute(nRS), substitute(nRG)),
                         family=substitute(sel[j]),  rand=rand_i)
      }, error = function(e) {
        cat("Error in iteration", j, "\n")
      })



      #try((gamlss::gamlssCV(formula=f1, sigma.fo=f2,
      #                            nu.fo=f3, tau.fo=f4,
      #                             data=as.data.frame(sample_data),
      #                             method = mixed(substitute(nRS), substitute(nRG)),
      #                              family=substitute(sel[j]),  rand=rand_i)  )
      #           , silent=T, outFile = getOption("try.outFile", default = stderr()))

      if (inherits(g3, "gamlssCV", which = FALSE)){
        a[[i, j]] <- c(CV(g3))
      }
    }
  }

  step3=as.data.frame(colMeans(a, na.rm=TRUE))
  colnames(step3)="Values"
  step3$Dist <- as.vector(sel)
  rm(sel2, envir = .GlobalEnv)
  names(select_v)=sel
  for (i in 1:length(names(select_v))) {
    names(select_v[[i]])=c("mu_f", "sigma_f", "nu_f", "tau_f")

  }
  select_v <- rmNullObs(select_v)

  return(list("step1"=sel, "step2"=select_v, "step3"=step3))

}

