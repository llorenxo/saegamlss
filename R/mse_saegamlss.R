#' Bootstrap Mean Square Error for SAE GAMLSS
#'
#' @description Compute the bootstrap MSE for the estimation of the mean and the HCR
#' @param est An object obtained with est_saegamlss()
#' @param loop Number of loop of bootstrap. Default is 200
#' @param l A vector, of maximum length equal to 4, in which each space is a function (the inverse of the link-function) to be applied to the corresponding parameter, i.e. mu, sigma, nu and tau
#' @param Iden  TRUE if the Identity link function is used for mu. Default is FALSE
#' @param cov1 A matrix or a data frame with covariates for the whole population used for mu. If an intercept is used the first columns have to be a vector of 1
#' @param cov2 A matrix or a data frame with covariates for the whole population used for sigma. If an intercept is used the first columns have to be a vector of 1
#' @param cov3 A matrix or a data frame with covariates for the whole population used for nu. If an intercept is used the first columns have to be a vector of 1
#' @param cov4 A matrix or a data frame with covariates for the whole population used for tau. If an intercept is used the first columns have to be a vector of 1
#' @param seed The seed. Default is 123
#'
#' @return an object of class "saegamlss_class" containing the MSE of the mean or/and HCR for each area
#'  and all the values returned by est_saegamlss()
#' @export
#' @import gamlss
#' @import dplyr
#' @import splitstackshape
#' @import utils
#'
#' @examples
#'
#' ##################
#' ###Using s_data###
#' ##################
#'
#'
#' est <- est_saegamlss(
#'   sample = s_data, nonsample = pop_data, y_dip="y",
#'   sa="sa", f1 = y ~ x1 + random(sa), f2 = y ~ x2 + random(sa),
#'   f3 = NULL, f4 = NULL, fdis = NO, R = 20,
#'   Dis = rNO,  param = "mean",
#'   tau.fix = NULL, nu.fix = NULL
#' )
#'
#' data <- rbind(est$input_var$origindata  %>%
#'               dplyr::select( dplyr::all_of(colnames(est$input_var$nonsample))),
#'               est$input_var$nonsample
#'                )
#'
#' x <- data.frame(rep(1, nrow(data)), "x1" = data$x1,
#'                 "x2" = data$x2)
#'
#'
#' MSE <- mse_saegamlss(
#'   est = est, loop = 2,
#'   l = c(identity, exp), Iden = TRUE,
#'   cov1 = x[,1:2], cov2 = x[,c(1,3)], cov3 = NULL,
#'   cov4 = NULL
#' )
#'
#' MSE$
#'
#'
#' # MSE for a self-defined parameter
#'
#' est <- est_saegamlss(
#'   sample = s_data, nonsample = pop_data, y_dip = "y",
#'   sa = "sa", f1 = y ~ x1 + random(sa), f2 = y ~ x2 + random(sa),
#'   f3 = NULL, f4 = NULL, fdis = NO, R = 20,
#'   Dis = rNO, param = function(x) (mean(x^2)),
#'   tau.fix = NULL, nu.fix = NULL
#'  )
#'
#'
#' MSE <- mse_saegamlss(
#'   est = est, loop = 2,
#'   l = c(identity, exp), Iden = TRUE,
#'   cov1 = x[,1:2], cov2 = x[,c(1,3)], cov3 = NULL,
#'   cov4 = NULL
#' )
#'
#'
#' MSE$est_mse$MSE_param
#'
#' @references Mori, L., & Ferrante, M. R. (2023). Small area estimation under unit-level generalized additive models for location, scale and shape. arXiv e-prints, arXiv-2302.
#'  Graf, M., Marin, J. M., & Molina, I. (2019). A generalized mixed model for skewed distributions applied to small area estimation. Test, 28(2), 565–597.
#' @author Lorenzo Mori and Maria Rosaria Ferrante

mse_saegamlss <- function(est, loop = 200, l, Iden = FALSE,
                          cov1, cov2 = NULL, cov3 = NULL, cov4 = NULL, seed = 123) {

  set.seed(seed)


  data <- rbind( est$input_var$origindata  %>%
                 dplyr::select( dplyr::all_of(colnames(est$input_var$nonsample))),
                 est$input_var$nonsample
                 )
  data <- data %>% dplyr::mutate(id = row_number())

  sa <- est$input_var$sa
  samplesize <- est$input_var$origindata %>% dplyr::count(sa )%>% dplyr::pull()/ sum(est$input_var$Ni)
  D <- est$input_var$D
  Dis <- est$input_var$Dis
  Ni <- est$input_var$Ni
  Ni2 <- est$input_var$Ni + est$input_var$ni
  np <- count_arguments(Dis)-1

  param <- est$input_var$param

  mse_p <- mse <- MSE <- rep(0, D)
  u <- 0


  gam <- est$input_var$fit


  for (i in 1:loop) {

    sigm <- gamlss::getSmo(gam, what = "mu")$sigb

    if (np > 1) {

      sigs <- gamlss::getSmo(gam, what = "sigma")$sigb

    } else {

      sigs <- 0
    }

    if (np > 2) {

      sign <- gamlss::getSmo(gam, what = "nu")$sigb
    } else {

      sign <- 0
    }

    if (np > 3) {

      sigt <- gamlss::getSmo(gam, what = "tau")$sigb

    } else {

      sigt <- 0
    }

    if (is.null(sigs)) sigs <- 0

    if (is.null(sign)) sign <- 0

    if (is.null(sigt)) sigt <- 0






    u1 <- stats::rnorm(D, 0, sigm) # re x
    u2 <- stats::rnorm(D, 0, sigs) # re x
    u3 <- stats::rnorm(D, 0, sign) # re x
    u4 <- stats::rnorm(D, 0, sigt) # re x

    ################### Error

    sm <- rep(u1, Ni2)
    ss <- rep(u2, Ni2)
    sn <- rep(u3, Ni2)
    st <- rep(u4, Ni2)



    mm <- cov1
    ms <- cov2
    mn <- cov3
    mt <- cov4

    llm <- length(gam$mu.coefficients)
    lls <- length(gam$sigma.coefficients)
    lln <- length(gam$nu.coefficients)
    llt <- length(gam$tau.coefficients)

    if (length(l) == 1) {
      l1 <- l[[1]]
    } else if (length(l) == 2) {
      l1 <- l[[1]]
      l2 <- l[[2]]
    } else if (length(l) == 3) {
      l1 <- l[[1]]
      l2 <- l[[2]]
      l3 <- l[[3]]
    } else {
      l1 <- l[[1]]
      l2 <- l[[2]]
      l3 <- l[[3]]
      l4 <- l[[4]]
    }

    mu_star <- l1(as.matrix(mm) %*% as.vector(gam$mu.coefficients[1:(llm - 1)]) + sm)

    if (np > 1) {

      def <- c("~", "y", "1")

      if (setequal(c(as.character(est$input_var$f2)), def) == TRUE) {

        if (Iden == TRUE) {

          sigma_star <- stats::fitted(gam, "sigma")[1] * gamlss::getSmo(gam)$sige

        } else {

          sigma_star <- stats::fitted(gam, "sigma")[1]

        }

      } else {
        if (Iden == TRUE) {

          if (is.na(as.character(gam$sigma.coefficients)[lls]) == TRUE) {
            sigma_star <- l2(as.matrix(ms) %*% as.vector(gam$sigma.coefficients[1:(lls - 1)]) + ss) * gamlss::getSmo(gam)$sige
          } else {
            sigma_star <- l2(as.matrix(ms) %*% as.vector(gam$sigma.coefficients[1:(lls)]) + ss) * gamlss::getSmo(gam)$sige
          }
        } else {
          if (is.na(as.character(gam$sigma.coefficients)[lls]) == TRUE) {
            sigma_star <- l2(as.matrix(ms) %*% as.vector(gam$sigma.coefficients[1:(lls - 1)]) + ss)
          } else {
            sigma_star <- l2(as.matrix(ms) %*% as.vector(gam$sigma.coefficients[1:(lls)]) + ss)
          }
        }
      }
    }

    if (np > 2) {

      def <- c("~", "y", "1")

      if (setequal(c(as.character(est$input_var$f3)), def) == TRUE) {

        nu_star <- stats::fitted(gam, "nu")[1]

      } else {

        if (is.na(as.character(gam$nu.coefficients)[lln]) == TRUE) {
          nu_star <- l3(as.matrix(mn) %*% as.vector(gam$nu.coefficients[1:(lln - 1)]) + sn)
        } else {
          nu_star <- l3(as.matrix(mn) %*% as.vector(gam$nu.coefficients[1:(lln)]) + sn)
        }

      }
    }

    if (np > 3) {
      def <- c("~", "y", "1")
      if (setequal(c(as.character(est$input_var$f4)), def) == TRUE) {
        tau_star <- stats::fitted(gam, "tau")[1]
      } else {
        if (is.na(as.character(gam$tau.coefficients)[llt]) == TRUE) {
          tau_star <- l4(as.matrix(mt) %*% as.vector(gam$tau.coefficients[1:(llt - 1)]) + st)
        } else {
          tau_star <- l4(as.matrix(mt) %*% as.vector(gam$tau.coefficients[1:(llt)]) + st)
        }
      }
    }

    ys1 <- array()
    if (np == 1) {

      ys1 <- Dis(sum(Ni2), mu = mu_star)

    } else if (np == 2) {

      ys1 <- Dis(sum(Ni2), mu = mu_star, sigma = sigma_star)

    } else if (np == 3) {

      ys1 <- Dis(sum(Ni2), mu = mu_star, sigma = sigma_star, nu = nu_star)

    } else if (np == 4) {

      if (is.null(est$input_var$tau.fix) == TRUE && is.null(est$input_var$nu.fix) == TRUE) {

        ys1 <- Dis(sum(Ni2), mu = mu_star, sigma = sigma_star, nu = nu_star, tau = tau_star)

      } else if (is.null(est$input_var$tau.fix) == TRUE) {

        ys1 <- Dis(sum(Ni2), mu = mu_star, sigma = sigma_star, nu = 1, tau = tau_star)

      } else {

        ys1 <- Dis(sum(Ni2), mu = mu_star, sigma = sigma_star, nu = nu_star, tau = 1)

      }
    }

    data$ys <- as.vector(ys1) # add dataset
    media_boot1 <- data.frame("sa" = c(1:D), "ysm" = as.vector(by(data$ys, data$sa, mean)))
    media_boot2 <- data.frame("sa" = c(1:D), "ysm" = as.vector(by(data$ys, data$sa, povinc)))
    if (is.function(est$input_var$param)) media_boot3 <- data.frame("sa" = c(1:D), "ysm" = as.vector(by(data$ys, data$sa, est$input_var$param)))


    data <- subset(data, select = which(!duplicated(names(data))))
    sample_boot <- splitstackshape::stratified(data, "sa", samplesize)

    s <- as.data.frame(sample_boot)
    ns <- subset(data, !(data$id %in% s$id))
    estMSE <- est_saegamlss(
      sample = s, nonsample = ns, y_dip=est$input_var$y,
      sa = est$input_var$sa,
      f1 = replace_term(est$input_var$f1, quote(y), quote(ys)),
      f2 = replace_term(est$input_var$f2, quote(y), quote(ys)),
      f3 = replace_term(est$input_var$f3, quote(y), quote(ys)),
      f4 = replace_term(est$input_var$f4, quote(y), quote(ys)),
      fdis = est$input_var$fdis,
      R = est$input_var$R,
      Dis = est$input_var$Dis,
      param = est$input_var$param,
      tau.fix = as.numeric(est$input_var$tau.fix),
      nu.fix = as.numeric(est$input_var$nu.fix),
      z = as.numeric(est$input_var$z)
    )

    dif1 <- rep(0, D)
    dif2 <- rep(0, D)
    dif3 <- rep(0, D)

    for (t in 1:D) {


      if ( !is.function(param) ) if (param == "both" | param == "mean") dif1[t] <- (media_boot1$ysm[t] - estMSE$estimates$Mean[t])^2


      if ( !is.function(param) ) if (param == "both" | param == "HCR")  dif2[t] <- (media_boot2$ysm[t] - estMSE$estimates$HCR[t])^2

      if ( is.function(param) ) dif3[t] <- (media_boot3$ysm[t] - estMSE$est$Parameter[t])^2


    }

    if (gam$iter < max(est$input_var$nRS,est$input_var$nCG)) {
      for (h in 1:D) {
        MSE[h] <- dif1[h] + MSE[h]
        mse[h] <- dif2[h] + mse[h]
        mse_p[h] <- dif3[h] + mse_p[h]

      }
      u <- u + 1
    } else {
      u <- u
    }
    message("Processing mse loop ", u, " of ", loop)
  }

  MSE <- MSE / u
  mse <- mse / u
  mse_p <- mse_p / u

  if (!is.function(param)){
  if (est$input_var$param == "both") {

    est_mse <- data.frame("Est_mean" = est$estimates$Mean,
                          "Est_HCR" = est$estimates$HCR,
                          "MSE_mean" = MSE, "MSE_HCR" = mse,
                          "Est_mean_SD" =  sqrt(MSE),
                          "Est_HCR_SD" =  sqrt(mse),
                          "CV_mean" =  sqrt(MSE) / abs(est$estimates$Mean),
                          "CV_HCR" = sqrt(mse)/  abs(est$estimates$HCR)
                          )


  } else if (est$input_var$param == "mean") {

    est_mse <-  data.frame("Est_mean" = est$estimates$Mean,
                           "MSE_mean" = MSE,
                           "Est_mean_SD" =  sqrt(MSE),
                           "CV_mean" =  sqrt(MSE) / abs(est$estimates$Mean)
                           )

  } else {

    est_mse <-  data.frame("Est_HCR" = est$estimates$HCR,
                           "MSE_HCR" = mse,
                           "Est_HCR_SD" =  sqrt(mse),
                           "CV_HCR" = sqrt(mse)/  abs(est$estimates$HCR)
                           )

  }
    } else{

    est_mse <-  data.frame("Est_param" = est$estimates$Parameter,
                           "MSE_param" = mse_p,
                           "Est_param_SD" =  sqrt(mse_p),
                           "CV_param" =  sqrt(mse_p) / abs(est$estimates$Parameter)
                            )

  }



  #rm(fdis, envir = .GlobalEnv)
  result <- list("est_mse" = est_mse, "estimates" = est, "replicates"=loop)
  attr(result, "class") <- "saegamlss_class"

  return(result)
}
