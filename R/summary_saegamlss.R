#' Summary of object of class "saegamlss_class"
#' @description
#' Summary for an object of class "saegamlss_class"
#'
#' @param object An object of class "saegamlss_class"
#' @param ... Additional parameters
#'
#' @return The summary of object of class "saegamlss_class"
#' @export
#' @examples
#'
#' #est_saegamlss
#'
#'dep.y <- data_gen(
#'   Ni = rep(10, 4), D = 4, M = 1, ty = "no", k = 1, b1 = 10,
#'   x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
#'   b4 = NULL, x4 = NULL, xh = NULL, Dis = NO,
#'   l = c(identity), sigma = 6, sigmah = NULL,
#'   sigmae = 2, costh = NULL,
#'   )
#'
#' data <- dep.y[[1]]
#' #
#' # sample data with a sample fraction of 0.1
#' #
#' library(splitstackshape)
#' # sample data
#' #
#' sample <- stratified(data, "sa", size = 0.1)
#' # nonsample data
#' #
#' nonsample <- subset(data, !(data$id%in%sample$id))
#' # estimate
#' est <- est_saegamlss(
#'   sample = sample, nonsample = nonsample,  y_dip="y",
#'   sa="sa", Ni = rep(10, 4),
#'   f1 = y ~ x1 + random(sa), f2 = NULL, f3 = NULL,
#'   f4 = NULL, fdis = NO, R = 200,
#'   Dis = rNO, param = "Mean",
#'   tau.fix = NULL, nu.fix = NULL
#' )
#' summary(est)
#'
#' #mse_saegamlss
#'
#' dep.y <- data_gen(
#'   Ni = rep(10, 4), D = 4, M = 1, ty = "no", k = 4, b1 = 100,
#'   x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
#'   x3 = NULL, b4 = NULL, x4 = NULL, xh = NULL,
#'   Dis, l = c(identity), sigma = 6, sigmah = NULL,
#'   sigmae = 22, costh = NULL
#' )
#' data <- dep.y[[1]]
#' #
#' # sample data with a sample fraction of 0.5
#' #
#' library(splitstackshape)
#' # sample data
#' #
#' sample <- stratified(data, "sa", size = 0.5)
#' # nonsample data
#' #
#' nonsample <- subset(data, !(data$id%in%sample$id))
#' # estimate
#' est <- est_saegamlss(
#'   sample = sample, nonsample = nonsample, y_dip="y",  sa="sa",
#'   Ni = rep(10, 4), f1 = y ~ x1 + random(sa),
#'   f2 = NULL, f3 = NULL, f4 = NULL, fdis = NO,
#'   R = 2, Dis = rNO, param = "Mean",
#'   tau.fix = NULL, nu.fix = NULL
#' )
#' #
#' # covariates
#' #
#' x <- data.frame(rep(1, nrow(data)), "x1" = data$x1)
#' #
#' # compute the MSE
#' #
#' MSE <- mse_saegamlss(
#'   est = est, loop = 2,
#'   l = c(identity), Iden = TRUE,
#'   data = data, cov1 = x, cov2 = NULL, cov3 = NULL,
#'   cov4 = NULL
#' )
#' summary(MSE)

summary.saegamlss_class <- function(object, ...){

  if (names(object[1])=="estimates"){

    summary(object$input_var$fit)

    cat("\n\n",

        "Estimated random-effect(s): \n\n")

        print(getSmo(object$input_var$fit)$coef)
        cat("\n\n")


        if (!is.null(object$estimates$ME)){

        cat("\n\n Estimated values of the Mean: \n\n")

        print(object$estimates$ME)

        cat("\n\n")

        }


        if (!is.null(object$estimates$HCR)){

          cat("Estimated values of the HCR: \n\n")

          print(object$estimates$HCR)
          cat("\n\n")

        }


        if (!is.null(object$estimates$Parameter)){
          cat("Estimated values of the self-defined parameter: \n\n")
         print(object$estimates$Parameter)

        }

  } else if (names(object[1])=="Gini" | names(object[1])=="Theil" | names(object[1])=="Atkinson"){

    summary(object$model)

    cat("\n\n",
        "Estimated random-effect(s): \n\n",
        getSmo(object$model)$coef, "\n\n",

        "Estimated values of the Gini: \n\n",

        if (!is.null(object$Gini)){
          object$Gini
        }else {
          print("Not estimated")
        }, "\n\n",
        "Estimated values of the Theil: \n\n",

        if (!is.null(object$Theil)){
          object$Theil
        } else {
          print("Not estimated")
        }, "\n\n",
        "Estimated values of the Atkinson: \n\n",

        if (!is.null(object$Atkinson)){
          object$Atkinson
        } else {
          print("Not estimated")
        }, "\n\n")

  } else if (names(object[[1]][1])=="P_Gini" |
             names(object[[1]][1])=="P_Theil" |
             names(object[[1]][1])=="P_Atkinson"){

    cat("Assumed distribution: \n\n")

    print(object[[2]])

    cat("\n\n")

    if("P_Gini" %in% names(object[[1]])){

      cat("Estimated value of Gini index: \n\n")

      print(object[[1]]$P_Gini)
      cat("\n\n")

      }

    if("P_Theil" %in% names(object[[1]])){
      cat("Estimated value of Theil index: \n\n")

      print(object[[1]]$P_Theil)
      cat("\n\n")

      }

    if("P_Atkinson" %in% names(object[[1]])){
      cat("Estimated value of Atkinson index: \n\n")

      print(object[[1]]$P_Atkinson)
      cat("\n\n")


      cat("Parameter: \n\n")


      print(object[[3]])



    }


  }   else if (names(object[1])=="Gini.MSE"|
               names(object[1])=="Theil.MSE"|
               names(object[1])=="Atkinson.MSE"){

    cat("\n\n",

        "Estimated MSE for Gini index: \n\n",

        if (!is.null(object$Gini.MSE)){
          object$Gini.MSE
        }else {
          print("Not estimated")
        }, "\n\n",
        "Estimated MSE for Theil index: \n\n",

        if (!is.null(object$Theil)){
          object$Theil.MSE
        } else {
          print("Not estimated")
        }, "\n\n",
        "Estimated MSE for Atkinson index: \n\n",

        if (!is.null(object$Atkinson)){
          object$Atkinson.MSE
        } else {
          print("Not estimated")
        }, "\n\n")
  } else if (names(object[1])=="step1") {

    cat("\n\n",
        "Step 1 distribution selection:\n\n",
        object$step1, "\n\n",
        "Step 2 variable selection: \n\n")

    for (i in 1:length(object$step1)) {
      print(object$step1[i])

      if (as.character(eval(parse(text = deparse(object$step2[i])))[[1]][2])=="NULL"){
        mu_f <- as.character(eval(parse(text = deparse(object$step2[i])))[[1]][1])
        cat(paste("mu_f=", mu_f, "\n "))
      } else if (as.character(eval(parse(text = deparse(object$step2[i])))[[1]][3])=="NULL"){
        mu_f <- as.character(eval(parse(text = deparse(object$step2[i])))[[1]][1])
        sigma_f <- as.character(eval(parse(text = deparse(object$step2[i])))[[1]][2])
        cat(paste("mu_f=", mu_f, " sigma_f=", sigma_f, "\n "))
      } else if (as.character(eval(parse(text = deparse(object$step2[i])))[[1]][2])=="NULL"){
        mu_f <- as.character(eval(parse(text = deparse(object$step2[i])))[[1]][1])
        sigma_f <- as.character(eval(parse(text = deparse(object$step2[i])))[[1]][2])
        nu_f <- as.character(eval(parse(text = deparse(object$step2[i])))[[1]][3])
        cat(paste("mu_f=", mu_f, " sigma_f=", sigma_f, " nu_f=", nu_f, "\n "))
      } else {
        mu_f <- as.character(eval(parse(text = deparse(object$step2[i])))[[1]][1])
        sigma_f <- as.character(eval(parse(text = deparse(object$step2[i])))[[1]][2])
        nu_f <- as.character(eval(parse(text = deparse(object$step2[i])))[[1]][3])
        tau_f <- as.character(eval(parse(text = deparse(object$step2[i])))[[1]][4])
        cat(paste("mu_f=", mu_f, " sigma_f=", sigma_f, " nu_f=", nu_f, " tau_f=", tau_f, "\n "))
      }
    }

    cat("\n\n",
        "Step 3 k-fold cross-validation: \n\n")
    object$step3

  } else if (names(object[1])=="dataset 1") {
    cat("\n\n")
    for (i in 1:length(object)) {
      print(paste("summary of population",i))

      print(summary((object[[i]])))
    }

  } else {

    cat("Summary of estimated MSE: \n\n")


        if (!is.null(object$est_mse$MSE_mean)){
          cat("Mean:  \n\n")
          print(summary(object$est_mse$MSE_mean))
          cat("\n\n")
        }

        if (!is.null(object$est_mse$MSE_HCR)){
          cat("HCR:  \n\n")
          print(summary(object$est_mse$MSE_HCR))
          cat("\n\n")
        }
    if (!is.null(object$est_mse$MSE_param)){
      cat("Self-defined parameter:  \n\n")
      print(summary(object$est_mse$MSE_param))
      cat("\n\n")
    }

  }
}

