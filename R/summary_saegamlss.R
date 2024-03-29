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
#'  dep.y <- data_gen(
#'   Ni = rep(10, 4), D = 4, M = 1, ty = "no", k = 1, b1 = 10,
#'   x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
#'   b4 = NULL, x4 = NULL, xh = NULL, Dis = NO,
#'   l = c(identity), sigma = 6, sigmah = NULL,
#'   sigmae = 2, costh = NULL, seed = 1234
#' )
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
#' nonsample <- nonsample(data = data, sample = sample, id = id)
#' # estimate
#' est <- est_saegamlss(
#'   sample = sample, nonsample = nonsample,
#'   D = 4, Ni = rep(10, 4), ni = rep(1, 4),
#'   f1 = y ~ x1 + random(sa), f2 = NULL, f3 = NULL,
#'   f4 = NULL, fdis = NO, nRS = 150, nCG = 150, R = 200,
#'   Dis = rNO, np = 2, param = NULL,
#'   seed = 1234, tau.fix = NULL, nu.fix = NULL
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
#'   sigmae = 22, costh = NULL, seed = 1234
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
#' nonsample <- nonsample(data = data, sample = sample, id = id)
#' # estimate
#' est <- est_saegamlss(
#'   sample = sample, nonsample = nonsample, D = 4,
#'   Ni = rep(10, 4), ni = rep(1, 4), f1 = y ~ x1 + random(sa),
#'   f2 = NULL, f3 = NULL, f4 = NULL, fdis = NO, nRS = 150,
#'   nCG = 150, R = 2, Dis = rNO, np = 2, param = NULL,
#'   seed = 1234, tau.fix = NULL, nu.fix = NULL
#' )
#' #
#' # covariates
#' #
#' x <- data.frame(rep(1, nrow(data)), "x1" = data$x1)
#' #
#' # compute the MSE
#' #
#' MSE <- mse_saegamlss(
#'   est = est, D = 4, Ni = rep(10, 4), loop = 2,
#'   l = c(identity), Dis = rNO, Iden = TRUE,
#'   samplesize = 0.1, data = data, cov1 = x, cov2 = NULL, cov3 = NULL,
#'   cov4 = NULL, seed = 1234
#' )
#'
#' summary(MSE)

summary.saegamlss_class <- function(object, ...){
  if (names(object[1])=="est"){
    summary(object$input_var$fit)

    cat("\n\n",
        "Estimated random-effect(s): \n\n",
        getSmo(object$input_var$fit)$coef, "\n\n",

        "Estimated values of the Mean: \n\n",

        if (!is.null(object$est$ME)){
          object$est$ME
        }else {
          print("Not estimated")
        }, "\n\n",
        "Estimated values of the HCR: \n\n",

        if (!is.null(object$est$HCR)){
          object$est$HCR
        } else {
          print("Not estimated")
        }, "\n\n")

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
  } else if (names(object[1])=="Gini.MSE" | names(object[1])=="Theil.MSE" | names(object[1])=="Atkinson.MSE"){

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
  }   else {
    summary(object$est$input_var$fit)

    cat("\n\n",

        "Estimated random-effect(s): \n\n",
        getSmo(object$est$input_var$fit)$coef, "\n\n",

        "Summary of estimated MSE (mean): \n\n",

        if (!is.null(object$est_mse$MSE_mean)){
          names(summary(object$est_mse$MSE_mean))
          summary(object$est_mse$MSE_mean)
        }else {
          print("Not estimated")
        }, "\n\n",
        "Summary of estimated MSE (HCR): \n\n",

        if (!is.null(object$est_mse$MSE_HCR)){
          names(summary(object$est_mse$MSE_HCR))
          summary(object$est_mse$MSE_HCR)
        } else {
          print("Not estimated")
        }, "\n\n")

  }
}

