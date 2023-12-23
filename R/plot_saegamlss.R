#' Plot of object of class "saegamlss_class
#'
#' @param obj An object of class "saegamlss_clas"
#' @param compare.Gini The MSE of a second estimator to be compared for the Gini index
#' @param compare.Theil The MSE of a second estimator to be compared for the Theil index
#' @param compare.Atkinson The MSE of a second estimator to be compared for the Atkinson index
#' @param compare.Mean The MSE of a second estimator to be compared for the Mean
#' @param compare.HCR The MSE of a second estimator to be compared for the HCR
#' @param ... Additional parameters
#'
#' @return Return a plot for object of class "saegamlss_class"
#' @export
#' @examples
#'
#' # est_saegamlss
#'
#'  dep.y <- data_gen(
#'   Ni = rep(10, 4), D = 4, M = 2, ty = "no", k = 1, b1 = 10,
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
#' plot(est)
#'
#' #mse_saegamlss
#'
#' dep.y <- data_gen(
#'   Ni = rep(10, 4), D = 4, M = 2, ty = "no", k = 4, b1 = 100,
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
#' plot(MSE, compare.Mean= runif(4,0, 2), compare.HCR= runif(4,2, 4))

plot.saegamlss_class <- function(obj, compare.Gini = NULL,
                                 compare.Theil = NULL,
                                 compare.Atkinson = NULL,
                                 compare.Mean = NULL,
                                 compare.HCR = NULL, ...){
  if (names(obj[1])=="est"){

    plot(obj$input_var$fit)

  } else if (names(obj[1])=="Gini" | names(obj[1])=="Theil" | names(obj[1])=="Atkinson"){

    plot(obj$model)

  } else if (names(obj[1])=="Gini.MSE" | names(obj[1])=="Theil.MSE" | names(obj[1])=="Atkinson.MSE"){

  if (is.null(compare.Gini) & is.null(compare.Theil) & is.null(compare.Atkinson)) stop(print("Error: an MSE to be compared is required"))

  Value <- Estimator <- NULL

  l <- length(obj$Gini.MSE)

  if ( l == 0 ) l <- length(obj$Theil.MSE)
  if ( l == 0 ) l <- length(obj$Atkinson.MSE)

  pp <- pp1 <- pp2 <- pp3 <- pp4 <- pp5 <- rep(FALSE, l)

  if (length(sapply(obj$Gini.MSE, function(x) is.null(x))) == l ) pp <- rep(TRUE, l)
  if (length(sapply(obj$Theil.MSE, function(x) is.null(x))) == l ) pp1 <- rep(TRUE, l)
  if (length(sapply(obj$Atkinson.MSE, function(x) is.null(x))) == l ) pp2 <- rep(TRUE, l)
  if (length(sapply(compare.Gini, function(x) is.null(x))) == l ) pp3 <- rep(TRUE, l)
  if (length(sapply(compare.Theil, function(x) is.null(x))) == l ) pp4 <- rep(TRUE, l)
  if (length(sapply(compare.Atkinson, function(x) is.null(x))) == l ) pp5 <- rep(TRUE, l)


  sub <- rep(NA, l)

  data <- data.frame(
      "Gini" = ifelse(!pp, sub , obj$Gini.MSE),
      "Theil" = ifelse(!pp1, sub, obj$Theil.MSE),
      "Atkinson" = ifelse(!pp2, sub, obj$Atkinson.MSE),
      "compared.Gini" = ifelse(!pp3, sub, compare.Gini),
      "compared.Theil" = ifelse(!pp4, sub, compare.Theil),
      "compared.Atkinson" = ifelse(!pp5, sub, compare.Atkinson)

    )

  data= data[, colnames(data)[apply(data, 2, function(x) any(!is.na(x)))]]


  sa = rep(1:l, each=ncol(data))

  data = as.data.frame(tidyr::pivot_longer(data, everything(), names_to = "Column", values_to = "Value"))

  colnames(data)[1]="Estimator"
  data$sa=sa
  ggplot2::ggplot(data,  ggplot2::aes(x= sa, y = Value, group = Estimator, color=Estimator)) +
  ggplot2::geom_line()

  }   else {

    if (is.null(compare.Mean) & is.null(compare.HCR)) stop(print("Error: an MSE to be compared is required"))

    Value <- Estimator <- NULL

    l <- length(obj$est_mse$MSE_mean)

    if ( l == 0 ) l <- length(obj$est_mse$MSE_HCR)

    pp <- pp1 <- pp2 <- pp3 <- rep(FALSE, l)

    if (length(sapply(obj$est_mse$MSE_mean, function(x) is.null(x))) == l ) pp <- rep(TRUE, l)
    if (length(sapply(obj$est_mse$MSE_HCR, function(x) is.null(x))) == l ) pp1 <- rep(TRUE, l)
    if (length(sapply(compare.Mean, function(x) is.null(x))) == l ) pp2 <- rep(TRUE, l)
    if (length(sapply(compare.HCR, function(x) is.null(x))) == l ) pp3 <- rep(TRUE, l)


    sub <- rep(NA, l)

    data <- data.frame(
      "Mean" = ifelse(!pp, sub , obj$est_mse$MSE_mean),
      "HCR" = ifelse(!pp1, sub, obj$est_mse$MSE_HCR),
      "compared.Mean" = ifelse(!pp2, sub, compare.Mean),
      "compared.HCR" = ifelse(!pp3, sub, compare.HCR)
    )

    data= data[, colnames(data)[apply(data, 2, function(x) any(!is.na(x)))]]


    sa = rep(1:l, each=ncol(data))

    data = as.data.frame(tidyr::pivot_longer(data, everything(), names_to = "Column", values_to = "Value"))

    colnames(data)[1]="Estimator"
    data$sa=sa
    ggplot2::ggplot(data,  ggplot2::aes(x= sa, y = Value, group = Estimator, color=Estimator)) +
      ggplot2::geom_line()

  }
}

