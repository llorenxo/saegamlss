#' Plot of object of class "saegamlss_class"
#' @description
#' Function to Plot object of class "saegamlss_class"
#'
#' @param x An object of class "saegamlss_class"
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
#' #est_saegamlss
#'
#'dep.y <- data_gen(
#'   Ni = rep(10, 4), D = 4, M = 1, ty = "no", k = 1, b1 = 10,
#'   x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
#'   b4 = NULL, x4 = NULL, xh = NULL, Dis = NO,
#'   l = c(identity), sigma = 6, sigmah = NULL,
#'   sigmae = 2, costh = NULL
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
#' nonsample <- subset(data, !(data$id%in%sample$id))
#' # estimate
#' est <- est_saegamlss(
#'   sample = sample, nonsample = nonsample, y_dip="y",
#'   sa="sa", Ni = rep(10, 4),
#'   f1 = y ~ x1 + random(sa), f2 = NULL, f3 = NULL,
#'   f4 = NULL, fdis = NO, R = 200,
#'   Dis = rNO, param = "Mean",
#'   tau.fix = NULL, nu.fix = NULL
#' )
#' plot(est)
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
#'   sample = sample, nonsample = nonsample, y_dip="y", sa="sa",
#'   Ni = rep(10, 4),  f1 = y ~ x1 + random(sa),
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
#'
#' plot(MSE, compare.Mean= runif(4,0, 2), compare.HCR= runif(4,2, 4))

plot.saegamlss_class <- function(x, compare.Gini = NULL,
                                 compare.Theil = NULL,
                                 compare.Atkinson = NULL,
                                 compare.Mean = NULL,
                                 compare.HCR = NULL, ...){
  if (names(x[1])=="estimates"){

    plot(x$input_var$fit)

  } else if (names(x[1])=="Gini" | names(x[1])=="Theil" | names(x[1])=="Atkinson"){

    plot(x$model)

  } else if (names(x[1])=="Gini.MSE" | names(x[1])=="Theil.MSE" | names(x[1])=="Atkinson.MSE"){

  if (is.null(compare.Gini) & is.null(compare.Theil) & is.null(compare.Atkinson)) stop(print("Error: an MSE to be compared is required"))

  Value <- Estimator <- NULL

  l <- length(x$Gini.MSE)

  if ( l == 0 ) l <- length(x$Theil.MSE)
  if ( l == 0 ) l <- length(x$Atkinson.MSE)

  pp <- pp1 <- pp2 <- pp3 <- pp4 <- pp5 <- rep(FALSE, l)

  if (length(sapply(x$Gini.MSE, function(x) is.null(x))) == l ) pp <- rep(TRUE, l)
  if (length(sapply(x$Theil.MSE, function(x) is.null(x))) == l ) pp1 <- rep(TRUE, l)
  if (length(sapply(x$Atkinson.MSE, function(x) is.null(x))) == l ) pp2 <- rep(TRUE, l)
  if (length(sapply(compare.Gini, function(x) is.null(x))) == l ) pp3 <- rep(TRUE, l)
  if (length(sapply(compare.Theil, function(x) is.null(x))) == l ) pp4 <- rep(TRUE, l)
  if (length(sapply(compare.Atkinson, function(x) is.null(x))) == l ) pp5 <- rep(TRUE, l)


  sub <- rep(NA, l)

  data <- data.frame(
      "Gini" = ifelse(!pp, sub , x$Gini.MSE),
      "Theil" = ifelse(!pp1, sub, x$Theil.MSE),
      "Atkinson" = ifelse(!pp2, sub, x$Atkinson.MSE),
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

  } else if (names(x[1])=="step1"){

    family <- as.character(x$step1[1])
    y <- x$y
    gamlss::histDist(y,family=family, data=x$sample_data)

  } else if (names(x[1])=="dataset 1"){

    for (i in 1:min(length(x), 5)) {
    nam <- colnames(x[[i]])

    graphics::par(mfrow=c(1,2))
    graphics::par(mar = c(5, 4, 4, 2) + 0.1)
    graphics::hist(x[[i]][,nam[1]],xlab="y", probability=TRUE, main="")

      corrplot::corrplot(stats::cor(x[[i]][,-c(length(nam)-1, length(nam))]),
                         method = "circle", type = "full", tl.cex = 0.7)

      graphics::mtext(paste("Histogram and Correlation Matrix - Population", i), side=3, line=-2, outer=TRUE)


      if (length(x)>1){
      cat("Press Enter for the next population")
      readline()
      }
      if (length(x)>5) print(paste("Omitted", length(x)-5, "populations"))
    }
    graphics::par(mfrow = c(1, 1))



  }  else {

    if (is.null(compare.Mean) & is.null(compare.HCR)) stop(print("Error: an MSE to be compared is required"))

    Value <- Estimator <- NULL

    l <- length(x$est_mse$MSE_mean)

    if ( l == 0 ) l <- length(x$est_mse$MSE_HCR)

    pp <- pp1 <- pp2 <- pp3 <- rep(FALSE, l)

    if (length(sapply(x$est_mse$MSE_mean, function(x) is.null(x))) == l ) pp <- rep(TRUE, l)
    if (length(sapply(x$est_mse$MSE_HCR, function(x) is.null(x))) == l ) pp1 <- rep(TRUE, l)
    if (length(sapply(compare.Mean, function(x) is.null(x))) == l ) pp2 <- rep(TRUE, l)
    if (length(sapply(compare.HCR, function(x) is.null(x))) == l ) pp3 <- rep(TRUE, l)


    sub <- rep(NA, l)

    data <- data.frame(
      "Mean" = ifelse(!pp, sub , x$est_mse$MSE_mean),
      "HCR" = ifelse(!pp1, sub, x$est_mse$MSE_HCR),
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


