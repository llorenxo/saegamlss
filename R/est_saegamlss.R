#' Monte Carlo estimation of mean and HCR based on SAE GAMLSS
#' @description A function to estimate the mean or/and the HCR     derived from unit-level small area estimation based on generalized additive models for location, scale and shape
#'
#' @param sample A dataset of sampled units.
#' @param nonsample  A dataset of non-sampled units. With covariates and small areas
#' @param y_dip The dependent variable
#' @param sa The name of the variable cotaining the areas
#' @param Ni 1xD vector containing the number of units (in population) for each area
#' @param ni 1xD vector containing the number of sampled-units for each area
#' @param f1 A formula object for fitting a model for the mu parameter, e.g.f1=y~x+random(sa)
#' @param f2 A formula object for fitting a model for the sigma parameter, e.g.f2=y~x+random(sa)
#' @param f3 A formula object for fitting a model for the nu parameter, e.g.f3=y~x+random(sa)
#' @param f4 A formula object for fitting a model for the tau parameter, e.g.f4=y~x+random(sa)
#' @param fdis The distribution family of the GAMLSS object
#' @param R 	Number of Monte-Carlo repetition. Default is 50
#' @param Dis Type of distribution in form of rDis where Dis is one of the distribution allowed by GAMLSS
#' @param np Number of parameters of the distribution. i.e. for the normal distribution np=2, for the GB2 distribution np=4
#' @param param The parameter to estimate, "mean" or "HCR". Default is both.
#' @param seed The seed
#' @param tau.fix A value to be fixed to 1 to obtain the Singh-Maddala distribution. Dist must be GB2 and np=4
#' @param nu.fix 	A value to be fixed to 1 to obtain the Dagum distribution. Dist must be GB2 and np=4
#' @param z 	The Poverty line. Default is equal to 0.6 of the dependent variable
#'
#' @return an object of class "saegamlss_class" which contains two list. The first list, named est has:
#'  ME 1xD vector with the estimate of the mean and
#'  HCR 1xD vector with the estimate of the HCR.
#'  The second, named input_var has:
#'  y the dependent variable
#'  sa the name of the variable cotaining the areas
#'  fit a GAMLSS object with all the components of the regression
#'  f1 a formula object used for fitting a model to the mu parameter
#'  f2 a formula object used for fitting a model to the sigma parameter
#'  f3 a  formula object used for fitting a model to the nu parameter
#'  f4 a formula object used for fitting a model to the tau parameter
#'  np the number of the parameters of the chosen distribution
#'  R the number of loop used to obtain the monte-carlo estimation
#'  D the number of area
#'  Ni 1xD vector containing the number of units (in population) for each area
#'  ni 1xD vector containing the number of sampled-units for each area
#'  originaldata The dataset of sampled-units
#'  fids The distribution used for the estimation
#'  z The poverty line
#' @export
#' @import gamlss
#' @import dplyr
#' @import splitstackshape
#' @examples # Generate data
#' #
#' dep.y <- data_gen(
#'   Ni = rep(10, 4), D = 4, M = 1, ty = "no", k = 1, b1 = 10,
#'   x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
#'   b4 = NULL, x4 = NULL, xh = NULL, Dis = NO,
#'   l = c(identity), sigma = 6, sigmah = NULL,
#'   sigmae = 2, costh = NULL, seed = 1234
#' )
#'
#' data <- dep.y[[1]]
#' #
#' # sample data with a sample fraction of 0.5
#' #
#' # sample data
#' #
#' sample <- stratified(data, "sa", size = 0.5)
#' # nonsample data
#' #
#' nonsample <- subset(data, !(data$id%in%sample$id))
#' # estimate
#' est_saegamlss(
#'   sample = sample, nonsample = nonsample, y_dip="y",
#'   sa="sa", Ni = rep(10, 4), ni = rep(1, 4),
#'   f1 = y ~ x1 + random(sa), f2 = NULL, f3 = NULL,
#'   f4 = NULL, fdis = NO, R = 200,
#'   Dis = rNO, np = 2, param = NULL,
#'   seed = 1234, tau.fix = NULL, nu.fix = NULL
#' )
#'
#' @references Mori, L., & Ferrante, M. R. (2023). Small area estimation under unit-level generalized additive models for location, scale and shape. arXiv e-prints, arXiv-2302.
#'  Graf, M., Marin, J. M., & Molina, I. (2019). A generalized mixed model for skewed distributions applied to small area estimation. Test, 28(2), 565–597.
#' @author Lorenzo Mori and Maria Rosaria Ferrante
#' @note With "object"$input_var$fit is possible to use all the classical function used by gamlss
#'  The Small Area have to be denoted with a number from 1 to D

est_saegamlss <- function(sample, nonsample, y_dip, sa, Ni, ni, f1, f2 = NULL,
                          f3 = NULL, f4 = NULL, fdis, R = NULL, Dis, np,
                          param = NULL, seed = NULL,
                          tau.fix = NULL, nu.fix = NULL, z = NULL) {

  y <- sample[[y_dip]]
  sa_n <- nonsample[[sa]]
  mixed <- NULL
  D <- length(Ni)

  if (is.null(seed)) seed <- 123

  if (is.null(param)) param <- c("both")


  if (is.null(R)) R <- 50

  #if (is.null(nRS)) nRS <- 150

  #if (is.null(nCG)) nCG <- 150

  if (is.null(f2)) f2 <- y ~ 1

  if (is.null(f3)) f3 <- y ~ 1

  if (is.null(f4)) f4 <- y ~ 1

  if (is.null(z)) z <- 0.6 * median(y, na.rm = TRUE)


  set.seed(seed)

  if (np == 1) {
    gam1 <- gamlss::gamlss(f1,
                           trace = F, family = substitute(fdis), data = sample, method = mixed(100, 100)
    )

    yns <- data.frame(
      "ys" = as.vector(gamlss::predictAll(gam1, data = sample, newdata = nonsample)$mu),
      "sa" = sa_n
    )

    ME <- rep(0, D)
    HCR <- rep(0, D)

    ranf <- gamlss::getSmo(gam1)$coef
    sige <- gamlss::getSmo(gam1)$sige
    sigb <- gamlss::getSmo(gam1)$sigb


    for (x in 1:R) {

      for (i in 1:D) {

        yns1 <- subset(yns, sa == i)
        lamb <- sigb / (sigb+sige/length(yns1$ys))
        yns1$ys <- yns1$ys - rep(ranf[i], length(yns1$ys))+lamb*(rep(ranf[i], length(yns1$ys)))
        samp <- as.vector(Dis(Ni[i] - ni[i], mu = yns1$ys))
        sp <- dplyr::pull(subset(sample, sa == i), f1[[2]])
        samp1 <- c(samp, sp)
        if (param == "both" | param == "mean") {
          ME[i] <- ME[i] + mean(samp1)
        }
        if (param == "both" | param == "HCR") {
          HCR[i] <- HCR[i] + povinc(samp1, z = z)
        }
      }
      message("Processing estimation loop ", x, " of ", R)
    }
    ME <- ME / x
    HCR <- HCR / x

  } else if (np == 2) {

    gam1 <- gamlss::gamlss(f1,
                           sigma.fo = f2, nu.fo = f3, tau.fo = f4,
                           trace = F, family = substitute(fdis), data = sample, method = mixed(100, 100)
    )

    yns <- data.frame(
      "ys" = as.vector(gamlss::predictAll(gam1, data = sample, newdata = nonsample)$mu),
      "yss" = as.vector(gamlss::predictAll(gam1, data = sample, newdata = nonsample)$sigma),
      "sa" = sa_n
    )
    ME <- rep(0, D)
    HCR <- rep(0, D)
    for (x in 1:R) {
      for (i in 1:D) {

        yns1 <- subset(yns, sa == i)
        samp <- as.vector(Dis(Ni[i] - ni[i], mu = yns1$ys, sigma = yns1$yss))
        sp <- dplyr::pull(subset(sample, sa == i), f1[[2]])
        samp1 <- c(samp, sp)
        if (param == "both" | param == "mean") {

          ME[i] <- ME[i] + mean(samp1)

        }
        if (param == "both" | param == "HCR") {

          HCR[i] <- HCR[i] + povinc(samp1, z = z)
        }

      }

      message("Processing estimation loop ", x, " of ", R)
    }

    ME <- ME / x
    HCR <- HCR / x

  } else if (np == 3) {

    gam1 <- gamlss::gamlss(f1,
                           sigma.fo = f2, nu.fo = f3, tau.fo = f4,
                           trace = F, family = substitute(fdis), data = sample, method = mixed(100, 100)
    )

    yns <- data.frame(
      "ys" = as.vector(gamlss::predictAll(gam1, data = sample, newdata = nonsample)$mu),
      "yss" = as.vector(gamlss::predictAll(gam1, data = sample, newdata = nonsample)$sigma),
      "ysn" = as.vector(gamlss::predictAll(gam1, data = sample, newdata = nonsample)$nu),
      "sa" = sa_n
    )
    ME <- rep(0, D)
    HCR <- rep(0, D)
    for (x in 1:R) {
      for (i in 1:D) {
        yns1 <- subset(yns, sa == i)
        samp <- as.vector(Dis(Ni[i] - ni[i], mu = yns1$ys, sigma = yns1$yss, nu = yns1$ysn))
        sp <- dplyr::pull(subset(sample, sa == i), f1[[2]])
        samp1 <- c(samp, sp)
        if (param == "both" | param == "mean") {
          ME[i] <- ME[i] + mean(samp1)
        }
        if (param == "both" | param == "HCR") {
          HCR[i] <- HCR[i] + povinc(samp1, z = z)
        }
      }
      message("Processing estimation loop ", x, " of ", R)
    }
    ME <- ME / x
    HCR <- HCR / x
  } else {
    gam1 <- gamlss::gamlss(f1,
                           sigma.fo = f2, nu.fo = f3, tau.fo = f4, tau.fixed = substitute(tau.fixed), nu.fixed = substitute(nu.fixed),
                           trace = F, family = substitute(fdis), data = sample, method = mixed(100, 100)
    )

    yns <- data.frame(
      "ys" = as.vector(gamlss::predictAll(gam1, data = sample, newdata = nonsample)$mu),
      "yss" = as.vector(gamlss::predictAll(gam1, data = sample, newdata = nonsample)$sigma),
      "ysn" = as.vector(gamlss::predictAll(gam1, data = sample, newdata = nonsample)$nu),
      "yst" = as.vector(gamlss::predictAll(gam1, data = sample, newdata = nonsample)$tau),
      "sa" = sa_n
    )

    ME <- rep(0, D)
    HCR <- rep(0, D)

    for (x in 1:R) {

      for (i in 1:D) {

        yns1 <- subset(yns, sa == i)
        samp <- as.vector(Dis(Ni[i] - ni[i], mu = yns1$ys, sigma = yns1$yss, nu = yns1$ysn, tau = yns1$yst))
        sp <- dplyr::pull(subset(sample, sa == i), f1[[2]])
        samp1 <- c(samp, sp)

        if (param == "both" | param == "mean") {

          ME[i] <- ME[i] + mean(samp1)

        }

        if (param == "both" | param == "HCR") {

          HCR[i] <- HCR[i] + povinc(samp1, z = z)

        }
      }

      message("Processing estimation loop ", x, " of ", R)

       }

     ME <- ME / x

     HCR <- HCR / x
  }

  if (param == "both") {

    estim <- list("ME" = ME, "HCR" = HCR)

  } else if (param == "mean") {

    estim <- list("ME" = ME)

    } else {

      estim <- list("HCR" = HCR)

      }

  input_var <- list(y=y_dip, sa=sa,
                    "fit" = gam1, "f1" = f1, "f2" = f2, "f3" = f3, "f4" = f4, "fdis" = fdis, Dis=Dis,
                    "np" = np, "nRS" = 100, "nCG" = 100, "R" = R, "D" = D, "Ni" = Ni, "ni" = ni,
                    "nu.fix" = nu.fix, "tau.fix" = tau.fix, "param" = param, "origindata" = sample, "z" = z

                    )

  result <- list("estimates" = estim, "input_var" = input_var)
  attr(result, "class") <- "saegamlss_class"
  return(result)
}
