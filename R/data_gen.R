#' Data generation
#'
#' @description A tool to generate a population with a dependent variable defined in terms of a given distribution and a set of (fixed) covariates.
#' The dependent variable can be generated for more than 100 distributions, each with a maximum of 4 parameters. Each parameter can be defined in terms of covariates and random effects.
#' When the link function is the identity function, the dependent variable may include an additional error term. For the normal distribution, it is possible to add heteroscedasticity.
#'
#' @param D Number of areas
#' @param Ni 1xD vector containing the number of units (in population) for each area
#' @param k A vector, maximum length equal to 4, of constants (intercepts) used to generate data, i.e. mu=k+b1*x
#' @param b1 A vector of regression coefficients used for mu
#' @param x1 A matrix/vector or a dataset where each column is a covariate used for mu
#' @param b2 A vector of regression coefficients used for sigma
#' @param x2 A matrix/vector or a dataset where each column is a covariate used for sigma
#' @param b3 A vector of regression coefficients used for nu
#' @param x3 A matrix/vector or a dataset where each column is a covariate used for nu
#' @param b4 A vector of regression coefficients used for tau
#' @param x4 A matrix/vector or a dataset where each column is a covariate used for tau
#' @param xh A vector of covariate used to generate heteroschedasticity
#' @param Dis Type of distribution in form of rDis where Dis is one of the distribution allowed by GAMLSS.
#' @param l A vector,  of maximum length equal to 4, in which each space is a function (the inverse of the link-function) to be applied to the corresponding parameter, i.e. mu, sigma, nu and tau
#' @param sigma A vector,  of maximum length equal to 4, of standard deviation used to generate random effects
#' @param sigmah  The standard deviation used to generate heteroschedasticity
#' @param sigmae The standard deviation used to generate the error term when the link-function is the identity one
#' @param ty To be specified equal to "no" if the generation process is mu=k+b1*x+e  or equal to "het" if the generation process have to include heteroschedasticity
#' @param costh A constant used to generate heteroschedasticity
#' @param id Unit id
#' @param M Number of replicates. Default is 1
#' @param seed The seed. Default is 123
#'
#' @import gamlss
#' @import dplyr
#' @import splitstackshape
#' @return An object of class "saegamlss" of length M. Each element of this list is a data frame that has: the dependent variable y, a set of fixed covariates, the small area column (sa) and the id column (id)
#' @export
#' @note The definition of the heteroschedastic term follows Ramirez-Aldana, R., & Naranjo, L. (2021).
#' @examples # Normal data (2 populations)
#'
#' data_gen(
#'   Ni = rep(10, 4), D = 4, ty = "no", k = 100, b1 = 4,
#'   x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
#'   x3 = NULL, b4 = NULL, x4 = NULL, xh = NULL,
#'   Dis=rNO, l = c(identity), sigma = 6, sigmah = NULL,
#'   sigmae = 22, costh = NULL, M = 2
#' )
#'
#' # Heteroschedastic Normal data
#'
#' data_gen(
#'   Ni = rep(10, 4), D = 4, ty = "het", k = 100, b1 = 4,
#'   x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
#'   x3 = NULL, b4 = NULL, x4 = NULL, xh = rnorm(40, 0, 1),
#'   Dis=rNO, l = c(identity), sigma = 6, sigmah = 1, sigmae = 22,
#'   costh = 0.1
#' )
#'
#' @references Mori, L., & Ferrante, M. R. (2024). Small area estimation under unit-level generalized additive models for location, scale and shape. Journal of Survey Statistics and Methodology, 1–37
#'  Ramirez-Aldana, R., & Naranjo, L. (2021). Random intercept and linear mixed models including heteroscedasticity in a logarithmic scale: Correction terms and prediction in the original scale. PloS one, 16(4), e0249910.
#'  Rigby, R. A., & Stasinopoulos, D. M. (2005). Generalized additive models for location, scale and shape. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(3), 507–554.
#'
#' @author Lorenzo Mori and Maria Rosaria Ferrante

data_gen <- function(Ni, D, k, b1, x1, b2 = NULL, x2 = NULL, b3 = NULL, x3 = NULL, b4 = NULL, x4 = NULL,
                     xh = NULL, Dis, l, sigma = NULL, sigmah = NULL, sigmae = NULL, ty = NULL, costh = NULL,
                     id = NULL, M = NULL, seed = 123) {

  set.seed(seed)

  list.data <- list()


  if (is.null(M))     M <- 1

  if (is.null(ty))  ty <- "other"

  if (is.null(costh) )  costh <- 1

  if (is.null(id)) id <- c(1:sum(Ni))

  if (!is.null(x1))   x1 <- x1 %>% as.data.frame() %>% stats::setNames(paste0("mu_x", 1:ncol(.)))

  if (!is.null(x2))  x2 <- x2 %>% as.data.frame() %>% stats::setNames(paste0("sigma_x", 1:ncol(.)))

  if (!is.null(x3))  x3 <- x3 %>% as.data.frame() %>% stats::setNames(paste0("nu_x", 1:ncol(.)))

  if (!is.null(x4))  x4 <- x4 %>% as.data.frame() %>% stats::setNames(paste0("tau_x", 1:ncol(.)))


  if (all(is.null(c(x2, x3, x4)))) {

    covariate <- as.data.frame(x1)


  } else if (all(is.null(c(x3, x4)))) {



    covariate <- as.data.frame(cbind(x1, x2))


  } else if (is.null(x4)) {



    covariate <- as.data.frame(cbind(x1, x2, x3))

  } else {


    covariate <- as.data.frame(cbind(x1, x2, x3, x4))

  }

  if (isTRUE(is.null(x2)))  x2 <- rep(0, sum(Ni))

  if (isTRUE(is.null(x3)))  x3 <- rep(0, sum(Ni))

  if (isTRUE(is.null(x4)))  x4 <- rep(0, sum(Ni))

  ym <- matrix(nrow = sum(Ni), ncol = M)

  for (i in 1:M) {

    message("Generating ", i, " of ", M)
    if (ty == "no") {
      sigma1 <- sigma[1]

      u1 <- rnorm(D, 0, sigma1) # re x


      rep1 <- rep(u1, Ni)

      repe <- rnorm(sum(Ni), 0, sigmae)


      b1 <- as.matrix(b1)
      x1 <- as.matrix(x1)

      ym[, i] <- k[1] + t(b1) %*% t(x1) + rep1 + repe

    } else if (ty == "het") {

      sigma1 <- sigma[1]

      u1 <- rnorm(D, 0, sigma1) # re x

      rep1 <- rep(u1, Ni)


      repe <- rnorm(sum(Ni), 0, sigmae)


      reph <- rep(rnorm(D, 0, sigmah), Ni)

      wj <- reph + xh * costh
      repw <- wj * repe

      ym[, i] <- k[1] + t(b1) %*% t(x1) + rep1 + repw

    } else {

      if (length(l) == 2) {
        l1 <- l[[1]]
        l2 <- l[[2]]


        sigma1 <- sigma[1]
        sigma2 <- sigma[2]


        u1 <- rnorm(D, 0, sigma1) # re x
        u2 <- rnorm(D, 0, sigma2) # re sig

        rep1 <- rep(u1, Ni)
        rep2 <- rep(u2, Ni)


        b1 <- as.matrix(b1)
        x1 <- as.matrix(x1)

        b2 <- as.matrix(b2)
        x2 <- as.matrix(x2)

        mu <- as.vector(l1(k[1] + t(b1) %*% t(x1) + rep1))
        sigmad <- as.vector(l2(k[2] + t(b2) %*% t(x2) + rep2))


        ym[, i] <- Dis(n = sum(Ni), mu, sigmad)

      } else if (length(l) == 3) {

        l1 <- l[[1]]
        l2 <- l[[2]]
        l3 <- l[[3]]


        sigma1 <- sigma[1]
        sigma2 <- sigma[2]
        sigma3 <- sigma[3]


        u1 <- rnorm(D, 0, sigma1) # re x
        u2 <- rnorm(D, 0, sigma2) # re sig
        u3 <- rnorm(D, 0, sigma3) # re x


        rep1 <- rep(u1, Ni)
        rep2 <- rep(u2, Ni)
        rep3 <- rep(u3, Ni)



        b1 <- as.matrix(b1)
        x1 <- as.matrix(x1)

        b2 <- as.matrix(b2)
        x2 <- as.matrix(x2)

        b3 <- as.matrix(b3)
        x3 <- as.matrix(x3)

        mu <- as.vector(l1(k[1] + t(b1) %*% t(x1) + rep1))
        sigmad <- as.vector(l2(k[2] + t(b2) %*% t(x2) + rep2))
        nu <- as.vector(l3(k[3] + t(b3) %*% t(x3) + rep3))

        ym[, i] <- Dis(sum(Ni), mu, sigmad, nu)

      } else {

        l1 <- l[[1]]
        l2 <- l[[2]]
        l3 <- l[[3]]
        l4 <- l[[4]]

        sigma1 <- sigma[1]
        sigma2 <- sigma[2]
        sigma3 <- sigma[3]
        sigma4 <- sigma[4]

        u1 <- rnorm(D, 0, sigma1)
        u2 <- rnorm(D, 0, sigma2)
        u3 <- rnorm(D, 0, sigma3)
        u4 <- rnorm(D, 0, sigma4)

        rep1 <- rep(u1, Ni)
        rep2 <- rep(u2, Ni)
        rep3 <- rep(u3, Ni)
        rep4 <- rep(u4, Ni)



        b1 <- as.matrix(b1)
        x1 <- as.matrix(x1)

        b2 <- as.matrix(b2)
        x2 <- as.matrix(x2)

        b3 <- as.matrix(b3)
        x3 <- as.matrix(x3)


        b4 <- as.matrix(b4)
        x4 <- as.matrix(x4)

        mu <- as.vector(l1(k[1] + t(b1) %*% t(x1) + rep1))
        sigmad <- as.vector(l2(k[2] + t(b2) %*% t(x2) + rep2))
        nu <- as.vector(l3(k[3] + t(b3) %*% t(x3) + rep3))
        tau <- as.vector(l4(k[4] + t(b4) %*% t(x4) + rep4))

        ym[, i] <- Dis(n = sum(Ni), mu, sigmad, nu, tau)
      }
    }

    data <- data.frame(
      "y" = ym[, i], covariate, "sa" = as.factor(sort(rep(seq(1:D), Ni))),
      "id" = id
    )

    list.data[[i]] <- data
    names(list.data)[i] <- paste("dataset", i)
  }

  attr(list.data, "class") <- "saegamlss"

  return(list.data)
}


