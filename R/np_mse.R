#' Non-parametric MSE for simplified SAE-GAMLSS
#' @description Compute the Non-parametric MSE for simplified SAE-GAMLSS for three possible different
#' indicators (Gini, Theil, Atkinson)
#'
#' @param est An object of class "saegamlss_class" obtained with \code{sa_p_index}
#' @param ncomp The number of components of each household
#' @param R The number of loops to be performed. Default is 200
#'
#'
#'
#' @return An object of class "saegamlss_class" containing the values of the MSE for each area and for each index
#'
#' @author Lorenzo Mori and Maria Rosaria Ferrante
#' @export
#'
#' @examples
#'
#' ##################
#' ###Using s_data###
#' ##################
#'
#' index_est <- sa_p_index(data = s_data, y = "y",
#'                         sa = "sa", fdis = "LOGNO",
#'                         sigma.f = TRUE, index = "all")
#'
#'
#' np <- np_mse(est = index_est, ncomp = "ncomp", R = 2)
#'
#' np$Gini.MSE
#'


np_mse <- function(est, ncomp,  R = 200){


  set.seed(est$index_est$seed)

  data <- est$input_var$data
  y <- est$index_est$y
  sa <- est$index_est$sa
  ncomp <- data[[ncomp]]
  fdis <- est$input_var$fdis
  index <- est$input_var$index
  epsilon <- est$input_var$epsilon

  f1 <- y ~1 + random(as.factor(sa))
  f2 <- f3 <- f4 <- y ~1

  if(isTRUE(est$index_est$sigma.f))f2 <- f1
  if(isTRUE(est$index_est$nu.f))f3 <- f1
  if(isTRUE(est$index_est$tau.f))f4 <- f1

  if (is.null(est$index_est$w)) data <- data %>% mutate(w = rep(1:nrow(data)))



  s2 <- data %>% dplyr::select("y", "sa", "ncomp", "w")

  s2$id <- c(1:nrow(s2))

  l_sa <- length(data %>% dplyr::distinct(sa) %>% dplyr::pull())


  dif3 <- rep(0, l_sa)
  dif4 <- rep(0, l_sa)
  dif5 <- rep(0, l_sa)


for (t2 in 1:R){

  message("Processing mse loop ", t2, " of ", R)

  s_b <- s2[1,]
  ss <- s2

  for (i in 1:l_sa){

    s_b_sa1 <- subset(ss, sa==i)
    s_b_sa1$id <- c(1: nrow(s_b_sa1))

    if(nrow(s_b_sa1)>1 ){

      s_b_sa <- s_b_sa1[sample(1:nrow(s_b_sa1), nrow(s_b_sa1)-1, replace=TRUE),]

      } else {

      s_b_sa <- s_b_sa1[sample(1:nrow(s_b_sa1), nrow(s_b_sa1), replace=TRUE),]

      }

    for (j in 1:nrow(s_b_sa)){

      s_b_sa [[j,5]] <- s_b_sa[[j,5]] * nrow(s_b_sa1) / nrow(s_b_sa) * rep(as.vector(table(s_b_sa$id)),
                                                                 as.vector(table(s_b_sa$id)))[j]
    }

    s_b_sa$w < -s_b_sa$w * sum(s_b_sa1$w) / sum(s_b_sa$w)

    s_b<-rbind(s_b, s_b_sa)

  }

  s_b <- s_b[-1,]
  s1 <- s_b
  s1 <- as.data.frame(lapply(s1, rep, s1$ncomp))

  cont<-gamlss.control(c.crit = 0.01, n.cyc = 15, mu.step = 1, sigma.step = 1, nu.step = 1,
                      tau.step = 1, gd.tol = Inf, iter = 0, trace = F, autostep = TRUE)


  ga<-gamlss::gamlss(f1, sigma.fo=f2, nu.fo=f3, tau.fo=f4, weights = w,
                      trace = FALSE, family = substitute(fdis), data = as.data.frame(s1),
                      method = RS(100))




  predi <- predictAll(ga, data=s1, newdata=s1)

  if (!is.null(predi$mu)) s1$mu_d <- predi$mu
  if (!is.null(predi$sigma)) s1$sigma_d <- predi$sigma
  if (!is.null(predi$nu)) s1$nu_d <- predi$nu
  if (!is.null(predi$tau)) s1$tau_d <- predi$tau


  if (length(s1$mu_d)<1)   s1$mu_d=s1$mu_d[1:nrow(s1)]

  if (length(s1$nu_d)<1)   s1$nu_d=s1$nu_d[1:nrow(s1)]

  if (length(s1$tau_d)<1)  s1$tau_d=s1$tau_d[1:nrow(s1)]

  if (length(s1$sigma_d)<1) s1$sigma_d=s1$sigma_d[1:nrow(s1)]



  atkinson_gamlss_boot <- 0
  gini_gamlss_boot <- 0
  theil_gamlss_boot <- 0


  for (i in 1:l_sa) {

    a <- subset(data, sa==i)


    if (index == "all" | index == "Gini") {

      gini_gamlss_boot <- rbind(gini_gamlss_boot, p_index(mu=s1$mu_d[1], sigma=s1$sigma_d[1],
                                           nu=s1$nu_d[1], tau=s1$tau_d[1], fdis= fdis, index="Gini", epsilon = epsilon )$index$P_Gini)

     }

    if (index == "all" | index == "Theil"){

       theil_gamlss_boot <- rbind(theil_gamlss_boot, p_index(mu=s1$mu_d[1], sigma=s1$sigma_d[1],
                                             nu=s1$nu_d[1], tau=s1$tau_d[1], fdis=fdis, index="Theil", epsilon = epsilon )$index$P_Theil)

    }

    if (index == "all" | index == "Atkinson"){

    atkinson_gamlss_boot <- rbind(atkinson_gamlss_boot, p_index(mu=s1$mu_d[1], sigma=s1$sigma_d[1],
                                                   nu=s1$nu_d[1], tau=s1$tau_d[1], fdis=fdis, index="Atkinson", epsilon = epsilon )$index$P_Atkinson)
     }

  }

  dif3 <- rbind(dif3, atkinson_gamlss_boot[-1])
  dif4 <- rbind(dif4, gini_gamlss_boot[-1])
  dif5 <- rbind(dif5, theil_gamlss_boot[-1])

  }

  dif3 <- as.data.frame(lapply(as.data.frame(dif3), as.numeric))
  dif4 <- as.data.frame(lapply(as.data.frame(dif4), as.numeric))
  dif5 <- as.data.frame(lapply(as.data.frame(dif5), as.numeric))


  if (index == "all"){

  result <- data.frame("Gini.Est" = est$Gini,
                       "Gini.MSE" = colMeans(dif4[-1,]),
                       "Gini.SD" = sqrt(colMeans(dif4[-1,])),
                       "Gini.CV" = sqrt(colMeans(dif4[-1,]))/abs(est$Gini),
                       "Theil.Est" = est$Theil,
                       "Theil.MSE" = colMeans(dif5[-1,]),
                       "Theil.SD" = sqrt(colMeans(dif5[-1,])),
                       "Theil.CV" = sqrt(colMeans(dif5[-1,]))/abs(est$Theil),
                       "Atkinson.Est" = est$Atkinson,
                       "Atkinson.MSE" = colMeans(dif3[-1,]),
                       "Atkinson.SD" = sqrt(colMeans(dif3[-1,])),
                       "Atkinson.CV" = sqrt(colMeans(dif3[-1,]))/abs(est$Atkinson)
                       )

  } else if (index=="Gini"){

  result <- data.frame("Gini.Est" = est$Gini,
                       "Gini.MSE" = colMeans(dif4[-1,]),
                       "Gini.SD" = sqrt(colMeans(dif4[-1,])),
                       "Gini.CV" = sqrt(colMeans(dif4[-1,]))/abs(est$Gini),
                      )

  } else if (index=="Theil"){

  result <- data.frame("Theil.Est" = est$Theil,
                        "Theil.MSE" = colMeans(dif5[-1,]),
                        "Theil.SD" = sqrt(colMeans(dif5[-1,])),
                        "Theil.CV" = sqrt(colMeans(dif5[-1,]))/abs(est$Theil),
                        )

  } else {

  result <- data.frame("Atkinson.Est" = est$Atkinson,
                       "Atkinson.MSE" = colMeans(dif3[-1,]),
                       "Atkinson.SD" = sqrt(colMeans(dif3[-1,])),
                       "Atkinson.CV" = sqrt(colMeans(dif3[-1,]))/abs(est$Atkinson)
                        )

  }


  rownames(result) <- levels(sa)

  result = list("est_mse" = result, "input_var" = est$input_var,
                "model" = est$model, "estimates_par" = est$estimates_par,
                "R" = R)

  attr(result, "class") <- "saegamlss_class"
  return(result)

}
