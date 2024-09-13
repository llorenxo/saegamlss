#' Random effects
#'
#' @description Extract random effects for object of class "saegamlss"
#'
#'
#'
#' @param object An object of class saegamlss
#' @param parameters The parameter(S) for which the random effects has to be extracted. Default is c("mu", "sigma", "nu", "tau")
#' @param ... Other parameters
#'
#' @return A dataset with the estimated random effects. The first column "sa" denotes the Samll Area
#' @note The function ranef is firstly imported from nlme and than adapted for object of class "saegamlss"
#' @export
#'
#'
#' @examples
#' index_est <- sa_p_index(sample = s_data, y = "y",
#'                         sa = "sa", fdis = "LOGNO",
#'                         sigma.f = TRUE, index = "all")
#'
#' ranef(object = index_est, parameters = c("mu", "sigma"))
#'
#'

ranef.saegamlss <- function (object, parameters = c("mu", "sigma", "nu", "tau"), ...) {


  if (names(object[1])=="dataset 1"){

    print("Population generated with data_gen(): no random effects estimated")

  } else if (names(object[1])=="step1"){


    d = object$sample %>% dplyr::select(gsub(".*random\\(([^)]+)\\).*", "\\1",
                                             as.character(object$add_res$dist_step2[[1]]$mu.formula)[length( as.character(object$add_res$dist_step2[[1]]$mu.formula))]
    ))

    data <- data.frame("sa" = d %>% dplyr::distinct() )

    rownames(data) <- data$sa  %>% as.character()

    for (i in 1:length(object$add_res$dist_step2)) {

      if ("mu" %in% parameters & !is.null(object$add_res$dist_step2[[i]]$mu.coefficient) )   data <- data %>% dplyr::mutate(!!paste(names(object$add_res$dist_step2)[i], "mu_raneff", sep = "_") :=
                                                                 gamlss::getSmo(object$add_res$dist_step2[[i]], what = "mu")$coef
                                                               %>% as.vector())

      if ("sigma" %in% parameters & !is.null(object$add_res$dist_step2[[i]]$sigma.coefficient)) data <-  data %>% dplyr::mutate(!!paste(names(object$add_res$dist_step2)[i], "sigma_raneff", sep = "_") :=
                                                                     gamlss::getSmo(object$add_res$dist_step2[[i]], what = "sigma")$coef
                                                                   %>% as.vector())
      if ("nu" %in% parameters & !is.null(object$add_res$dist_step2[[i]]$nu.coefficient)) data <- data %>%  dplyr::mutate(!!paste(names(object$add_res$dist_step2)[i], "nu_raneff", sep = "_") :=
                                                                  gamlss::getSmo(object$add_res$dist_step2[[i]], what = "nu")$coef
                                                                %>% as.vector())
      if ("tau" %in% parameters & !is.null(object$add_res$dist_step2[[i]]$tau.coefficient)) data <-  data %>% dplyr::mutate(!!paste(names(object$add_res$dist_step2)[i], "tau_raneff", sep = "_") :=
                                                                   gamlss::getSmo(object$add_res$dist_step2[[i]], what = "tau")$coef
                                                                 %>% as.vector())


    }

    return (ran_eff = data)


  } else if (names(object[[1]][1])=="P_Gini" |
             names(object[[1]][1])=="P_Theil" |
             names(object[[1]][1])=="P_Atkinson"){

  print("Parametric inequality indices obtained with p_index(): no random effects estimated")

    } else{

  data <- data.frame("sa" = c(1:length(object$estimates[[1]])))
  rownames(data) <- rownames(object$estimates)

  results <- object %>% as.list() %>% purrr::pluck("input_var")


  if ("mu" %in% parameters & !is.null(results$fit$mu.coefficient)) data <- data %>% dplyr::mutate("mu_raneff" = gamlss::getSmo(results$fit, what = "mu")$coef %>% as.vector())
  if ("sigma" %in% parameters  & !is.null(results$fit$sigma.coefficient)) data <-  data %>% dplyr::mutate("sigma_raneff" = gamlss::getSmo(results$fit, what = "sigma")$coef  %>% as.vector())
  if ("nu" %in% parameters  & !is.null(results$fit$nu.coefficient)) data <- data %>%  dplyr::mutate("nu_raneff" = gamlss::getSmo(results$fit, what = "nu")$coef %>% as.vector())
  if ("tau" %in% parameters  & !is.null(results$fit$tau.coefficient) ) data <-  data %>% dplyr::mutate("tau_raneff" = gamlss::getSmo(results$fit, what = "tau")$coef %>% as.vector())

  return (ran_eff = data)

  }
}


