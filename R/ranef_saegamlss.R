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
#' @return A dataset with the estimated random effects
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

  data <- data.frame("sa" = c(1:length(object$estimates[[1]])))
  rownames(data) <- rownames(object$estimates)

  results <- object %>% as.list() %>% purrr::pluck("input_var")


  if ("mu" %in% parameters) data <- data %>% dplyr::mutate("mu_raneff" = gamlss::getSmo(results$fit, what = "mu")$coef %>% as.vector())
  if ("sigma" %in% parameters) data <-  data %>% dplyr::mutate("sigma_raneff" = gamlss::getSmo(results$fit, what = "sigma")$coef  %>% as.vector())
  if ("nu" %in% parameters) data <- data %>%  dplyr::mutate("nu_raneff" = gamlss::getSmo(results$fit, what = "nu")$coef %>% as.vector())
  if ("tau" %in% parameters) data <-  data %>% dplyr::mutate("tau_raneff" = gamlss::getSmo(results$fit, what = "tau")$coef %>% as.vector())

  return (ran_eff = data)
}


