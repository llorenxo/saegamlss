#' Extract random effects for object of class "saegamlss_class"
#'
#' @param model An object of class saegamlss_class
#' @param parameters the distribution parameters from which random effects are extracted. Default is  c("mu", "sigma", "nu", "tau")
#'
#' @return A dataset with the estiamted random effects
#' @export
#'
#' @examples
#'
#' index_est <- sa_p_index(sample = s_data, y = "y",
#'                         sa = "sa", fdis = "LOGNO",
#'                         sigma.f = TRUE, index = "all")
#'
#' ran_eff(model = index_est, parameters = c("mu", "sigma"))
#'
ran_eff <- function (model, parameters = c("mu", "sigma", "nu", "tau")) {

  data <- data.frame("sa" = c(1:length(model$estimates[[1]])))
  rownames(data) <- rownames(model$estimates)

  results <- model %>% as.list() %>% purrr::pluck("input_var")


  if ("mu" %in% parameters) data <- data %>% dplyr::mutate("mu_raneff" = gamlss::getSmo(results$fit, what = "mu")$coef %>% as.vector())
  if ("sigma" %in% parameters) data <-  data %>% dplyr::mutate("sigma_raneff" = gamlss::getSmo(results$fit, what = "sigma")$coef  %>% as.vector())
  if ("nu" %in% parameters) data <- data %>%  dplyr::mutate("nu_raneff" = gamlss::getSmo(results$fit, what = "nu")$coef %>% as.vector())
  if ("tau" %in% parameters) data <-  data %>% dplyr::mutate("tau_raneff" = gamlss::getSmo(results$fit, what = "tau")$coef %>% as.vector())

  data = data[,-1]
  return (ran_eff = data)
}

