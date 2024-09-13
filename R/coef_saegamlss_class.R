#' Estimated regression parameters
#'
#' @description Extract Estimated regression parameters for object of class "saegamlss"
#'
#'
#' @param object An object of class saegamlss
#' @param parameters The parameter(S) for which the Estimated regression parameters has to be extracted. Default is c("mu", "sigma", "nu", "tau")
#' @param ... Other parameters
#'
#' @return A dataset with the estimated Estimated regression parameters
#' @note The function ranef is firstly imported from nlme and than adapted for object of class "saegamlss"
#' @export
#'
#'
#' @examples
#' index_est <- sa_p_index(sample = s_data, y = "y",
#'                         sa = "sa", fdis = "LOGNO",
#'                         sigma.f = TRUE, index = "all")
#'
#' coef(object = index_est, parameters = c("mu", "sigma"))
#'
#' est <- est_saegamlss(
#'   sample = s_data, nonsample = pop_data, y_dip="y",
#'   sa="sa", f1 = y ~ x1 + random(sa), f2 = y ~ x2 + random(sa),
#'   f3 = NULL, f4 = NULL, fdis = NO, R = 2,
#'   Dis = rNO, param = "both",
#'   tau.fix = NULL, nu.fix = NULL
#' )
#'
#' coef(object = est, parameters = c("mu", "sigma"))
#'

coef.saegamlss <- function (object, parameters = c("mu", "sigma", "nu", "tau"), ...) {

  results <- object %>% as.list() %>% purrr::pluck("input_var")


  if ("mu" %in% parameters) data_mu <- data.frame("Variables" = gsub("\\(([^)]+)\\)", "\\1", names(results$fit$mu.coefficients)), "Est_mu" = results$fit$mu.coefficients) %>% na.omit()
  if ("sigma" %in% parameters) data_sigma <- data.frame("Variables" = gsub("\\(([^)]+)\\)", "\\1", names(results$fit$sigma.coefficients)), "Est_sigma" = results$fit$sigma.coefficients) %>% na.omit()
  if ("nu" %in% parameters)  data_nu <- data.frame("Variables" = gsub("\\(([^)]+)\\)", "\\1", names(results$fit$nu.coefficients)), "Est_nu" = results$fit$nu.coefficients) %>% na.omit()
  if ("tau" %in% parameters)  data_tau <- data.frame("Variables" = gsub("\\(([^)]+)\\)", "\\1", names(results$fit$tau.coefficients)), "Est_tau" = results$fit$tau.coefficients) %>% na.omit()

  if (exists("data_mu")) {
    final_data <- data_mu
  }
  if (exists("data_sigma")) {
    final_data <- if (!is.null(final_data)) {
      full_join(final_data, data_sigma, by = "Variables")
    } else {
      data_sigma
    }
  }
  if (exists("data_nu")) {
    final_data <- if (!is.null(final_data)) {
      full_join(final_data, data_nu, by = "Variables")
    } else {
      data_nu
    }
  }
  if (exists("data_tau")) {
    final_data <- if (!is.null(final_data)) {
      full_join(final_data, data_tau, by = "Variables")
    } else {
      data_tau
    }
  }
  return (ran_eff = final_data)
}


