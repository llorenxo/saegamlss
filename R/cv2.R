#' Create the non-sample dataset
#' @description CV2
#'
#' @param est_MSE coefficient fo variation
#'
#' @return CV2 data frame with the original cv and cv2
#'
#' @export
#'
#' @author Lorenzo Mori and Maria Rosaria Ferrante
#'
#' @example
#'
#' index_est <- sa_p_index(sample = s_data, y = "y",
#'                         sa = "sa", fdis = "LOGNO",
#'                         sigma.f = TRUE, index = "all")
#'
#'
#' np <- np_mse(est = index_est, ncomp = "ncomp", R = 2)
#'
#'
#' @references Coefficient of variation: the second-order alternative (Tarald O. Kvålseth, 2017)
#'
#'
cv2= function(est_MSE){

 if ( class(est_MSE) != "saegamlss_class") error("est_MSE has to be of class saegamlss_class")



    cv2 =  est_MSE$est_mse %>%
      dplyr::mutate(dplyr::across(dplyr::ends_with("cv"),
                                  .fns = ~ (.^2 / (1 + .^2))^0.5,
                                  .names = "{.col}2"))


  attr(result, "class") <- "saegamlss_class"

 return(cv2)

}
