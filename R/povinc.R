#' Head Count Ratio
#' @description Compute the Head Count Ration
#' @param y  A vector with numerical values. Usually income or consumption data.
#' @param z The poverty line. Default is 0.6*median(y)
#' @return The Head Count Ratio
#' @export
#'
#' @examples y <- runif(100, 0, 1000)
#' povinc(y)
#' @author Lorenzo Mori and Maria Rosaria Ferrante

povinc <- function(y, z = NULL) {
  if (is.null(z)) {
    z <- 0.6 * median(y)
  }

  result <- mean(y < z)
  return(result)
}
