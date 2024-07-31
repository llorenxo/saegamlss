#Not exp


povinc <- function(y, z = NULL) {
  if (is.null(z)) {
    z <- 0.6 * median(y)
  }

  result <- mean(y < z)
  return(result)
}
