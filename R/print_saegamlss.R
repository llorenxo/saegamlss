#' Print of x of class "saegamlss_class"
#' @description
#' Print for an x of class "saegamlss_class"
#'
#' @param x An x of class "saegamlss_class"
#' @param ... Additional parameters
#'
#' @return The Print of x of class "saegamlss_class"
#' @export
#' @examples
#'
#' #est_saegamlss
#'
#'dep.y <- data_gen(
#'   Ni = rep(10, 4), D = 4, M = 1, ty = "no", k = 1, b1 = 10,
#'   x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
#'   b4 = NULL, x4 = NULL, xh = NULL, Dis = NO,
#'   l = c(identity), sigma = 6, sigmah = NULL,
#'   sigmae = 2, costh = NULL
#' )
#'
#'print(dep.y)

print.saegamlss_class <- function(x, ...){
  if (names(x[1])=="estimates"){
    cat("SAE-GAMLSS based on:", "\n\n")
    print(x$input_var$fit)



  } else if (names(x[1])=="Gini" | names(x[1])=="Theil" | names(x[1])=="Atkinson"){
    print("ok")
  } else if (names(x[1])=="Gini.MSE"
             | names(x[1])=="Theil.MSE"
             | names(x[1])=="Atkinson.MSE"){
    print("ok")

  } else if (names(x[1])=="step1") {

    cat("\n\n",
        "Step 1 distribution selection:\n\n",
        x$step1, "\n\n",
        "Step 2 variable selection: \n\n")

    for (i in 1:length(x$step1)) {
      print(x$step1[i])

      if (as.character(eval(parse(text = deparse(x$step2[i])))[[1]][2])=="NULL"){
        mu_f <- as.character(eval(parse(text = deparse(x$step2[i])))[[1]][1])
        cat(paste("mu_f=", mu_f, "\n "))
      } else if (as.character(eval(parse(text = deparse(x$step2[i])))[[1]][3])=="NULL"){
        mu_f <- as.character(eval(parse(text = deparse(x$step2[i])))[[1]][1])
        sigma_f <- as.character(eval(parse(text = deparse(x$step2[i])))[[1]][2])
        cat(paste("mu_f=", mu_f, " sigma_f=", sigma_f, "\n "))
      } else if (as.character(eval(parse(text = deparse(x$step2[i])))[[1]][2])=="NULL"){
        mu_f <- as.character(eval(parse(text = deparse(x$step2[i])))[[1]][1])
        sigma_f <- as.character(eval(parse(text = deparse(x$step2[i])))[[1]][2])
        nu_f <- as.character(eval(parse(text = deparse(x$step2[i])))[[1]][3])
        cat(paste("mu_f=", mu_f, " sigma_f=", sigma_f, " nu_f=", nu_f, "\n "))
      } else {
        mu_f <- as.character(eval(parse(text = deparse(x$step2[i])))[[1]][1])
        sigma_f <- as.character(eval(parse(text = deparse(x$step2[i])))[[1]][2])
        nu_f <- as.character(eval(parse(text = deparse(x$step2[i])))[[1]][3])
        tau_f <- as.character(eval(parse(text = deparse(x$step2[i])))[[1]][4])
        cat(paste("mu_f=", mu_f, " sigma_f=", sigma_f, " nu_f=", nu_f, " tau_f=", tau_f, "\n "))
      }
    }

    cat("\n\n",
        "Step 3 k-fold cross-validation: \n\n")
    x$step3

  } else if (names(x[1])=="dataset 1") {
    cat("\n\n")
    for (i in 1:min(length(x), 5)) {
      print(paste("Print of population",i))
      print((as.data.frame(x[i])))
      if(length(x)>5) print(paste("Omitted", length(x)-5, "populations"))
    }

  } else {
    cat("The MSE values are based on:", x$replicates, "replicates", "\n\n")
    if (!is.null(x$est_mse$MSE_mean)) cat("The mean MSE for the Mean is:",mean(x$est_mse$MSE_mean), "\n\n")
    if (!is.null(x$est_mse$MSE_HCR)) cat("The mean MSE for the HCR is:",mean(x$est_mse$MSE_HCR), "\n\n")

  }
}

