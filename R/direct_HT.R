
#' Direct Horvitz-Thompson estimatorof the mean
#'
#' @description
#' Compute the the Horvitz-Thompson Estimator for a finite population mean/proportion based on sample data
#' collected from a complex sampling design for Small Area.
#'
#' @param sample The data-set of sampled units containing y, sa, and pi
#' @param y The name of the responde variable
#' @param sa The name of the Small Areas
#' @param pi The name of the first order inclusion probability
#' @param pi2 A square matrix of the joint inclusion probabilities. Needed for the "LinHT" variance estimator
#' @param N A numeric value of the population size for each area. If NULL, it is estimated with the sum of the inverse of the pis
#' @param var_est A logical indicating whether or not to compute a variance estimator. Default is TRUE
#' @param var_method The method to use when computing the variance estimator. Options are a Taylor linearized technique: "LinHB"= Hajek-Berger estimator, "LinHH" = Hansen-Hurwitz estimator, "LinHTSRS" = Horvitz-Thompson estimator under simple random sampling without replacement, and "LinHT" = Horvitz-Thompson estimator or a resampling technique: "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement. The default is "bootstrapSRS"
#' @param B The number of bootstrap samples if computing the bootstrap variance estimator. Default is 100
#' @param fpc Default to TRUE, logical for whether or not the variance calculation should include a finite population correction when calculating the "LinHTSRS" or the "SRSbootstrap" variance estimator
#'
#' @importFrom mase horvitzThompson
#'
#' @return A list with the Horvitz-Thompson estimate of the mean, the variance, the Standard Deviation and the Coefficient of Variation for each Area
#' @export
#'
#' @examples
#'
#' ##################
#' ###Using s_data###
#' ##################
#'
#' #Using the default bootstrap
#'
#' dir = direct_HT(y = "y", sa = "sa", pi = "pi", sample = s_data)
#'
#' dir$results_HT
#'
#' #Using HH linearization
#' dir = direct_HT(y = "y", sa = "sa", pi = "pi", sample = s_data,
#'                 var_method = "LinHH")
#'
#' dir$results_HT
#'




direct_HT <- function(sample, y, sa, pi, pi2 = NULL, N = NULL, var_est = TRUE, var_method = "bootstrapSRS",
                      B = 100, fpc=TRUE){

  if (var_method == "LinHT") print("With LinHT option the estimated variance is not guaranteed to be greater than 0.")


  sa_name <- sa
  sa <- unique.numeric_version(sample[[sa]])
  D <- length(unique.numeric_version(sa))
  if (is.null(N)) N = by(1/sample[[pi]], sample[[sa_name]], sum)

 for (i in 1:D) {

  df3 <- subset(sample, sample[[sa_name]] == sa[i])

  d <- mase::horvitzThompson(y = as.vector(df3[[y]]),
                            pi =  as.vector(df3[[pi]]),
                            pi2 = pi2,
                            N = as.numeric(N[i]),
                            var_est = var_est,
                            var_method = var_method,
                            B = B,
                            fpc = fpc,
                            message =FALSE)

  if(i==1){

    results <- c(d$pop_mean, d$pop_mean_var, sqrt(d$pop_mean_var),  sqrt(d$pop_mean_var)/d$pop_mean )

    } else{

    results <- rbind(results, c(d$pop_mean, d$pop_mean_var, sqrt(d$pop_mean_var),  sqrt(d$pop_mean_var)/abs(d$pop_mean )))

    }

  }

  rownames(results) <- as.character(sa); results  <- results[order(rownames(results)), ]

  colnames(results) <- c("Est_mean", "Est_var", "Est_SD", "CV")

  results = list("results_HT" = results, "type" = var_method, "B" = B)

  attr(results, "class") <- "saegamlss_class"

  return(results)

}
