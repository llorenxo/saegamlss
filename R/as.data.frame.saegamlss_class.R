#' as.data.frame
#'
#' @description transforms in data.frame the x of class "saegamlss".
#'
#' @param x any R object of class "saegamlss"
#' @param row.names NULL (Default) or a character vector giving the row names for the data frame. Missing values are not allowed
#' @param optional logical. If TRUE, setting row names and converting column names (to syntactic names: see make.names) is optional
#' @param n The population to be transformed in a data.frame. To be defined only for data_gen(). Default is 1
#' @param ... Other parameters
#'
#' @return A dataset for each population that has been generated
#' @seealso [data_gen()]
#' @export
#'
#' @examples
#'
#' data_p = data_gen(
#'   Ni = rep(10, 4), D = 4, ty = "no", k = 100, b1 = 4,
#'   x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
#'   x3 = NULL, b4 = NULL, x4 = NULL, xh = NULL,
#'   Dis=rNO, l = c(identity), sigma = 6, sigmah = NULL,
#'   sigmae = 22, costh = NULL, M = 2
#'   )
#'
#' as.data.frame(x = data_p)

as.data.frame.saegamlss <- function(x, row.names = NULL, optional = FALSE, n = 1, ...){

  if (names(x[1])=="estimates"){

   data <- data.frame("sa" = rownames(x$estimates) %>% as.factor())

   data <- as.data.frame(x = x$estimates, row.names = row.names, optional = optional)

  } else if (names(x[1])=="Gini" |
             names(x[1])=="Theil" |
             names(x[1])=="Atkinson"){


  data <- data.frame("sa" = names(x[[1]]) %>% as.factor())

  if (!is.null(x$Gini)) data <- data %>% dplyr::mutate("Gini" = x$Gini)

  if (!is.null(x$Theil)) data <- data %>% dplyr::mutate("Theil" = x$Theil)

  if (!is.null(x$Atkinson)) data <- data %>% dplyr::mutate("Atkinson" = x$Atkinson)


  data <- as.data.frame(x = data, row.names = row.names, optional = optional)

  }   else if ("Gini.MSE" %in% names(x[[1]]) |
               "Theil.MSE" %in% names(x[[1]])|
               "Atkinson.MSE" %in% names(x[[1]])){

    data <- data.frame("sa" = rownames(x$est_mse) %>% as.factor(),
                       x$est_mse)

    data <- as.data.frame(x = data, row.names = row.names, optional = optional)

  } else if (names(x[1])=="step1") {

    data <- data.frame(x$step3)

    data <- as.data.frame(x = data, row.names = row.names, optional = optional)


  } else if (names(x[[1]][1])=="P_Gini" |
             names(x[[1]][1])=="P_Theil" |
             names(x[[1]][1])=="P_Atkinson"){

    data <- data.frame("c"=1)

    if (!is.null(x$index$P_Gini)) data <- data %>% dplyr::mutate("Gini" = x$index$P_Gini)

    if (!is.null(x$index$P_Theil)) data <- data %>% dplyr::mutate("Theil" = x$index$P_Theil)

    if (!is.null(x$index$P_Atkinson)) data <- data %>% dplyr::mutate("Atkinson" = x$index$P_Atkinson)

    data <- data[,-1]

    data <- as.data.frame(x = data, row.names = row.names, optional = optional)


  } else if (names(x[1])=="dataset 1") {

    data <- as.data.frame(x = x[[n]], row.names = row.names, optional = optional)


  } else {


    data <- data.frame("sa" = rownames(x$est_mse), x$est_mse)
    data <- as.data.frame(x = data, row.names = row.names, optional = optional)

  }

  return(data)

  }



