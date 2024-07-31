#' Create the non-sample dataset
#' @description Subset the non-sample part of the population
#' @param data Population dataset
#' @param sample Sample dataset
#' @param id id variable, must be the same in both the datasets
#'
#' @return A dataset of non-sampled units

nonsample <- function(data, sample, id) {
  nonsample <- subset(data, !(data$id %in% sample$id))
  return("nonsample" = as.data.frame(nonsample))
}
