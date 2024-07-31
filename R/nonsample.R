#Not Exp


nonsample <- function(data, sample, id) {
  nonsample <- subset(data, !(data$id %in% sample$id))
  return("nonsample" = as.data.frame(nonsample))
}
