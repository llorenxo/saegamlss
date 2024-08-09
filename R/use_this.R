##################################
###Minor not exported functions###
##################################

# Replace term in functions
#
# @param f1 a function
# @param old an old term to be replaced
# @param new the new term
#
# @return a function

replace_term <- function(f1, old, new) {
  n <- length(f1)
  if (n > 1) {
    for (i in 1:n) f1[[i]] <- Recall(f1[[i]], old, new)

    return(f1)
  }

  if (f1 == old) new else f1
}

# is.NullOb
#
# @param x a list
#
# @return a list

is.NullOb <- function(x) {
  is.null(x) | all(sapply(x, is.null))
}

# Recursively step down into list, removing all such objects
#
# @param x a list
#
# @return a list

rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}


## Count argoment for np in est_saegamlss
# @param func a function

count_arguments <- function(func) {
  if (!is.function(func)) {
    stop("Input must be a function")
  }
  length(formals(func))
}





