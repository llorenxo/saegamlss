% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonsample.R
\name{nonsample}
\alias{nonsample}
\title{Create the non-sample dataset}
\usage{
nonsample(data, sample, id)
}
\arguments{
\item{data}{Population dataset}

\item{sample}{Sample dataset}

\item{id}{id variable, must be the same in both the datasets}
}
\value{
A dataset of non-sampled units
}
\description{
Subset the non-sample part of the population
}
\examples{
# sample data with a sample fraction of 0.1
library(splitstackshape)
data <- data.frame("y" = runif(100, 0, 100), "id" = c(1:100), "sa" = rep(1:5, 20))
# sample data
sample <- stratified(data, "sa", size = 0.1)
# nonsample data
nonsample <- nonsample(data = data, sample = sample, id = id)

}
\references{
Mori, L., & Ferrante, M. R. (2023). Small area estimation under unit-level generalized additive models for location, scale and shape. arXiv e-prints, arXiv-2302.
}
\author{
Lorenzo Mori and Maria Rosaria Ferrante
}
