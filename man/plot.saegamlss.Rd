% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_saegamlss.R
\name{plot.saegamlss}
\alias{plot.saegamlss}
\title{Plot function}
\usage{
\method{plot}{saegamlss}(
  x,
  compare.Gini = NULL,
  compare.Theil = NULL,
  compare.Atkinson = NULL,
  compare.Mean = NULL,
  compare.HCR = NULL,
  compare.param = NULL,
  xlim = c(0, 100),
  ylim = c(0, 1),
  ...
)
}
\arguments{
\item{x}{An object of class "saegamlss"}

\item{compare.Gini}{The MSE of a second estimator to be compared for the Gini index}

\item{compare.Theil}{The MSE of a second estimator to be compared for the Theil index}

\item{compare.Atkinson}{The MSE of a second estimator to be compared for the Atkinson index}

\item{compare.Mean}{The MSE of a second estimator to be compared for the Mean}

\item{compare.HCR}{The MSE of a second estimator to be compared for the HCR}

\item{compare.param}{The MSE of a second estimator to be compared for the self-defined parameter}

\item{xlim}{A vector with the x-axis limits}

\item{ylim}{A vector with the y-axis limits}

\item{...}{Additional parameters}
}
\value{
Return a plot for object of class "saegamlss"
}
\description{
Function to Plot object of class "saegamlss"
}
\examples{

#est_saegamlss

dep.y <- data_gen(
  Ni = rep(10, 4), D = 4, M = 1, ty = "no", k = 1, b1 = 10,
  x1 = rnorm(40, 0, 1), b2 = NULL, x2 = NULL, b3 = NULL,
  b4 = NULL, x4 = NULL, xh = NULL, Dis = NO,
  l = c(identity), sigma = 6, sigmah = NULL,
  sigmae = 2, costh = NULL
)

data <- dep.y[[1]]

# sample data with a sample fraction of 0.1

library(splitstackshape)

#sample data

sample <- stratified(data, "sa", size = 0.1)

# nonsample data
#
nonsample <- subset(data, !(data$id\%in\%sample$id))

#estimate

est <- est_saegamlss(
  sample = sample, nonsample = nonsample, y_dip="y",
  sa="sa", f1 = y ~ mu_x1 + random(sa),
  f2 = NULL, f3 = NULL,
  f4 = NULL, fdis = NO, R = 200,
  Dis = rNO, param = "Mean",
  tau.fix = NULL, nu.fix = NULL
)

plot(est)



}
