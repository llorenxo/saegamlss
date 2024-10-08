% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/np_mse.R
\name{np_mse}
\alias{np_mse}
\title{Non-parametric MSE for simplified SAE-GAMLSS}
\usage{
np_mse(est, ncomp, R = 200)
}
\arguments{
\item{est}{An object of class "saegamlss" obtained with \code{sa_p_index}}

\item{ncomp}{The number of components of each household}

\item{R}{The number of loops to be performed. Default is 200}
}
\value{
An object of class "saegamlss" containing the values of the MSE for each area and for each index
}
\description{
Compute the Non-parametric MSE for simplified SAE-GAMLSS for three possible different
indicators (Gini, Theil, Atkinson).
}
\examples{

##################
###Using s_data###
##################

index_est <- sa_p_index(sample = s_data, y = "y",
                        sa = "sa", fdis = "LOGNO",
                        sigma.f = TRUE, index = "all")


np <- np_mse(est = index_est, ncomp = "ncomp", R = 2)

np$est_mse

}
\author{
Lorenzo Mori and Maria Rosaria Ferrante
}
