% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p_index.R
\name{p_index}
\alias{p_index}
\title{Distribution based Gini, Theil and Atkinson index estimation}
\usage{
p_index(
  mu,
  sigma = NULL,
  nu = NULL,
  tau = NULL,
  fdis,
  index = "all",
  epsilon = 1
)
}
\arguments{
\item{mu}{The estimated value of mu}

\item{sigma}{The estimated value of sigma}

\item{nu}{The estimated value of nu}

\item{tau}{The estimated value of tau}

\item{fdis}{The assumed distribution. Options are: GB2 (Generalized Beta of 2-type), GAMMA (Gamma), EXP (Exponential), LOGNO (Log-Normal), PARETO (Pareto), WEI (Weibull)}

\item{index}{The index to be estimated ("Gini", "Theil" or "Atkinson"). Default is all}

\item{epsilon}{The value for the poverty aversion parameter. Default is  1}
}
\value{
The estimated index (indicators)
}
\description{
Provide the estimated values for Gini, Theil, and Atkinson indices using the closed-form expression for an assumed distribution.
}
\examples{


#Using GB2 distribution

p_index(mu=5, sigma=2, nu=1, tau=1, fdis="GB2", index="Gini", epsilon=2)
p_index(mu=5, sigma=2, nu=1, tau=1, fdis="GB2", index="Theil", epsilon=2)
p_index(mu=5, sigma=2, nu=1, tau=1, fdis="GB2", index="Atkinson", epsilon=2)
p_index(mu=5, sigma=2, nu=1, tau=1, fdis="GB2", index="all", epsilon=2)

#Using Exponential distribution

p_index(mu = 10, sigma=0.8, fdis="EXP", index="Gini", epsilon=2)
p_index(mu = 10, sigma=0.8, fdis="EXP", index="Theil", epsilon=2)
p_index(mu = 10, sigma=0.8, fdis="EXP", index="Atkinson", epsilon=2)
p_index(mu = 10, sigma=0.8, fdis="EXP", index="all", epsilon=2)

}
\author{
Lorenzo Mori and Maria Rosaria Ferrante
}
