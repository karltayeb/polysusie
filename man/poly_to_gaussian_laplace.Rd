% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gauss_poly.R
\name{poly_to_gaussian_laplace}
\alias{poly_to_gaussian_laplace}
\title{Laplace approximation for polynomial density}
\usage{
poly_to_gaussian_laplace(p, lower = -10, upper = 10)
}
\arguments{
\item{p}{vector of polynomial coefficients}

\item{lower}{lower bound to search for mode}

\item{upper}{upper bound to serach for mode}
}
\value{
a list(mu = mu, var=var) for a normal distribution
  representing the Laplace approximation
}
\description{
Make a quadratic approximation to the polynomial log-density at the mode
}
