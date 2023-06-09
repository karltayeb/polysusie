% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poly_ops.R
\name{construct_shift_matrix}
\alias{construct_shift_matrix}
\title{Make shift matrix}
\usage{
construct_shift_matrix(moments)
}
\arguments{
\item{moments}{moments of y}
}
\description{
Generate matrix that maps coefficients of a polynomial
f(x + y) (represented by coefficients p) to coefficients of
f2(x) = E_{p(y)}[f(x+y)]
}
