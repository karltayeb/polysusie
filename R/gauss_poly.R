# functions for finding gaussian approximation to log-polynomial density

#' Polynomial to Gaussian via Laplace approximation
#'
#' convert (unnormalized) polynomial density to gaussian approximation
#' by taking the mode of the polynomial density and taking a laplace approximation
#' NOTE: this isn't necessarily the gaussian that minimizes KL(q || p)?
#' @param p the coefficients of a polynomial in increasing order p = c(p0, p1, ..., pK)
#' @param poly boolean if true return polynomial representation of gaussian approximation
#' @returns a gaussian centered on the mode of exp(p)
poly_to_gaussian_mode <- function(p, poly = F){
  p <- rev(p)
  d <- pracma::polyder(p)
  d2 <- pracma::polyder(d)

  roots <- Re(pracma::polyroots(d)$root)
  mu <- roots[which.max(pracma::polyval(p, roots))]
  var <- - 1 / pracma::polyval(d2, mu)
  if(poly){
    res <- gaussian_to_poly(mu, var, 2)
  } else{
    res <- list(mu=mu, var=var)
  }
  return(res)
}

#' Convert polynomial to Gaussian
#'
#' @param p the coefficients of a polynomial in increasing order p = c(p0, p1, ..., pK)
#' @returns a gaussian approximation to p
poly_to_gaussian <- function(p, poly=F){
  M <- length(p) - 1

  if(M > 2){
    g <- poly_to_gaussian_mode(p)
  } else{
    # read the mean and variance off of coefficients
    var = 1/ (-2 * p[3])
    mu = p[2] * var
    g <- list(mu=mu, var = var)
  }
  return(g)
}

#' Represent normal log-density as a polynomial
gaussian_to_poly <- function(mu, var, M=2){
  prior_p <- c(
    -0.5 * (mu^2/var +log(2 * pi * var)),
    mu / var,
    -0.5 / var
  )
  prior_p <- c(prior_p, rep(0, M-3+1))
  return(prior_p)
}
