# Functions for approximating moments of polynomial exponenital families with importance sampling.
# We can use coefficients of the unnormalized log-density directcly.


#' Normalize log weights
#'
#' Convert unnormalized log probabilities into probabilities
#'
#' @param logw unnormalized log-probabilities
#' @return normalized probabilities that sum to 1
normalizelogweights <- function (logw) {
  c <- max(logw)
  w <- exp(logw - c)
  return(w/sum(w))
}


#' Log mean
#'
#' This function computes log(mean(w)) in a numerically stable way
#' @param logw = log(w), w a vector
#' @return log(mean(w))
logmean <- function (logw) {
  c <- max(logw)
  w <- exp(logw - c)
  return(c + log(mean(w)))
}

#' Compute polynomial moments via importance sampling
#'
#' @param p vector of coefficients
#' @param K moments to computes
#' @param left left integration boundary
#' @param right right integration boundary
#' @param n number of points for IS
#' @export
compute_moments_polynomial_is <- function(p, K, left=-10, right=10, n=64){
  # log density function
  f <- function(x){
    return(pracma::polyval(rev(p), x))
  }
  x <- seq(left, right, length.out=n)
  logw <- f(x)
  w <- normalizelogweights(logw)
  moments <- purrr::map_dbl(1:K, ~sum((x^.x)*w))
  return(moments)
}

