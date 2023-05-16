lse <- function(x){
  matrixStats::logSumExp(x)
}

sigmoid <- function(x){
  return(1/(1 + exp(-x)))
}

softmax <- function(x){
  return(exp(x - lse(x)))
}

compute_kl_gaussian <- function (mu, var, mu0 = 0, var0 = 1){
  kl <- 0.5 * (log(var0) - log(var) + var/var0 + (mu - mu0)^2/var0 - 1)
  return(kl)
}

compute_kl_categorical <- function (alpha, pi=rep(1/length(alpha), length(alpha))){
  kl <- sum(alpha * (log(alpha) - log(pi)), na.rm = T)
  return(kl)
}

compute_kl_polynomial <- function(q, p){
  f <- function(x){
    qx <- evaluate_poly(q, x)
    px <- evaluate_poly(p, x)
    return(exp(qx) * (qx - px))
  }
  kl <- cubature::hcubature(f, -20, 20)$integral
  return(kl)
}

#' Compute KL divergence from normal (quadratic) to polynomial distributions
#'
#' Compute KL divergence between normal distribution and a distribution with
#' polynomial log-density. Faster that `compute_kl_polynomial` because we
#' can analytically compute moments undert the normal distribution.
#' @param mu mean of normal distribution
#' @param var variance of normal distribution
#' @param p vector of polynomial coefficients
#' @param normalize whether to normalize the polynomial density--
#'   note: can optimize without normalizing, but can return negative values
#' @return the KL divergence KL[N(mu, var) | p]
compute_kl_normal_poly <- function(mu, var, p, normalize=F){
  if(normalize){
    p <- polysusie:::normalize_polynomial_log_density(p)
  }
  K = length(p) - 1
  moments <- c(1, compute_normal_moments2(mu, var, K))
  kl <- - 0.5 * log(2 * pi * exp(1) * var) - sum(moments * p)
  return(kl)
}

#' Make polynomial function
#'
#' Take a vector of coefficients and make a vectorized
#' function that evaluates the polynomial
#' @param p a vector of coefficients in increasing degree e.g. c(c0, c1, ... c_k)
#' @return a function `f` such that f(x)
#'  evaluates the polynomial with coefficients `p` at x
polyf <- function(p){
  f <- function(x){pracma::polyval(rev(p), x)}
  return(f)
}

is_monotone <- function (v) {
  return(all(tail(v, -1) - head(v, -1) >= 0))
}
