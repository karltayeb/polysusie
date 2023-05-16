# Functions for computing moments and normalizing polynomial log density

evaluate_poly <- function(p, x){
  pracma::polyval(rev(p), x)
}

plot_poly <- function(p, range=c(-5, 5)){
  f <- function(x) {evaluate_poly(p, x)}
  plot(f, xlim = range)
}

evaluate_exp_poly <- function(p, x){
  exp(evaluate_poly(p, x))
}

make_exp_poly <- function(p){
  f <- function(x){evaluate_exp_poly(p,x)}
  return(f)
}

plot_exp_poly <- function(p, range=c(-5, 5)){
  f <- make_exp_poly(p)
  plot(f, xlim = range)
}


#' Normalize polynomial log density
#'
#' update intercept so that \int \exp \{f(x) \} dx = 1
#' computed numerically via quadrature.
#' where f(x) is a polynomial with coefficients given by argument `p`
#' @param p polynomial coefficients ordered by increasing degree e.g. (p0, ... pM)
normalize_polynomial_log_density <- function(p){
  f <- make_exp_poly(p)
  Z <- cubature::hcubature(f, -20, 20)$integral
  p[1] <- p[1] - log(Z)
  return(p)
}


make_moment_function <- function(p, k){
  f <- function(x){x^k * evaluate_exp_poly(p,x)}
  return(f)
}

#' Compute moment
#'
#' Compute the moments of a distribution with polynomial log-density
#' @param p polynomial coefficients ordered by increasing degree e.g. (p0, ... pM)
#'  Assumes p is normalized!
#' @param k the moment to compute
#' @returns E[X^k] = \int x^k exp(f(x)) dx
compute_moment_polynomial <- function(p, k, range=c(-20,20)){
  f <- make_moment_function(p, k)
  cubature::hcubature(f, -20, 20)
  mu_k <- cubature::hcubature(f, range[1], range[2])$integral
  return(mu_k)
}

compute_moments_polynomial <- function(p, K, range=c(-20, 20)){
  mu <- purrr::map_dbl(1:K, ~compute_moment_polynomial(p, .x, range))
  return(c(1, mu))
}

#' Compute normal moments
#'
#' Compute first `k` moments of a normal distribution (including 0th moment)
#' @param mu mean
#' @param var variance
#' @param k number of moments to copmute
compute_normal_moments <- function(mu, var, k){
  return(purrr::map_dbl(0:k, ~actuar::mnorm(.x, mu, sqrt(var))))
}

#' Compute normal moments
#'
#' Recursively compute moments of normal distribution using reccurence relation
#' M_k = mu M_{k-1} + var * (k-1) * M_{k-2}
#' @param mu mean parameter of normal distribution
#' @param var variance of normal distribution
#' @param K number of moments to compute
#' @return a length K vector with the first K (un-centered) moments of N(mu, var)
compute_normal_moments2 <- function(mu, var, K){
  moments <- rep(0, K)
  moments[1] <- mu
  moments[2] <- mu^2 + var
  for (k in 3:K){
    moments[k] <- mu * moments[k-1] + (k-1) * var * moments[k-2]
  }
  return(moments)
}

#' compute k moments for psi = xb, b ~ N(mu, var)
compute_psi_moments <- function(x, b_moments){
  M <- length(b_moments) - 1
  psi_moments <- do.call(cbind, purrr::map(0:M, ~ (x^.x) * b_moments[.x + 1]))
  return(psi_moments)
}
