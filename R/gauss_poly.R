# functions for finding gaussian approximation to log-polynomial density

#' Laplace approximation for polynomial density
#'
#' Make a quadratic approximation to the polynomial log-density at the mode.
#' This version is a bit slower for high degree polynomials,
#' but it doesnt require specifying a range to optimize on.
#' @param p vector of polynomial coefficients
#' @param lower lower bound to search for mode
#' @param upper upper bound to serach for mode
#' @return a list(mu = mu, var=var) for a normal distribution
#'   representing the Laplace approximation
poly_to_gaussian_laplace2 <- function(p, poly = F){
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


#' Laplace approximation for polynomial density
#'
#' Make a quadratic approximation to the polynomial log-density at the mode
#' @param p vector of polynomial coefficients
#' @param lower lower bound to search for mode
#' @param upper upper bound to serach for mode
#' @return a list(mu = mu, var=var) for a normal distribution
#'   representing the Laplace approximation
poly_to_gaussian_laplace <- function(p, lower=-10, upper=10){
  # find mode
  f <- polyf(p)
  mu <- optimise(f, lower = lower, upper= upper, maximum = T)$maximum

  # compute variance at mode-- Laplace approximation
  var <- p %>% rev() %>%
    pracma::polyder() %>%
    pracma::polyder() %>%
    pracma::polyval(mu) %>%
    {-1/.}

  res <- list(mu=mu, var=var)
  return(res)
}

#' Find Gaussian q that minimizes KL[q || p]
#'
#' Minimize KL divergence for fixed-form variantional inference
#' Use a Laplace approximation as initialization, and then
#' perform coordinate ascent on mean and variance
#' @param p vector of polynomial coefficients
#' @param tol tolerance on change in KL divergence between iterations
#' @param max_iter maximum numbr of iterations
#' @return a list(mu=mu, var=var) with parameters of normal distribution
minimize_kl_normal_poly <- function(p, tol=1e-5, max_iter=10){
  q_init <- poly_to_gaussian_laplace2(p)
  mu <- q_init$mu
  var <- q_init$var
  kls <- compute_kl_normal_poly(mu, var, p)

  # set optimization range based on Laplace approximation
  mu_l <- mu - 10 * var
  mu_u <- mu + 10 * var

  var_l <- var/100
  var_u <- var*100

  # iterate
  for(i in 2:max_iter){
    # optimize mu on approximation interval
    f_mu <- function(mu){purrr::map_dbl(mu, ~compute_kl_normal_poly(.x, var, p))}
    mu <- optimise(f_mu, lower=mu_l, upper=mu_u, maximum = F)$minimum

    # optimize var
    f_var <- function(var){purrr::map_dbl(var, ~compute_kl_normal_poly(mu, .x, p))}
    var <- optimise(f_var, lower=var_l, upper=var_u, maximum = F)$minimum

    kls <- c(compute_kl_normal_poly(mu, var, p), kls)

    converged <- abs(diff(head(kls, 2))) < tol
    if(converged){
      break
    }
  }
  return(list(mu=mu, var=var, kl=kls[1], iter=i, converged=converged))
}

#' Convert polynomial to Gaussian
#'
#' @param p the coefficients of a polynomial in increasing order p = c(p0, p1, ..., pK)
#' @returns a gaussian approximation to p
poly_to_gaussian <- function(p, poly=F){
  M <- length(p) - 1

  if(M > 2){
    g <- minimize_kl_normal_poly(p)
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

