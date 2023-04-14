# Functions for dealing with quadratic log-density (Gaussian)
# approximation of polynomial log-density

#' convert (unnormalized) polynomial density to gaussian approximation
#' p the coefficients of a polynomial in increasing order p = c(p0, p1, ..., pK)
poly_to_gaussian <- function(p){
  p <- rev(p)
  d <- pracma::polyder(p)
  d2 <- pracma::polyder(d)

  #f <- function(x){polyval2(p, x)}
  #mu <- optimize(f, interval = c(-100, 100), maximum = T)$maximum
  roots <- Re(pracma::polyroots(d)$root)
  mu <- roots[which.max(pracma::polyval(p, roots))]
  var <- - 1 / pracma::polyval(d2, mu)
  return(list(mu=mu, var=var))
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

#' compute k moments for psi = xb, b ~ N(mu, var)
compute_psi_moments <- function(x, mu, var, k){
  normal_moments <- compute_normal_moments(mu, var, k)
  psi_moments <- do.call(cbind, purrr::map(0:k, ~ (x^.x) * normal_moments[.x + 1]))
}


#' Make shift matrix
#'
#' Generate matrix that maps coefficients of a polynomial
#' f(x + y) (represented by coefficients p) to coefficients of
#' f2(x) = E_{p(y)}[f(x+y)]
#' @param moments moments of y
construct_shift_matrix <- function(moments){
  # construct map
  K <- length(moments) - 1
  M <- matrix(nrow= K+1, ncol=K+1)
  for(j in 0:K){
    for(k in 0:K){
      M[j+1, k+1] <- choose(k, j) * moments[[max(k-j+1, 1)]]
    }
  }
  return(M)
}


a <- rnorm(10)
moments <- rnorm(10)
sgn <- c(1, -1)

M1 <- make_shift_matrix(moments)
M2 <- make_shift_matrix(moments*sgn)

