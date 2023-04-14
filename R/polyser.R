# Implement single effect regression with polynomial approximation

# for each covariate:
#   1. "scale" the polynomial f(psi_l | gamma=j) = f(x_j*b_l) = f2(b)
#   2. sum the coefficients across observations + prior
#   3. option: project polynomial approximation to gaussian?
# posterior approximation q:
#   1. compute BFs + normalize to compute alpha = q(gamma)
#   2. compute moments for future updates E[psi_l^k] for k = 1, ..., M

#' Polynomial approximate univariate regression
#'
#' Fit a univariate regression with polynomial approximation to the likelihood
#'
#' @param A an n x M matrix of coefficients--
#'    i-th row gives coefficients for polynomial approximating
#'    log p(y_i | psi) as a function of the linear predictor psi = xb
#' @param x a n vector of covariates
#' @param prior_variance the prior variance of the effect
polynomial_approximate_univariate_regression <- function(A, x, prior_coef=NULL){
  # TODO: how to include intercept?
  M <- ncol(A) - 1
  Ascale <- purrr::map(1:nrow(A), ~ A[.x,] * x[.x]^(0:M))
  post_coef_M <- colSums(do.call(rbind, Ascale)) + prior_coef

  # project to Gaussian/quadratic
  post_gauss <- poly_to_gaussian(post_coef_M)

  # compute moments E[psi^k] k = 0,...,M
  # result is a n x (M+1) matrix
  moments <- compute_psi_moments(x, post_gauss$mu, post_gauss$var, M)

  # compute ELBO
  loglik <- sum(A * moments)
  prior_var <- -prior_coef[3]/0.5
  kl <- with(post_gauss, normal_kl(mu, var, 0, prior_var))  # KL[q || p], analytic if q, p Gaussian
  elbo <- loglik - kl

  # summarize posterior
  post <- list(
    coef = coef,
    q = post_gauss,
    moments = moments,
    elbo = elbo
  )
  return(post)
}

polynomial_approximate_ser <- function(A, X, prior_variance){
  # 1. compute prior_coef from prior_variance
  M <- ncol(A) - 1
  prior_coef = make_prior_p(0, prior_var, M)

  # 2. compute polynomial approximate posterior for each column of X
  post <- purrr::map(1:ncol(X), ~polynomial_approximate_univariate_regression(A, X[, .x], prior_coef=prior_coef))

  # 3. compute alpha, NOTE: assuming uniform prior over effect variables
  elbos <- purrr::map_dbl(1:length(post), ~ post[[.x]]$elbo)
  alpha <- softmax(elbos)
  elbo <- mean(elbos)
  moments <- 0
  for(i in 1:length(post)){
    moments <- moments + post[[i]]$moments * alpha[[i]]
  }

  mu <- purrr::map_dbl(1:length(post), ~ post[[.x]]$q$mu)
  var <- purrr::map_dbl(1:length(post), ~ post[[.x]]$q$var)

  post2 <- list(
    mu = mu,
    var = var,
    alpha = alpha,
    elbo = elbo,
    moments = moments,
    post = post
  )
  return(post2)
}


