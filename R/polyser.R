# Implement single effect regression with polynomial approximation

compute_loglik_coef <- function(A, x){
  loglik_coef <- 0
  n = nrow(A)
  M <- ncol(A) - 1
  loglik_coef <- colSums(A * outer(x, seq(0, M), `^`))
  return(loglik_coef)
}

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
#' @export
polynomial_approximate_univariate_regression <- function(A, x, prior_coef=NULL, gaussian=T){
  # TODO: how to include intercept?
  M <- ncol(A)-1
  loglik_coef <- compute_loglik_coef(A, x)

  if(gaussian){
    # project to Gaussian
    post_gauss <- poly_to_gaussian(loglik_coef + prior_coef)
    post_coef <- gaussian_to_poly(post_gauss$mu, post_gauss$var, M)
    prior_var <- 1 / (-2 * prior_coef[3])
    kl <- with(post_gauss, compute_kl_gaussian(mu, var, 0, prior_var))  # KL[q || p], analytic if q, p Gaussian
    b_moments <- compute_normal_moments(post_gauss$mu, post_gauss$var, M)
  } else{
    # normalize log density. compute moments + kl
    post_coef <- normalize_polynomial_log_density(loglik_coef + prior_coef)
    kl <- compute_kl_polynomial(post_coef, prior_coef)
    b_moments <- compute_moments_polynomial(post_coef, M)
  }

  loglik <- sum(loglik_coef * b_moments)
  elbo <- loglik - kl

  post <- list(
    q = post_gauss,
    post_coef = post_coef,
    b_moments = b_moments,
    elbo = elbo,
    kl = kl
  )
  return(post)
}


compute_kl_polynomial_approximate_ser <- function(ser){
  mu <- 0
  var <- 1
}

#' @export
polynomial_approximate_ser <- function(A, X, prior_variance){
  # 1. compute prior_coef from prior_variance
  M <- ncol(A) - 1
  prior_coef = gaussian_to_poly(0, prior_variance, M)

  # 2. compute polynomial approximate posterior for each column of X
  post <- purrr::map(1:ncol(X), ~polynomial_approximate_univariate_regression(A, X[, .x], prior_coef=prior_coef))

  # 3. compute alpha, NOTE: assuming uniform prior over effect variables
  elbos <- purrr::map_dbl(1:length(post), ~ post[[.x]]$elbo)
  alpha <- softmax(elbos)
  elbo <- mean(elbos)

  moments <- 0
  for(i in 1:length(post)){
    moments <- moments + compute_psi_moments(X[, i], post[[i]]$b_moments) * alpha[[i]]
  }

  mu <- purrr::map_dbl(1:length(post), ~ post[[.x]]$q$mu)
  var <- purrr::map_dbl(1:length(post), ~ post[[.x]]$q$var)
  kl <- sum(alpha * purrr::map_dbl(1:length(post), ~ post[[.x]]$kl)) + compute_kl_categorical(alpha)

  post2 <- list(
    mu = mu,
    var = var,
    alpha = alpha,
    elbo = elbo,
    moments = moments,
    post = post,
    kl = kl
  )
  return(post2)
}


