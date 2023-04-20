# Impliment polynomial approximation SuSiE

polynomial_susie_compute_elbo <- function(A, sers, intercept){
  loglik <- sum(A[,1])
  kl <- sum(purrr::map_dbl(sers, ~.x$kl)) + intercept$kl
  return(loglik - kl)
}

#' @export
polynomial_approximate_susie <- function(A, X,
                                         prior_variance, L=5, fit_intercept=T,
                                         max_iter=100, tol=1e-8){
  # 0. Setup
  M <- ncol(A)-1
  shifter <- construct_shift_matrix_generator(M) # make shift function
  ones <- rep(1, nrow(A))
  intercept_coef <- gaussian_to_poly(0, 1e8, M) # for intercept model
  A2 <- A  # these will have the "residualized" coefficients

  # 1. Fit SERs, first pass
  if(fit_intercept){
    intercept <- polynomial_approximate_univariate_regression(A2, ones, intercept_coef)
    A2 <- add_intercept(A2, intercept$b_moments, shifter)
  } else{
    intercept <- list(q=list(mu=0, var=0), kl=0)
  }

  sers <- list()
  for(l in 1:L){
    sers[[l]] <- polynomial_approximate_ser(A2, X, prior_variance)
    A2 <- add_coef2(A2, sers[[l]]$moments, shifter = shifter)
  }

  # 2. iterate
  elbos <- -Inf # polynomial_susie_compute_elbo(A2, sers)
  i = 1
  while(i < max_iter){
    for(l in 1:L){
      # update SER
      A2 <- sub_coef2(A2, sers[[l]]$moments, shifter = shifter)
      sers[[l]] <- polynomial_approximate_ser(A2, X, prior_variance)
      A2 <- add_coef2(A2, sers[[l]]$moments, shifter = shifter)

      # update intercept
      if(fit_intercept){
        A2 <- sub_intercept(A2, intercept$b_moments, shifter)
        intercept <- polynomial_approximate_univariate_regression(A2, ones, intercept_coef)
        A2 <- add_intercept(A2, intercept$b_moments, shifter)
      }
    }
    # track ELBO
    elbo <- polynomial_susie_compute_elbo(A2, sers, intercept)
    elbos <- c(elbos, elbo)
    if(diff(tail(elbos, 2)) < tol){
      message('converged')
      break
    }
    i <- i + 1
  }

  # 3. extract posterior
  alpha <- do.call(rbind, purrr::map(1:length(sers), ~ sers[[.x]]$alpha))
  mu <- do.call(rbind, purrr::map(1:length(sers), ~ sers[[.x]]$mu))
  var <- do.call(rbind, purrr::map(1:length(sers), ~ sers[[.x]]$var))

  # summarize
  cs <- get_all_cs(alpha)

  res <- list(alpha = alpha,
              mu = mu,
              var = var,
              intercept = intercept$q,
              elbos = elbos,
              sers=sers,
              A2 = A2,
              converged = (diff(tail(elbos, 2)) < tol),
              cs=cs)
  return(res)
}

