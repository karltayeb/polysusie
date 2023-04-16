# Impliment polynomial approximation SuSiE

# add functions
add_coef <- function(a, mu){
  return((construct_shift_matrix(mu) %*% a)[,1])
}
add_coef2 <- function(A, moments){
  A2 <- purrr::map(1:nrow(A), ~add_coef(A[.x,], moments[.x,]))
  return(do.call(rbind, A2))
}

# subtract functions
sub_coef <- function(a, mu){
  return(solve(construct_shift_matrix(mu), a))
}
sub_coef2 <- function(A, moments){
  A2 <- purrr::map(1:nrow(A), ~sub_coef(A[.x,], moments[.x,]))
  return(do.call(rbind, A2))
}

polynomial_approximate_susie <- function(A, X, prior_variance, L=5, max_iter=5){
  # 1. Fit SERs, first pass
  sers <- list()
  A2 <- A  # these will have the "residualized" coefficients
  for(l in 1:L){
    sers[[l]] <- polynomial_approximate_ser(A2, X, prior_variance)
    A2 <- add_coef2(A2, sers[[l]]$moments)
  }

  # 2. iterate
  for(i in 1:(max_iter - 1)){
    for(l in 1:L){
      A2 <- sub_coef2(A2, sers[[l]]$moments)
      sers[[l]] <- polynomial_approximate_ser(A2, X, prior_variance)
      A2 <- add_coef2(A2, sers[[l]]$moments)
    }
  }

  # 3. extract posterior
  alpha <- do.call(rbind, purrr::map(1:length(sers), ~ sers[[.x]]$alpha))
  mu <- do.call(rbind, purrr::map(1:length(sers), ~ sers[[.x]]$mu))
  var <- do.call(rbind, purrr::map(1:length(sers), ~ sers[[.x]]$var))

  res <- list(alpha = alpha,
              mu = mu,
              var = var)
  return(res)
}

