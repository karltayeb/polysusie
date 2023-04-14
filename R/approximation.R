make_approximation <- function(f, R, n, plot=F){
  p <- rev(pracma::polyApprox(f, -R, R, n =n)$p)
  if(plot){
    S <- R + 2
    x <- seq(-S, S, by=0.1)
    plot(f, -S, S)
    lines(x, polyval2(p, x), col='red', lty='dotted')
    abline(v=-R); abline(v=R)
  }
  return(p)
}

loglik0 <- function(psi){
  log(sigmoid(-psi))
}

#' get approximate polynomial representation of the data y
bernoulli_poly_approx <- function(y, R, k){
  n <- length(y)
  p0 <- make_approximation(loglik0, R, k)

  # for y=1 flip the sign of odd coefficients (note: 0 indexing)
  p1 <- p0
  p1[seq(2, length(p0), by=2)] <- p1[seq(2, length(p0), by=2)] * -1

  m <- matrix(nrow = n, ncol = k + 1)
  for(i in 1:length(y)){
    if(y[i] == 0){
      m[i,] <- p0
    } else{
      m[i,] <- p1
    }
  }
  return(m)
}

#' f(x + c) = f2(x)
shift_polynomial <- function(p, c){
  # construct map
  K <- length(p) - 1
  M <- matrix(nrow= K+1, ncol=K+1)
  for(j in 0:K){
    for(k in 0:K){
      M[j+1, k+1] <- choose(k, j) * c**(k - j)
    }
  }

  coef_new <- (M %*% p)[, 1]
  return(coef_new)
}

# change back to original scale
# f(bx) = f2(x)
scale_polynomial <- function(p, b){
  K <- length(p) - 1
  coef_new <- p * sapply(0:K, function(k) b**k)
  return(coef_new)
}

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


#' Make shift matrix
#'
#' Generate matrix that maps coefficients of a polynomial
#' f(x + y) (represented by coefficients p) to coefficients of
#' f2(x) = E_{p(y)}[f(x+y)]
#' @param moments moments of y
make_shift_matrix <- function(moments){
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

#' Transform coefficients of a polynomial f(x + y) (represented by coefficients p)
#' to coefficients of f2(x) = E_{p(y)}[f(x+y)]
#' @param p K+1 coefficients of a degree-k polynomial
#' @param moments moments of y (including E[y^0] = 1)
shift_polynomial3 <- function(p, moments){
  M <- make_shift_matrix(moments)
  p_new <- (M %*% p)[, 1]
  return(p_new)
}

compute_normal_moments <- function(mu, var, k){
  return(purrr::map_dbl(0:k, ~actuar::mnorm(.x, mu, sqrt(var))))
}


#' update q with polynomial approximation of arbitrary degree
polynomial_update3 <- function(m, X, prior_p, q){
  K <- ncol(m) - 1
  p <- ncol(X)
  n <- nrow(X)
  for(j in 1:p){
    m_tilde <- m
    for(k in (1:p)[-j]){
      moments <- compute_psi_moments(X[, k], q[[k]]$mu, q[[k]]$var, K)
      m_tilde <- do.call(rbind, lapply(1:n, function(i) shift_polynomial3(
        m_tilde[i,], moments[i,])))
    }

    # scale-- new polynomial in terms of b_j
    m_hat <- do.call(rbind, lapply(1:n, function(i) scale_polynomial(
      m_tilde[i,], X[i, j])))

    # compute posterior polynomial
    m_post <- colSums(m_hat) + prior_p[[j]]

    # find gaussian approximation
    q[[j]] <- poly_to_gaussian(m_post)
    q[[j]]$m_post <- m_post
  }

  return(q)
}

#' Represent normal log-density as a polynomial
make_prior_p <- function(mu, var, M){
  prior_p <- c(
    -0.5 * (mu^2/var +log(2 * pi * var)),
    mu / var,
    -0.5 / var
  )

  prior_p <- c(prior_p, rep(0, M-3+1))
  return(prior_p)
}

logistic_polynomial_approximation <- function(y, X, R, K=2){
  # observations in polynomial coeeficients
  m <- bernoulli_poly_approx(y, R, K)
  q <- list()
  prior_p <- list()
  for(j in 1:p){
    prior_p[[j]] <- c(c(0, 0, -0.5), rep(0, K-2)) # extend polynomial to agree with m
    q[[j]] <- list(mu = 0, var=1) # initialize normal posterior
  }

  # iteratively update
  param_history <- list()
  param_history[[1]] <- q
  for(i in 1:50){
    q <- polynomial_update3(m, X, prior_p, q)
    param_history[[i+1]] <- q
  }
  return(param_history)
}
