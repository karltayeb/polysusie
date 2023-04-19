# Operations on polynomial coefficients

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

#' Construct coef matrix
#'
#' Upper triangular matrix with entires M[j+1, k+1] = (k choose j) k >= 0
construct_coef_matrix <- function(K){
  # construct map
  M <- matrix(nrow= K+1, ncol=K+1)
  for(j in 0:K){
    for(k in 0:K){
      M[j+1, k+1] <- choose(k, j)
    }
  }
  return(M)
}

#' Construct shift matrix generator
#'
#' Returns a function for constructing the shift matrix for a polynomial
#' of degree M
#' we construct this function to avoid repeated callse to `construct_coef_matrix`
#' with a fixed coef matrix we can construct the shift matrix by
#'    1. constructing a toeplitz matrix of the moments "moment matrix"
#'    2. an elementwise multiplication of the moment matrix and coef matrix
construct_shift_matrix_generator <- function(M){
  coef_mat <- construct_coef_matrix(M)
  f <- function(moments){
    return(stats::toeplitz(moments) * coef_mat)
  }
  return(f)
}

powers <- function(x, M){
  #cumprod(c(1, rep(x, M)))
  x^(0:M)
}

mpowers <- memoise::memoise(powers)

scale_coefficients <- function(a, x, M){
  return(a * powers(x, M))
}

scale_coefficients2 <- function(A, x){
  M <- ncol(A) - 1
  purrr::map(1:nrow(A), ~ scale_coefficients(A[.x,], x[.x], M))
}

# add functions
add_coef <- function(a, mu, shifter = construct_shift_matrix){
  return((shifter(mu) %*% a)[,1])
}
add_coef2 <- function(A, moments, shifter = construct_shift_matrix){
  A2 <- purrr::map(1:nrow(A), ~add_coef(A[.x,], moments[.x,], shifter))
  return(do.call(rbind, A2))
}

add_intercept <- function(A, moments, shifter=construct_shift_matrix){
  return(t(shifter(moments) %*% t(A)))
}

# subtract functions
sub_coef <- function(a, mu, shifter = construct_shift_matrix){
  return(backsolve(shifter(mu), a))
}
sub_coef2 <- function(A, moments, shifter = construct_shift_matrix){
  A2 <- purrr::map(1:nrow(A), ~sub_coef(A[.x,], moments[.x,], shifter))
  return(do.call(rbind, A2))
}

sub_intercept <- function(A, moments, shifter=construct_shift_matrix){
  shift <- shifter(moments)
  return(t(backsolve(shift, t(A))))
}


