# functions for approximating specific likelihoods e.g. logistic + sigmoid link

#' Make approximation
#'
#' make polynomial approximation on the interval [left, right]
#' @param f a 1d function
#' @param left boundary of approximation
#' @param right boundary of approximation
#' @param M degree of the approximation
#' @param plot boolean to create a plot of the function
#' @export
make_approximation <- function(f, left=-4, right=4, M=2, plot=F){
  p <- rev(pracma::polyApprox(f, left, right, n=M)$p)
  if(plot){
    x <- seq(left-2, right+2, by=0.1)
    plot(f, left-2, right+2)
    lines(x, polyval2(p, x), col='red', lty='dotted')
    abline(v=left); abline(v=right)
  }
  return(p)
}

loglik0 <- function(psi){
  log(sigmoid(-psi))
}

#' Bernoulli + sigmoid polynomial approximation
#'
#' @param y observation
#' @param left left boundary or approximation
#' @param right boundardy of approximation
#' @param M degree of approximation
#' @returns coefficients for polynomial approximation f(\psi) \approd log p(y | \psi)
bernoulli_poly_approx <- function(y, left, right, M){
  n <- length(y)
  p0 <- make_approximation(loglik0, left, right, M)

  # for y=1 flip the sign of odd coefficients (note: 0 indexing)
  p1 <- p0
  p1[seq(2, length(p0), by=2)] <- p1[seq(2, length(p0), by=2)] * -1

  m <- matrix(nrow = n, ncol = M + 1)
  for(i in 1:length(y)){
    if(y[i] == 0){
      m[i,] <- p0
    } else{
      m[i,] <- p1
    }
  }
  return(m)
}

