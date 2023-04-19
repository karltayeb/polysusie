# functions for approximating specific likelihoods e.g. logistic + sigmoid link

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

