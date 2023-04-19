lse <- matrixStats::logSumExp

sigmoid <- function(x){
  return(1/(1 + exp(-x)))
}

softmax <- function(x){
  return(exp(x - lse(x)))
}

compute_kl_gaussian <- function (mu, var, mu0 = 0, var0 = 1){
  kl <- 0.5 * (log(var0) - log(var) + var/var0 + (mu - mu0)^2/var0 - 1)
  return(kl)
}

compute_kl_categorical <- function (alpha, pi=rep(1/length(alpha), length(alpha))){
  kl <- sum(alpha * (log(alpha) - log(pi)), na.rm = T)
  return(kl)
}

compute_kl_polynomial <- function(q, p){
  f <- function(x){
    qx <- evaluate_poly(q, x)
    px <- evaluate_poly(p, x)
    return(exp(qx) * (qx - px))
  }
  kl <- cubature::hcubature(f, -20, 20)$integral
  return(kl)
}


is_monotone <- function (v) {
  return(all(tail(v, -1) - head(v, -1) >= 0))
}
