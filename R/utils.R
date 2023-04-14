lse <- matrixStats::logSumExp

sigmoid <- function(x){
  return(1/(1 + exp(-x)))
}

softmax <- function(x){
  return(exp(x - lse(x)))
}

normal_kl <- function (mu, var, mu0 = 0, var0 = 1){
  kl <- 0.5 * (log(var0) - log(var) + var/var0 + (mu - mu0)^2/var0 - 1)
  return(kl)
}
