compare_with_jj <- function(){
  sim <- logisticsusie::sim_ser()
  A <- bernoulli_poly_approx(sim$y, 4, 2)
  ser <- polynomial_approximate_ser(A, sim$X, 1)
  ser2 <- with(sim, logisticsusie::binser(X, y))

  plot(log(ser2$alpha), log(ser$alpha)); abline(0, 1, col='red')
  plot(ser2$mu, ser$mu); abline(0, 1, col='red')
}


test_polysusie <- function(){
  sim <- logisticsusie::sim_susie()

  tictoc::tic()
  A <- bernoulli_poly_approx(sim$y, 4, 2)
  susie <- polynomial_approximate_susie(A, sim$X, prior_variance=1, L=5)
  cs <- logisticsusie:::get_all_cs(susie$alpha)
  tictoc::toc()

  tictoc::tic()
  A <- bernoulli_poly_approx(sim$y, 4, 6)
  susie <- polynomial_approximate_susie(A, sim$X, prior_variance=1, L=5)
  cs <- logisticsusie:::get_all_cs(susie$alpha)
  tictoc::toc()

  tictoc::tic()
  A <- bernoulli_poly_approx(sim$y, 4, 10)
  susie <- polynomial_approximate_susie(A, sim$X, prior_variance=1, L=5)
  cs <- logisticsusie:::get_all_cs(susie$alpha)
  tictoc::toc()

  susie2 <- with(sim, logisticsusie::binsusie(X, y, L=5, estimate_prior_variance=F))
  plot(susie2$sers[[1]]$mu, susie$mu[1,]); abline(0, 1, col='red')

  ser2 <- with(sim, logisticsusie::binser(X, y))

  plot(log(ser2$alpha), log(ser$alpha)); abline(0, 1, col='red')
  plot(ser2$mu, ser$mu); abline(0, 1, col='red')
}
