test_polysusie <- function(){
  sim <- logisticsusie::sim_susie()
  fit <- with(sim, logistic_polysusie(X, y, L=3, M=10))

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

  tictoc::tic()
  A <- bernoulli_poly_approx(sim$y, 4, 14)
  susie <- polynomial_approximate_susie(A, sim$X, prior_variance=1, L=5)
  cs <- logisticsusie:::get_all_cs(susie$alpha)
  tictoc::toc()


  susie2 <- with(sim, logisticsusie::binsusie(X, y, L=5, estimate_prior_variance=F))
  plot(susie2$sers[[1]]$mu, susie$mu[1,]); abline(0, 1, col='red')
}

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
