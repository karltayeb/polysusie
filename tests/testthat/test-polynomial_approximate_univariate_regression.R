test_that("polynomial q and gaussian q agree for M=2", {
  sim <- logisticsusie::sim_susie()
  A <- bernoulli_poly_approx(sim$y, -4, 4, 2)
  prior_coef <- gaussian_to_poly(0, 1, 2)

  uni1 <- polynomial_approximate_univariate_regression(A, sim$X[,1], prior_coef = prior_coef)
  uni2 <- polynomial_approximate_univariate_regression(A, sim$X[,1], prior_coef = prior_coef, gaussian = F)

  expect_equal(uni1$post_coef, uni2$post_coef)
  expect_equal(uni1$b_moments, uni2$b_moments)
  expect_equal(uni1$kl, uni2$kl)
  expect_equal(uni1$elbo, uni2$elbo)
})

