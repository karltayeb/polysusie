test_that("gaussian_to_poly and poly_to_gaussian are inverses", {
  # poly -> gauss -> poly
  p <- normalize_polynomial_log_density(c(1, 40, -400))
  gauss <- poly_to_gaussian(p)
  p2 <- with(gauss, gaussian_to_poly(mu, var, 2))
  expect_equal(p, p2)

  # gauss -> poly -> gauss
  gauss <- list(mu = -0.1, var=2.3)
  p <- with(gauss, gaussian_to_poly(mu, var, 2))
  gauss2 <- poly_to_gaussian(p)
  expect_equal(gauss, gauss2)
})
