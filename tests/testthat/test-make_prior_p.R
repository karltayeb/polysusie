test_that("gaussian_to_poly is normalized", {
  prior_coef <- gaussian_to_poly(0, 1, 2)
  prior_coef_norm <- normalize_polynomial_log_density(prior_coef)
  expect_equal(prior_coef, prior_coef_norm)

  prior_coef <- gaussian_to_poly(1, 0.1, 2)
  prior_coef_norm <- normalize_polynomial_log_density(prior_coef)
  expect_equal(prior_coef, prior_coef_norm)
})
