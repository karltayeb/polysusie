# compute kl divergence between two normals analytically, and compare to
# our kl computation using quadrature on the log polynomial densities for q and p
test_that("polynomial kl agrees with normal kl", {
  q <- gaussian_to_poly(1, 0.5, 2)
  p <- gaussian_to_poly(0, 1, 2)
  kl1 <- compute_kl_polynomial(q, p)
  kl2 <- compute_kl_gaussian(1, 0.5, 0, 1)
  expect_equal(kl1, kl2)
})
