logistic_polysusie<- function(X, y, L, prior_variance=1,
                                left=-5, right=5, M=2,
                                max_iter=100, tol=1e-5){
  A <- bernoulli_poly_approx(y, left, right, M)
  fit <- polynomial_approximate_susie(A, X,
                                      prior_variance=prior_variance,
                                      L=L, max_iter = max_iter, tol=tol)
  return(fit)
}
