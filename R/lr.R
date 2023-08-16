#' @export
logistic_polysusie <- function(X, y, L, prior_variance=1,
                                left=5, right=5, center=0, M=2,
                                max_iter=100, tol=1e-5, gaussian=T){
  if(!gaussian){
    warning('Non-Gaussian variational approximation not yet implimented')
  }

  # center approximations at ridge regression predictions
  if(center == 'ridge'){
    ridge <- glmnet::cv.glmnet(X, y, family='binomial')
    center <- drop(predict(ridge, X, s='lambda.min'))
  }
  A <- bernoulli_poly_approx(y, center - left, center + right, M)
  fit <- polynomial_approximate_susie(A, X,
                                      prior_variance=prior_variance,
                                      L=L, max_iter = max_iter, tol=tol)
  return(fit)
}
