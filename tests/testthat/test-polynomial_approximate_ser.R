compare_with_jj <- function(){
  sim <- logisticsusie::sim_ser()
  A <- bernoulli_poly_approx(sim$y, 4, 2)
  prior_coef = gaussian_to_poly(0, 1, 2)
  uni <- polynomial_approximate_univariate_regression(A, sim$X[,1], prior_coef)

  ser <- polynomial_approximate_ser(A, sim$X, 1)
  ser2 <- with(sim, logisticsusie::binser(X, y, estimate_prior_variance = F))

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

  tictoc::tic()
  A <- bernoulli_poly_approx(sim$y, 4, 14)
  susie <- polynomial_approximate_susie(A, sim$X, prior_variance=1, L=5)
  cs <- logisticsusie:::get_all_cs(susie$alpha)
  tictoc::toc()


  susie2 <- with(sim, logisticsusie::binsusie(X, y, L=5, estimate_prior_variance=F))
  plot(susie2$sers[[1]]$mu, susie$mu[1,]); abline(0, 1, col='red')
}


compare_solvers <- function() {
  mu <- rnorm(5) + 1
  a <- rnorm(5) + 1
  s1 <- backsolve(construct_shift_matrix(mu), a)
  s2 <- solve(construct_shift_matrix(mu), a)
  testthat::expect_equal(s1, s2)
}


check_elbo_monotone <- function(){
  sim <- logisticsusie::sim_susie()

  tictoc::tic()
  A <- bernoulli_poly_approx(sim$y, 4, 2)
  susie <- polynomial_approximate_susie(A, sim$X, prior_variance=1, L=5, max_iter = 5, fit_intercept = T)
  plot(susie$elbos)
  cs <- logisticsusie:::get_all_cs(susie$alpha)
  tictoc::toc()


  susie <- polynomial_approximate_susie(A, sim$X, prior_variance=1, L=5, max_iter = 1, fit_intercept = F)

  A2_init <- susie$A2
  sers_init <- susie$sers
  prior_variance <- 1
  X <- sim$X

  A2_history <- list()
  A2 <- A2_init
  sers <- sers_init
  elbos <- rep(1, 5)
  lls <- rep(1, 5)
  kls <- rep(1, 5)
  for(l in 1:5){
    A2 <- sub_coef2(A2, sers[[l]]$moments, shifter = shifter)
    sers[[l]] <- polynomial_approximate_ser(A2, X, prior_variance)
    A2 <- add_coef2(A2, sers[[l]]$moments, shifter = shifter)
    A2_history[[l]] <- A2
    lls[l] <- sum(A2[, 1])
    kls[l] <- sum(purrr::map_dbl(sers, ~.x$kl))
    elbos[l] <- polynomial_susie_compute_elbo(A2, sers, list(kl = 0))
  }


  # testthat::is_true(logisticsusie:::.monotone(susie$elbos))
}


profile_polysusie <- function(){
  sim <- logisticsusie::sim_susie(p=50)
  A <- bernoulli_poly_approx(sim$y, 4, 2)
  susie <- polynomial_approximate_susie(A, sim$X, prior_variance=1, L=5)

  profvis::profvis({
    susie <- polynomial_approximate_susie(A, sim$X, prior_variance=1, L=5, max_iter = 20)
  })
}

