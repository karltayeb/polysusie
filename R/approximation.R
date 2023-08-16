# functions for approximating specific likelihoods e.g. logistic + sigmoid link

#' Make approximation
#'
#' make polynomial approximation on the interval [left, right]
#' @param f a 1d function
#' @param left boundary of approximation
#' @param right boundary of approximation
#' @param M degree of the approximation
#' @param plot boolean to create a plot of the function
#' @export
make_approximation <- function(f, left=-4, right=4, M=2, plot=F){
  p <- rev(pracma::polyApprox(f, left, right, n=M)$p)
  if(plot){
    x <- seq(left-2, right+2, by=0.1)
    plot(f, left-2, right+2)
    lines(x, polyval2(p, x), col='red', lty='dotted')
    abline(v=left); abline(v=right)
  }
  return(p)
}

loglik0 <- function(psi){
  log(sigmoid(-psi))
}

loglik1 <- function(psi){
  psi + log(sigmoid(-psi))
}

#' Bernoulli + sigmoid polynomial approximation
#'
#' @param y observation
#' @param left left boundary or approximation
#' @param right boundardy of approximation
#' @param center midpoint of approximation region
#' @param M degree of approximation
#' @returns coefficients for polynomial approximation f(\psi) \approd log p(y | \psi)
bernoulli_poly_approx <- function(y, left, right, M){
  n <- length(y)

  # make polynomial for each unique combination of (y, left, right)

  # put data into a tibble
  dat <- tibble::tibble(y = y) %>%
    dplyr::mutate(left = left, right = right)

  # compute unique coefficients for y==0
  coef0 <- dat %>%
    unique() %>%
    dplyr::filter(y == 0) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coef = list(make_approximation(loglik1, left, right, M)))

  # compute unique coefficients for y==1
  coef1 <- dat %>%
    unique() %>%
    dplyr::filter(y == 1) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coef = list(make_approximation(loglik1, left, right, M)))

  # concatenate results
  coef <- rbind(coef0, coef1)

  # map coefficents to observations, returns a n x (M+1) matrix
  m <- dat %>%
    dplyr::inner_join(coef) %>%
    {do.call(rbind, .$coef)}
  return(m)
}


#' Bernoulli + sigmoid polynomial approximation 2
#'
#' @param y observation
#' @param left left boundary or approximation
#' @param right boundardy of approximation
#' @param center midpoint of approximation region
#' @param M degree of approximation
#' @returns coefficients for polynomial approximation f(\psi) \approd log p(y | \psi)
bernoulli_poly_approx <- function(y, left, right, M){
  n <- length(y)

  # make polynomial for each unique combination of (y, left, right)
  intervals <- tibble::tibble(left = left, right = right) %>%
    dplyr::distinct() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      coef0 = list(make_approximation(loglik0, left, right, M)),
      coef1 = list(make_approximation(loglik1, left, right, M))
    ) %>%
    dplyr::ungroup()

  m <- dplyr::tibble(y=y, left=left, right=right) %>%
    dplyr::left_join(intervals) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coef = list((1 - y) * coef0 + y * coef1)) %>%
    {do.call(rbind, .$coef)}

  return(m)
}
