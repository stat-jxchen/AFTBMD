#' Censoring Rate Control
#'
#' This function finds a number \eqn{\tau} to control the censoring rate.
#'
#' @param c.r a number, the censoring rate you want to control.
#' @param n the number of observations you want to generate.
#' @param d the number of covariates you want to generate.
#' @param rho the number to control the correlation between two adjacent columns of generated data.
#' @param fun_list a list, whose elements are functions to be applied to covariates specified by \code{col_index} one by one.
#' @param col_index a vector, which specifies the indexes of covariates corresponding to \code{fun_list} in the data matrix.
#' @param tol a number, the tolerance error of \code{c.r} when searching for \eqn{\tau}.
#' @param N.max a number, which is used to control the number of observations (\code{N.max} * \code{n}) for searching for \eqn{\tau}.
#' @return a number  \eqn{\tau}  to control the censoring rate \code{c.r}.
#' @details Assume the logarithm of the survival time \eqn{T_i} following an
#' additive accelerated failure time model:
#' \deqn{T_i = \sum_{j=1}^pg_j(Z_{ij})+\varepsilon_i,}
#' where \eqn{Z_{ij}} denotes the value of \eqn{jth} covariate of
#' \eqn{ith} individual. \code{col_index} specifies which covariates have
#' a functional effect on survival time while \code{fun_list} specifies this
#' functions. Assume the censoring time follow the uniform distribution
#' \eqn{U[0,\tau]}, This function is used to find \eqn{\tau} to control
#' the censoring rate by using bisection method.
#' @examples
#' f1 <- function(x) return(x)
#' f2 <- function(x) return(-x)
#' f3 <- function(x) return(x^2)
#' find_tau(c.r = 0.5, n = 200, d = 30, rho = 0.25,
#'          fun_list = list(f1, f2, f3), col_index = 1:3)
#' @export
find_tau <- function(c.r, n, d, rho, fun_list, col_index,
                     tol = 1e-5, N.max = 30){
  NN <- N.max * n
  x <- AR(NN, d, rho)
  wn <- rnorm(NN, 0, 1)
  y <- apply_fun(x, fun_list, col_index) + wn
  calculate.rate <- function(tau){
    # Generate NN random numbers in 0~tau
    c <- log(runif(NN, 0, tau))
    delta <- ifelse(y <= c, 1, 0)
    c.rate <- 1 - mean(delta)
  }
  c.rate <- 1
  low <- 0
  high <- 1e+9
  while(abs(c.rate - c.r) > tol & low < high){
    tau <- (low + high)/2
    c.rate <- calculate.rate(tau)
    if(c.rate < c.r){
      high <- tau
    }else{
      low <- tau
    }
  }
  return(tau)
}
