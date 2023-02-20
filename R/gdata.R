#' Generate Simulation Data
#'
#' This function generates simulation survival data.
#'
#' @param tau a number,which is used to generate unif(0, tau) random number.
#' @param n the number of observations you want to generate.
#' @param d the number of covariates you want to generate.
#' @param rho the number to control the covariance between two adjacent columns of generated data.
#' @param fun_list a list, whose elements are functions to be applied to covariates specified by \code{col_index} one by one.
#' @param col_index a vector, which specifies the indexes of covariates corresponding to \code{fun_list} in the data matrix.
#' @return A list containing:\tabular{ll}{
#'    \code{y} \tab A numeric vector recording logarithm survival time. \cr
#'    \tab \cr
#'    \code{x} \tab The covariates matrix. \cr
#'    \tab \cr
#'    \code{delta} \tab A numeric vector that consists of status indicator, normally 0=alive, 1=dead. \cr
#' }
#' @export
#' @examples
#' f1 <- function(x) return(x)
#' f2 <- function(x) return(-x)
#' f3 <- function(x) return(x^2)
#' gdata(tau = 1, n = 200, d = 30, rho = 0.25,
#'       fun_list = list(f1, f2, f3), col_index = 1:3)
gdata <- function(tau, n, d, rho, fun_list, col_index){
  x <- AR(n, d, rho)
  wn <- rnorm(n, 0, 1)
  y1 <- apply_fun(x, fun_list, col_index) + wn
  c <- log(runif(n, 0, tau))
  y <- pmin(y1, c)
  delta <- as.numeric(y1 <= c)
  # y1 is actually unknown, only y and delta can be seen in actual observation
  return(list(y = y, x = x, delta = delta))
}
