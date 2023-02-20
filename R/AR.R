#' Auto-Regression Time Series
#'
#' This function generates data for auto-regression time series.
#'
#' @param n the number of rows of the result matrix, denotes the number of obversations.
#' @param d the number of columns of the result matrix, denotes the number of covariates.
#' @param rho the number to control the correlation between two adjacent columns of generated data.
#' @return A matrix with \code{n} rows and \code{d} columns.
#' @details The result matrix is a design matrix and the rows
#' denote observations and the columns denote covariates. The covariates
#' generated from a \code{d}-dimension normal distribution with mean zero
#' and covariance matrix \eqn{\Sigma = (\delta_{ij})},
#' where \eqn{\delta_{ij} = \rho^{i-j}}. But since we controled
#' the range between (-1, 1), result obtained by using \code{\link[stats]{cov}}
#' is not approximately \eqn{\rho}.
#' @examples
#' AR(200,30,0.25)
#'
#' @export
AR <- function(n, d, rho){
  wn.sd <- sqrt(1 - rho^2)
  #  wn means white noise
  wn.mat <- matrix(rnorm(n * d), nrow = n, ncol = d) * wn.sd
  Xmat <- matrix(0, nrow = n, ncol = d)
  Xmat[,1] <- rnorm(n)
  for(j in 2:d)
    Xmat[, j] <- rho * Xmat[, j-1] + wn.mat[, j]
  # control the range between (-1, 1)
  Xmat <- matrix(pmax(-1, pmin(1, Xmat)), nrow = n, ncol = d)
  return(Xmat)
}
