#' Select the Optimal Estimator
#'
#' This function selects the optimal estimator under the EBIC criterion.
#'
#' @param theta.mat a three-dimensional array of estimated values of coefficients under each pair of tuning parameters.
#' The structure refers to \code{Details}
#' @param n the number of observations.
#' @param Loss a matrix containing the value of loss function under each pair of tuning parameters.The structure refers to \code{Details}
#' @return The optimal estimator under the EBIC.
#' @details \code{theta[i,j,]} contains the estimator under the \eqn{ith} value of the first tuning parameter and
#' the \eqn{jth} value of the second parameter.\code{Loss[i,j]} is the weighted least squares loss function value when the estimator is \code{theta[i,j,]}.
#' The EBIC is extended Bayesian information criterion proposed by Chen and Chen
#' (2008) combining the ideas suggested by Gray (1992) to calculate the degree of freedom.{TODO}
#' @note See the example of \code{\link{BMD}} function for an example.
#' @export

selection <- function(theta.mat, n, Loss){
  type <- theta.mat != 0
  pp <- dim(theta.mat)[3]
  ddf <- apply(type, c(1,2), sum)
  EBIC <- 2 * Loss+log(n)*ddf/n+2*log(choose(pp, ddf))/n
  EBIC.min <- min(EBIC)
  idj <- t(which(EBIC == EBIC.min, arr.ind = TRUE))
  result <- list(theta_est = theta.mat[idj[1], idj[2], ], id = idj)
  return(result)

}
