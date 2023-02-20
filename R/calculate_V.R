#' Hessian Matrix
#'
#' This function is used to calculate the Hessian matrix of specified loss function.
#'
#' @param BS a matrix, the columns of which should be the predicted value
#' of the orthogonal B-spline base at the corresponding point for continuous
#' variables or standardized data for discrete variables.
#' @param KM a vector, the Kaplan-Meier weights vector.
#' @param n the number of observations.
#' @return The Hessian matrix of weighted least squares loss function respect to \eqn{\theta}.
#' @export
#' @details The weighted least squares loss function is defined as:
#' \deqn{\ell_n(\theta) = \frac{1}{2}\sum_{i = 1}^nw_i[Y_{(i)}-\sum_{j = 1}^pZ_{(i)j}\sum_{k=1}^{q_n}\theta_{jk}B_k(Z_{(i)j})]^2,}
#' where \eqn{w_i} is the Kaplan-Meier weight for the \eqn{ith} order statistics \eqn{Y_{(i)}} of \eqn{Y_i's}, \eqn{Z_{(i)j}} is the
#' value of \eqn{jth} covariate correspond to \eqn{Y_{(i)}}.
#' @note See the example of \code{\link{BMD}} function for an example.
calculate_V <- function(BS, KM, n){
  W <- diag(KM, n, n)
  V <- t(BS) %*% W %*% BS
  # b <- t(BS) %*% (KM * ny)
  return(V)
}
