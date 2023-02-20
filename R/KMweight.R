#' Kaplan-Meier Weights.
#'
#' This function is used to calculate Kaplan-Meier weights vector.
#'
#' @param data a list, in which \code{x} is the covariates matrix, \code{y} is the logarithm survival time, and \code{delta} is the status indicator, normally 0=alive, 1=dead.
#' @return A vector denotes Kaplan-Meier weights vector.
#' @export
#' @examples
#' f1 <- function(x) return(x)
#' f2 <- function(x) return(-x)
#' f3 <- function(x) return(x^2)
#' mydata <- gdata(tau = 1, n = 200, d = 30, rho = 0.25,
#'                 fun_list = list(f1, f2, f3), col_index = 1:3)
#' KMweight(mydata)
KMweight <- function(data){
  y <- data$y
  delta <- data$delta
  delta <- delta[order(y)]
  n <- length(delta)
  KM <- rep(1/n, n)
  for(i in 1:(n-1)){
    if(delta[i] == 0){
      for(j in (i+1):n){
        KM[j] <- KM[j] + KM[i]/(n-i)
      }
      KM[i] <- 0
    }
  }
  if(delta[n] == 0){
    # finally all KM weights may add up to less than 1
    KM[n] <- 0
  }
  # KM weight corresponds to the weight after y is arranged from small to large
  return(KM)
}
