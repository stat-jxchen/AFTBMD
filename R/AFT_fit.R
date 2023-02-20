#' Fit Additive AFT Model
#'
#' This function fits an additive model for given standardized data with BMD algorithm.
#'
#' @param data a list containing:\tabular{ll}{
#'    \code{y} \tab a numeric vector recording logarithm survival time. \cr
#'    \tab \cr
#'    \code{x} \tab the standardized covariates matrix with all entries between -1 and 1.\cr
#'    \tab \cr
#'    \code{delta} \tab a numeric vector that consists of status indicator, normally 0=alive, 1=dead. \cr
#' }
#' @param mset a vector, the number of intervals for each continuous variable to expand the B-spline.
#' @param pset a vector, the degree of the B-spline expanded for each continuous variable.
#' @param lambda1 a number, penalty parameters for linear effects.
#' @param lambda2 a number, penalty parameters for nonlinear effects.
#' @param lambda1s a vector, penalty parameters for linear effects, used to select the best hyperparameter.
#' @param lambda2s a vector, penalty parameters for nonlinear effects, used to select the best hyperparameter.
#' @param penalty a character, which indicates the penalty function, and can be one of "glasso", "gscad" and "gmcp".
#' @param M a number, the maximum times of iteration of BMD algorithm.
#' @param epsilon eps a number, the tolerance error.
#' @return a list containing:\tabular{ll}{
#'    \code{ny} \tab a vector, the logarithm survival time after sorting in ascending order.\cr
#'    \tab \cr
#'    \code{nx} \tab a vector, x sorted by rows corresponding to \code{ny}.\cr
#'    \tab \cr
#'    \code{BS} \tab a matrix, the columns of which should be the predicted value
#'     of the orthogonal B-spline base at the corresponding point for continuous
#'     variables or standardized data for discrete variables.\cr
#'    \tab \cr
#'    \code{theta_est} \tab a vector, the estimated values of parameters.\cr
#'    \tab \cr
#'    \code{best_lambda} \tab a vector, the best value for \code{lambda1} and \code{lambda2} according to EBIC criterion.\cr
#'    \tab \cr
#'    \code{group_re} \tab a list, the estimated values of parameters for each group.\cr
#'    \tab \cr
#'    \code{model_struct} \tab a vector, identify the effect type of each component.\cr}
#' @examples
#' mset <- rep(4,30)
#' pset <- rep(3,30)
#' lam_gscad1 <- seq(0.03, 0.04, 0.002)
#' lam_gscad2 <- seq(0.03, 0.04, 0.002)
#' re <- AFT_fit(SimuData, mset, pset, lambda1s = lam_gscad1, lambda2s = lam_gscad2,
#'               penalty = "gscad", M = 300)
#' # coefficients estimate
#' re$theta_est
#' # best pair of regulation parameters
#' re$best_lambda
#' # model
#' re$model_structure
#' @export
AFT_fit <- function(data, mset, pset, lambda1, lambda2, lambda1s,
                    lambda2s, penalty, M = 300, epsilon = 1e-4){

  y <- data$y
  x <- data$x
  # data preparation
  ny <- sort(y)
  nx <- x[order(y), ]
  KM <- KMweight(data)
  # record some key values
  d <- ncol(x)
  n <- nrow(x)
  continue_num <- length(mset)
  dummy_num <- d - continue_num
  dfset <- c(mset + pset, rep(1,dummy_num))
  # construct orthogonal Bsplines for continuous variable
  BS <- numeric(0)
  for(i in 1:continue_num){
    p <- pset[i]
    df <- dfset[i]
    order <- p + 1
    nIknots <- df - order
    Boundary.knots <- range(-1, 1)
    knots <- seq.int(from = (-1), to = 1,
                     length.out = nIknots + 2)[-c(1, nIknots + 2)]
    BS1 <- splines2::bSpline(nx[, i], degree = p, knots = knots,
                             Boundary.knots = Boundary.knots, intercept = TRUE)
    BS2 <- nx[, i] * bar_bs(BS1, d)
    BS <- cbind(BS, BS2)
  }
  # if there exists dummy variables, combine them
  if(dummy_num != 0){
    BS <- cbind(BS,data$x[,-c(1:continue_num)])
  }

  # calculate the max eigenvalue of corresponding submatrix for each group
  h <- calculate_V(BS, KM, n)
  hj <- hj2 <- rep(0, d)
  gr_loc <- c(0,cumsum(dfset)) # used to find the location of groups
  for(i in 1:continue_num){
    start <- gr_loc[i]+1
    hj[i] <- h[start,start]
    start <- gr_loc[i]+2
    end <- gr_loc[i+1]
    h1 <- h[start:end, start:end]
    hj2[i] <- max(eigen(h1)$value)
  }
  if(dummy_num != 0){
    hj[(continue_num+1):d] <- diag(h)[gr_loc[continue_num+2]:gr_loc[length(gr_loc)]]
  }
  # fit
  if(missing(lambda1s) | missing(lambda2s)){
    if(missing(lambda1) | missing(lambda2)){
      stop('Arguments "lambda1s" or "lambda2s" are missing, then both "lambda1"
           and "lambda2" must be passed in simultaneously.')
    }else{
      theta_est <- BMD(BS, KM, ny, hj, hj2, dfset, M, lambda1, lambda2,
                   penalty, epsilon)
      best_lambda <- NULL
    }
  }else{
    len1 <- length(lambda1s)
    len2 <- length(lambda2s)
    theta_mat <- array(0, c(len1, len2, sum(dfset)))
    Loss <- matrix(0, len1, len2)
    for(j1 in 1:len1){
      for(j2 in 1:len2){
        theta_mat[j1,j2,] <- BMD(BS, KM, ny, hj, hj2, dfset, M, lambda1s[j1],
                                  lambda2s[j2], penalty, epsilon)
        Loss[j1,j2] <- sum((ny - BS %*% theta_mat[j1, j2, ])^2 * KM)/2
      }
    }
    theta_est_obj <- selection(theta_mat, n, Loss)
    theta_est <- theta_est_obj$theta_est
    id <- theta_est_obj$id
    best_lambda <- c(lambda1s[id[1]], lambda2s[id[2]])
  }

  group_re <- split(theta_est, rep(1:length(dfset),dfset))
  model_struc <- character(d)
  for (i in 1:d) {
    if(sum(group_re[[i]] != 0) == 0){
      model_struc[i] = "unimportant"
    }else if(sum(group_re[[i]][-1] != 0) == 0){
      model_struc[i] = "linear"
    }else{
      model_struc[i] = "nonlinear"
    }
  }
  re <- list(ny = ny,
             nx = nx,
             BS = BS,
             theta_est = theta_est,
             best_lambda = best_lambda,
             group_re = group_re,
             model_struc = model_struc)
  return(re)
}




