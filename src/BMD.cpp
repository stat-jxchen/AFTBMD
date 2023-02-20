#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// select the larger
double mymax(double x, double y){
  return (x < y) ? y : x;
}

// l2-norm of x
double l2_norm(NumericVector x){
  double re = 0;
  for(int i = 0;i<x.size();i++){
    re += pow(x[i],2);
  }
  re = sqrt(re);
  return re;
}

// l1-norm of x
double l1_norm(NumericVector x){
  double re = 0;
  for(int i = 0;i<x.size();i++){
    re += fabs(x[i]);
  }
  return re;
}

// soft threshold function
NumericVector S(NumericVector x, double lambda){
  NumericVector re(x.size());
  double l2_x = l2_norm(x);
  if((l2_x>=1e-15)&(l2_x<DBL_MAX)){
    re = mymax(1 - lambda/l2_x, 0) * x;
  }
  return re;
}

// glasso penalty
NumericVector glasso(NumericVector theta, NumericVector u, double h, double lambda){
  int l = theta.size();
  NumericVector cj(l),re(l);
  cj = theta-u/h;
  re = S(cj,lambda/h);
  return re;
}

// dot product of two vector
double dot(NumericVector c, NumericVector d){
  int nc = c.size();
  double e = 0;
  for (int i = 0; i < nc; i++){
    e = e + c[i] * d[i];
  }
  return e;
}

// calculate the gradient of weighted least squares loss function
NumericVector grad(NumericMatrix BS, NumericVector KM, NumericVector ny,
                   NumericVector theta, int start, int end, bool first){
  NumericVector temp(BS.nrow());
  int n = KM.size();
  if(first){
    NumericVector re(1);
    for(int i = 0;i<n;i++){
      if(KM[i] != 0){
        temp = BS(i,_);
        re[0] = re[0] - KM[i]*(ny[i]-dot(temp,theta))*temp[start];
      }
    }
    return re;
  }else{
    NumericVector re(end-start);
    for(int i = 0;i<n;i++){
      if(KM[i] != 0){
        temp = BS(i,_);
        re = re - KM[i]*(ny[i]-dot(temp,theta))*temp[Range(start+1,end)];
      }
    }
    return re;
  }
}

// gmcp penalty
NumericVector gmcp(NumericVector theta, NumericVector u, double h, double lambda){
  int l = theta.size();
  NumericVector cj(l),re(l);
  cj = theta-u/h;
  double gamma = mymax(3.7, 2/h);
  if(l2_norm(cj) <= lambda * gamma){
    re = S(h/(h-1/gamma) * cj, lambda/(h-1/gamma));
  }else{
    re = cj;
  }
  return re;
}

// gscad penalty
NumericVector gscad(NumericVector theta, NumericVector u, double h, double lambda){
  int l = theta.size();
  NumericVector cj(l),re(l);
  cj = theta-u/h;
  double gamma = 3.7;
  if(l2_norm(cj) <= lambda + lambda/h){
    re = S(cj, lambda/h);
  }else if(l2_norm(cj) <= lambda * gamma){
    re = (h * (gamma - 1) - lambda * gamma/l2_norm(cj)) / (h * gamma - h - 1) * cj;
  }else{
    re = cj;
  }
  return re;
}


//' Blockwise Majorization Descent Algorithm
//'
//' This function returns the optimal solution vector according to the BMD algorithm.
//'
//' @param BS a matrix, the columns of which should be the predicted value
//' of the orthogonal B-spline base at the corresponding point for continuous
//' variables or standardized data for discrete variables.
//' @param KM a vector, the Kaplan-Meier weights vector.
//' @param ny a vector, the logarithm survival time after sorting.
//' @param h1 a vector, the first entry of Hessian matrix according to every group.
//' @param h2 a vector, the max eigenvalue of Hessian matrix according to every group except the first entry.
//' In particular, if a group only contains a single column, such as discrete variable, the
//' corresponding element should be assign to 0.
//' @param group a vector, the columns of \code{BS} in each group.
//' @param M a number, the maximum times of iteration.
//' @param lambda1 a number, the penalty factor corresponding to the first entry of every group.
//' @param lambda2 a number, the penalty factor corresponding to remaining entries of every group.
//' @param penalty a character, which indicates the penalty function, and can be one of "glasso", "gscad" and "gmcp".
//' @param eps a number, the tolerance error.
//' @return The optimal solution vector which minimizes the weighted least squares loss function
//' with specific penalty.
//' @details (ToDo) Refer to the method in initial paper and the method in this paper.
//' @export
//' @examples
//' x <- SimuData$x
//' n <- nrow(x)
//' d <- ncol(x)
//' y <- SimuData$y
//' M <- 200 # maximum iter time
//' m <- 4
//' p <- 3 # the degree of the piecewise polynomial
//' df <- m + p # degree of freedom
//' dfset <- rep(df, d) # record the df of every group
//' nIknots <- m - 1 # the number of knots
//' Boundary.knots <- range(-1, 1)
//' knots <- seq.int(from = (-1), to = 1, length.out = nIknots + 2)[-c(1, nIknots + 2)]
//' ny <- sort(y)
//' nx <- x[order(y), ]
//' KM <- KMweight(SimuData)
//' BS <- numeric(0)
//' for(i in 1:d){
//'   BS1 <- splines2::bSpline(nx[, i], degree = p, knots = knots,
//'   Boundary.knots = Boundary.knots, intercept = TRUE)
//'   BS2 <- nx[, i] * bar_bs(BS1, d)
//'   BS <- cbind(BS, BS2)}
//'
//' h <- calculate_V(BS, KM, n)
//' hj <- hj2 <- rep(0, d)
//' gl <- c(0,cumsum(dfset))
//' for(i in 1:d){
//'   start <- gl[i]+1
//'   hj[i] <- h[start,start]
//'   start <- gl[i]+2
//'   end <- gl[i+1]
//'   h1 <- h[start:end, start:end]
//'   hj2[i] <- max(eigen(h1)$value)}
//' # rough tuning parameter grid
//' lambda1 <- seq(0.01, 0.014, 0.002)
//' lambda2 <- seq(0.01, 0.014, 0.002)
//' lg.lambda1 <- length(lambda1)
//' lg.lambda2 <- length(lambda2)
//' # storage the tuning results
//' theta_est <- array(0, c(lg.lambda1, lg.lambda2, sum(dfset)))
//' lf_value <- matrix(0, lg.lambda1, lg.lambda2)
//' for(j1 in 1:lg.lambda1){
//'   for(j2 in 1:lg.lambda2){
//'     theta_est[j1,j2,] <- BMD(BS, KM, ny, hj, hj2, dfset, M, lambda1[j1],
//'                              lambda2[j2], "gscad", 1e-4)
//'     lf_value[j1,j2] <- sum((ny - BS %*% theta_est[j1, j2, ])^2 * KM)/2
//'   }
//' }
//' # select the best estimator and record the best regularization parameters pair
//' theta_obj <- selection(theta_est, n, lf_value)
//' theta_hat <- theta_obj$theta_est
//' id <- theta_obj$id
//' lambda_pair <- c(lambda1[id[1]], lambda2[id[2]])
// [[Rcpp::export]]
NumericVector BMD(SEXP BS, SEXP KM, SEXP ny, SEXP h1, SEXP h2, SEXP group,
                      int M, double lambda1, double lambda2, String penalty,
                      double eps){

  NumericMatrix bs(BS);
  NumericVector km(KM), nyy(ny), hj1(h1), hj2(h2);
  NumericVector theta(bs.ncol()), theta_old(bs.ncol());
  NumericVector u1(1);
  NumericVector g(group);
  int gl = g.size();
  NumericVector gcum(gl+1);
  gcum[0] = 0;
  for(int i = 0;i<gl;i++){
    gcum[i+1] = gcum[i] + g[i];
  }
  int start, end;
  double diff = 1;
  int k = 0;
  if(penalty == "glasso"){
    while((diff > eps) & (k < M)){
      for(int i = 0; i < gl; i++){
        start = gcum[i];
        end = gcum[i+1]-1;
        if((theta[start] != 0) || (k == 0)){
          u1 = grad(bs, km, nyy, theta, start, end, true);
          theta[Range(start,start)] = glasso(theta[Range(start,start)], u1, hj1[i], lambda1);
        }
        if(g[i] != 1){
          if((theta[start+1] != 0) || (k == 0)){
            NumericVector u2(g[i]-1);
            u2 = grad(bs, km, nyy, theta, start, end, false);
            theta[Range(start+1,end)] = glasso(theta[Range(start+1,end)],u2,hj2[i],lambda2);
          }
        }
      }
      diff = l2_norm(theta-theta_old);
      theta_old = clone(theta);
      k++;
    }
  }else if(penalty == "gmcp"){
    while((diff > eps) & (k < M)){
      for(int i = 0; i < gl; i++){
        start = gcum[i];
        end = gcum[i+1]-1;
        if((theta[start] != 0) || (k == 0)){
          u1 = grad(bs, km, nyy, theta, start, end, true);
          theta[Range(start,start)] = gmcp(theta[Range(start,start)], u1, hj1[i], lambda1);
        }
        if(g[i] != 1){
          if((theta[start+1] != 0) || (k == 0)){
            NumericVector u2(g[i]-1);
            u2 = grad(bs, km, nyy, theta, start, end, false);
            theta[Range(start+1,end)] = gmcp(theta[Range(start+1,end)],u2,hj2[i],lambda2);
          }
        }
      }
      diff = l2_norm(theta-theta_old);
      theta_old = clone(theta);
      k++;
    }
  }else if(penalty == "gscad"){
    while((diff > eps) & (k < M)){
      for(int i = 0; i < gl; i++){
        start = gcum[i];
        end = gcum[i+1]-1;
        if((theta[start] != 0) || (k == 0)){
          u1 = grad(bs, km, nyy, theta, start, end, true);
          theta[Range(start,start)] = gscad(theta[Range(start,start)], u1, hj1[i], lambda1);
        }
        if(g[i] != 1){
          if((theta[start+1] != 0) || (k == 0)){
            NumericVector u2(g[i]-1);
            u2 = grad(bs, km, nyy, theta, start, end, false);
            theta[Range(start+1,end)] = gscad(theta[Range(start+1,end)],u2,hj2[i],lambda2);
          }
        }
      }
      diff = l2_norm(theta-theta_old);
      theta_old = clone(theta);
      k++;
    }
  }
  return theta;
}





