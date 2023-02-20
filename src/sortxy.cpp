#include <Rcpp.h>
using namespace Rcpp;

// sort x by the order of y
NumericVector Rcpp_sort(NumericVector x, NumericVector y) {
  IntegerVector idx = seq_along(x) - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});
  return x[idx];
}


//' Sorted Estimation
//'
//' This function returns the sorted estimation at \code{pointx}.
//'
//' @param x a vector consisting of the covariate values.
//' @param y a vector consisting of the estimated response values corresponding to \code{x}.
//' @param pointx a vector consisting of the values of the covariates to be predicted.
//' @return the sorted estimation at \code{pointx}.
//' @details Denote \eqn{x_{(i)}} as the \eqn{ith} order statistics of \eqn{x},
//' and \eqn{y_{(i)}} is \eqn{y} value corresponding to \eqn{x_{(i)}}.
//' The estimator of \eqn{y} value of \eqn{pointx_j} is \eqn{y_{(i)}} if
//' \eqn{x_{(i)}\leq pointx_j < x_{(i+1)}}. If \eqn{x_{(n)}\leq pointx_j},
//' then the estimator is \eqn{y_{(n)}} and if \eqn{pointx_j < x_{(1)}} then the
//' estimator is \eqn{y_{(1)}}.
//' @export
//' @examples
//' sortxy(x = rnorm(100), y = rnorm(100), pointx = rnorm(100))
// [[Rcpp::export]]
NumericVector sortxy(SEXP x, SEXP y, SEXP pointx){
  NumericVector xx(x), yy(y), p(pointx);
  int n = xx.size();
  int l = p.size();
  NumericVector re(l),xx1(n),yy1(n);
  NumericVector xxc = clone(xx);
  xx1 = xxc.sort();
  yy1 = Rcpp_sort(yy,xx);
  for(int i=0;i<l;i++){
    if(p[i]<xx1[0]){
      re[i] = yy1[0];
    }else if(p[i]>=xx1[n-1]){
      re[i] = yy1[n-1];
    }else{
      for(int j=1;j<n;j++){
        if((p[i]>=xx1[j-1]) & (p[i]<xx1[j])){
          re[i] = yy1[j-1];
        }
      }
    }
  }
  return re;

}





