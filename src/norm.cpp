#include <Rcpp.h>
using namespace Rcpp;


//' Calculate \eqn{L2}-norm
//'
//' This function calculates the \eqn{L_2}-norm of a vector.
//'
//' @param x a vector.
//' @return the \eqn{L_2}-norm of x.
//' @export
// [[Rcpp::export]]
double l2_norm(SEXP x){
  NumericVector  xx(x);
  double re = 0;
  for(int i = 0;i<xx.size();i++){
    re += pow(xx[i],2);
  }
  re = sqrt(re);
  return re;
}

//' Calculate \eqn{L1}-norm
//'
//' This function calculates the \eqn{L_1}-norm of a vector.
//'
//' @param x a vector.
//' @return the \eqn{L_1}-norm of x.
//' @export
// [[Rcpp::export]]
double l1_norm(SEXP x){
  NumericVector  xx(x);
  double re = 0;
  for(int i = 0;i<xx.size();i++){
    re += fabs(xx[i]);
  }
  return re;
}




