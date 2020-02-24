#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


//' @title Perform matrix convolution
//' @param A Matrix A
//' @param B Matrix B
//' @return A convolved with B
// [[Rcpp::export]]
arma::mat convolve3cpp(arma::mat A, arma::mat B) {
  
  arma::mat D = arma::conv2(A, B, "same");
  return D;
}



