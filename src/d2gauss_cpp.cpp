#include <Rcpp.h>
using namespace Rcpp;


//' @title Generate 2D Gaussian distribution
//' @param x description of first dimensions (e.g. ppm values)
//' @param y description of second dimensions (e.g. Hz values)
//' @param mu distribution's center position using coordinate in x and y dimension
//' @param sigma distribution's spread (standard deviation) in x and y dimension
//' @return Numeric vector of coordinates
// [[Rcpp::export]]
Rcpp::NumericVector d2gauss_cpp(NumericVector x, NumericVector y, NumericVector mu, NumericVector sigma) {
  static const double pi = 3.141592653589;
  
  double sprod = sigma(0) * sigma(1);
  
  NumericVector expo = (-1)*((pow((x-mu(0)), 2.0)/(2*pow(sigma(0), 2.0))) + (pow((y-mu(1)), 2.0)/(2*pow(sigma(1), 2.0))));
  
  double basis = (1/(2*pi*sprod));
  
  NumericVector out = basis * exp(expo);
  
  return out;
}