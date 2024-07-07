// [[Rcpp::depends(RcppArmadillo)]]
// my_functions.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;

// Function to calculate lambda
double get_lambda(std::vector<int> node, const arma::mat &Theta) {
  double lambda {0};
  if (node.size() == 1) {
    int i = node.at(0) - 1;
    lambda = exp(Theta(i,i));
  } else {
    int i = node.back() - 1;
    node.pop_back();
    double temp {0};
    for (int j : node) {
      temp += Theta(i,j - 1);
    }
    lambda = exp(Theta(i,i) + temp);
  }
  return lambda;
}

// Rcpp export wrapper
// [[Rcpp::export]]
double get_lambda_wrapper(SEXP nodeSEXP, SEXP ThetaSEXP) {
  std::vector<int> node = Rcpp::as<std::vector<int>>(nodeSEXP);
  const arma::mat& Theta = Rcpp::as<arma::mat>(ThetaSEXP);
  return get_lambda(node, Theta);
}