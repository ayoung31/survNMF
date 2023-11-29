#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat f_grad_rcpp(arma::mat& X, arma::mat& W, arma::mat& H, arma::mat& beta, 
             float alpha, arma::mat& y, arma::mat& delta, float theta, int j){
  int k = H.n_rows;
  arma::mat grad(k,1);
  NumericVector p1 (k);
  NumericVector p2 (k);
  p1 = W.t()*W*H.col(j-1) - W.t()*X.col(j-1);
  p2 = as_scalar(1 - (exp(beta.t()*H.col(j-1))/sum((y>=y[j-1]) % exp(beta.t()*H).t()))) * beta;
  
  grad = theta*((1-alpha)*p1-alpha*delta[j-1]*p2);
  
  return grad;
}