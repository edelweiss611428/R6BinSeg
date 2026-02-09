#include "Cost.h"

// [[Rcpp::depends(RcppArmadillo)]]

// ========================================================
// Utility
// ========================================================

arma::mat getCumSumCpp(const arma::mat& X){

  int nr = X.n_rows;
  int nc = X.n_cols;

  arma::mat cumsumMat(nr + 1, nc, arma::fill::zeros);
  cumsumMat.rows(1, nr) = arma::cumsum(X, 0);

  return cumsumMat;
}

// ========================================================
// Cost_L2
// ========================================================

Cost_L2::Cost_L2(const arma::mat& inputMat){

  csX   = getCumSumCpp(inputMat);
  csXsq = getCumSumCpp(arma::pow(inputMat, 2));

  nr = inputMat.n_rows;
}

double Cost_L2::eval(int start, int end) const {

  if(start >= end-1) return 0.0;

  int len = end - start;

  return arma::sum(csXsq.row(end) - csXsq.row(start)) -
    std::pow(arma::norm(csX.row(end) - csX.row(start), 2), 2) / len;
}


// ========================================================
// RCostClass
// ========================================================

RCostClass::RCostClass(Function cost_fun, int n)
  : cost_fun_(cost_fun), n_(n) {}

double RCostClass::eval(int start, int end) const {
  return as<double>(cost_fun_(start, end));
}

// ========================================================
// Rcpp module
// ========================================================

RCPP_MODULE(cost_module){

  class_<CostBase>("CostBase");

  class_<Cost_L2>("Cost_L2")
    .derives<CostBase>("CostBase")
    .constructor<arma::mat>()
    .method("eval", &Cost_L2::eval)
    .method("size", &Cost_L2::size);

  class_<RCostClass>("RCostClass")
    .derives<CostBase>("CostBase")
    .constructor<Function, int>()
    .method("eval", &RCostClass::eval)
    .method("size", &RCostClass::size);
}
