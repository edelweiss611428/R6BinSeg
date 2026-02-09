#ifndef COST_H
#define COST_H

#include <RcppArmadillo.h>

using namespace Rcpp;

// ============================================================
// Base class
// ============================================================

class CostBase {
public:
  virtual ~CostBase() {}
  virtual double eval(int start, int end) const = 0;
  virtual int size() const = 0;   //
};

RCPP_EXPOSED_CLASS(CostBase)

  // ============================================================
  // Multivariate L2 cost
  // ============================================================

  class Cost_L2 : public CostBase {

  public:
    arma::mat csX;
    arma::mat csXsq;
    int nr;

    Cost_L2(const arma::mat& inputMat);

    double eval(int start, int end) const override;

    int size() const override { return nr; }   //
  };

RCPP_EXPOSED_CLASS(Cost_L2)


class RCostClass : public CostBase {

private:
  Function cost_fun_;
  int n_;

public:
  RCostClass(Function cost_fun, int n);

  double eval(int start, int end) const override;

  int size() const override { return n_; }
};

RCPP_EXPOSED_CLASS(RCostClass)
#endif
