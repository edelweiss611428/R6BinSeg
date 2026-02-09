#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include "Cost.h"


IntegerVector binSegPredCpp(const IntegerVector& bkps,
                            const NumericVector& cost,
                            double penalty = 0){

  // # nocov start
  if(penalty <0){
    stop("Penalty should be non-negative!"); //tested in R
  }
  // # nocov end

  NumericVector penCost = clone(cost);

  for(int i = 0; i < cost.size(); i++){
    penCost[i] = cost[i] + penalty*i;
  }

  int minIdx = which_min(penCost);

  if(minIdx == 0){
    return IntegerVector();
  }

  return bkps[Range(0, minIdx-1)];

}


// ========================================================
//                  Utility: Segment class
// ========================================================

struct Segment {
  int start;
  int end;
  bool valid;
  int cp;
  double gain;
  double lErr; //error of left segment
  double rErr; //error of right segment
  double err; //total error

  // Max-heap: higher gain has higher priority
  bool operator<(const Segment& other) const {
    return gain < other.gain;
  }
};



// ========================================================
//                Utility: miniOptHeapCpp
// ========================================================

inline Segment miniOptHeapCpp(const CostBase& costModule, int start, int end, int minLen,
                              int minSize = 1,  int jump = 1, double totalErr = -1) {

  int len = end - start;

  if(totalErr < 0){
    totalErr = costModule.eval(start, end);
  }

  if(len < minLen){

    return Segment{start, end, true, start,
                   -std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()};

  } else if(len == minLen){

    int cp = start + len/2;
    double lErr = costModule.eval(start, cp);
    double rErr = costModule.eval(cp, end);
    double err = lErr + rErr;
    return Segment{start, end, true, cp,
                   totalErr - err, //gain
                   lErr,
                   rErr,
                   err};
  }

  double minErr = std::numeric_limits<double>::infinity();
  int cp;
  int tempCp = start;
  double err;
  double lErr;
  double rErr;
  double minlErr;
  double minrErr;

  auto allBkps = arma::regspace<arma::ivec>(start, jump, end); //all breakpoinbts
  allBkps = allBkps(arma::find((allBkps - start >= minSize) % (end - allBkps >= minSize)));
  //(start, End]
  for(int i = 0; i < allBkps.n_elem; i++){

    tempCp = allBkps(i);
    lErr = costModule.eval(start,tempCp);
    rErr = costModule.eval(tempCp,end);
    err = lErr + rErr;

    if(err < minErr){
      minErr = err;
      minlErr = lErr;
      minrErr = rErr;
      cp = tempCp;
    }
  }

  return Segment{start, end, true, cp,
                 totalErr - minErr, //gain
                 minlErr,
                 minrErr,
                 minErr};
}



class binSegCpp {

private:
  const CostBase& cost_;
  int minSize;
  int jump;
  int minLen;
  int nSamples;

public:
  IntegerVector bkpsVec;
  NumericVector costVec;

  binSegCpp(SEXP costSEXP, int nSamples_, int minSize_, int jump_)
    : cost_(Rcpp::as<const CostBase&>(costSEXP)),
      minSize(minSize_),
      jump(jump_),
      nSamples(nSamples_) {

    if(minSize < 1) stop("minSize >= 1");
    if(jump < 1) stop("jump >= 1");

    int k = std::ceil(double(minSize) / jump);
    minLen = 2 * k * jump;

    if(nSamples < minLen)
      stop("Too few observations");
  }

  // ================= fit =================

  void fit(){

    int nr = nSamples;
    int maxNRegimes = std::floor(double(nr) / minSize);

    NumericVector cost(maxNRegimes);
    IntegerVector changePoints(maxNRegimes - 1);

    double initCost = cost_.eval(0, nr);

    Segment seg0 =
      miniOptHeapCpp(cost_, 0, nr, minLen, minSize, jump, initCost);

    cost[0] = initCost;

    std::priority_queue<Segment> heap;
    heap.push(seg0);

    int idx = 0;
    int nRegimes = 2;

    while(nRegimes <= maxNRegimes){

      Segment best = heap.top();
      heap.pop();

      if(best.gain < 0) break;

      changePoints[idx] = best.cp;
      cost[idx+1] = cost[idx] - best.gain;
      idx++;

      heap.push(miniOptHeapCpp(cost_, best.start, best.cp,
                               minLen, minSize, jump, best.lErr));

      heap.push(miniOptHeapCpp(cost_, best.cp, best.end,
                               minLen, minSize, jump, best.rErr));

      nRegimes++;
    }

    bkpsVec = changePoints[Range(0, idx-1)];
    costVec = cost[Range(0, idx)];
  }

  // ================= predict =================

  IntegerVector predict(double penalty){
    return binSegPredCpp(bkpsVec, costVec, penalty);
  }

  double eval(int start, int end){
    return cost_.eval(start,end);
  }
};



RCPP_EXPOSED_CLASS(binSegCpp)

  RCPP_MODULE(binseg_module){

    class_<binSegCpp>("binSegCpp")
      .constructor<SEXP,int,int,int>()   // cost, nSamples, minSize, jump
      .method("fit", &binSegCpp::fit)
      .method("predict", &binSegCpp::predict)
      .method("eval", &binSegCpp::eval)
      .field("bkpsVec", &binSegCpp::bkpsVec)
      .field("costVec", &binSegCpp::costVec);
  }
