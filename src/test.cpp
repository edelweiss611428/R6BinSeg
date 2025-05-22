#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
inline arma::mat getCumsumCpp(const arma::mat& X) {

  int nr = X.n_rows;
  int nc = X.n_cols;

  // Create a matrix with one extra row (initialized to zeros)
  arma::mat cumsumMat(nr + 1, nc, arma::fill::zeros);

  // Compute cumulative sum starting from the second row
  cumsumMat.rows(1, nr) = arma::cumsum(X, 0);

  return cumsumMat;

}



class Cost {
private:
  arma::mat X;
  arma::mat csX; //cumsum(X)
  arma::mat csXsq; //cumsum(Xsq)

public:
  int nr;

  Cost(const arma::mat& inputMat) { //initialise a Cost object
    X = inputMat;
    csX = getCumsumCpp(inputMat);
    csXsq = getCumsumCpp(arma::pow(inputMat, 2));
    nr = X.n_rows;
  }


  double effEvalCpp(int start, int end) const {  //use precomputation

    if(start == end - 1){
      return(0.0);
    }

    int len = end - start;
    double sumErrXsq =  arma::sum(csXsq.row(end) - csXsq.row(start));
    double sqErrsumX =  std::pow(arma::norm(csX.row(end) - csX.row(start),2), 2);

    return sumErrXsq - sqErrsumX/len;

  }

};

RCPP_MODULE(costModule) {  //to use the class Cost in R
  class_<Cost>( "Cost")
  .constructor<arma::mat>()
  .method( "effEvalCpp", &Cost::effEvalCpp)
  ;
}
// [[Rcpp::export]]
SEXP createCostObj(const arma::mat& X) {
  Rcpp::XPtr<Cost> ptr(new Cost(X), true);
  return ptr;
}


inline List miniOptCpp(const Cost& Xnew, const int& start, const int& end) {

  int len = end - start;

  if(len == 1){
    return List::create(
      Named("valid") = false
    );
  } else if(len == 2){
    return List::create(
      Named("valid") = true,
      Named("err") = 0,
      Named("lErr") = 0,
      Named("rErr") = 0,
      Named("start") = start,
      Named("cp") = start+1,
      Named("end") = end
    );
  }

  double minErr = std::numeric_limits<double>::infinity();
  int cp;
  int tempCp;
  double err;
  double lErr;
  double rErr;
  double minlErr;
  double minrErr;

  for(int i = 0; i < (len-1); i++){
    tempCp = start+i+1;
    lErr = Xnew.effEvalCpp(start,tempCp);
    rErr = Xnew.effEvalCpp(tempCp,end);
    err = lErr + rErr;

    if(err < minErr){
      minErr = err;
      minlErr = lErr;
      minrErr = rErr;
      cp = tempCp;
    }
  }

  return List::create(
    Named("valid") = true,
    Named("err") = minErr,
    Named("lErr") = minlErr,
    Named("rErr") = minrErr,
    Named("start") = start,
    Named("cp") = cp,
    Named("end") = end
  );

}


// [[Rcpp::export]]
List binSegCpp(Rcpp::XPtr<Cost> Xptr, const int& maxNRegimes) {

  Cost& Xnew = *Xptr;
  int nr = Xnew.nr;
  List cpd0 = miniOptCpp(Xnew, 0, nr);
  if(maxNRegimes > nr){
    stop("The maximum number of regimes must be less than or equal to the number of observations.");;
  }

  if(nr == 1){
    stop("There is no changepoint as there is only one observation!");
  } else if (maxNRegimes == 2){
    return cpd0;
  }


  IntegerVector changePoints(maxNRegimes-1);
  IntegerMatrix regimes(2, maxNRegimes);
  IntegerVector savedCps(maxNRegimes);

  LogicalVector visited(maxNRegimes, false);
  NumericVector currentErrs(maxNRegimes);
  NumericVector tempErrs(maxNRegimes);
  NumericVector leftErrs(maxNRegimes);
  NumericVector rightErrs(maxNRegimes);
  NumericVector gains(maxNRegimes);


  changePoints[0] = Rcpp::as<int>(cpd0["cp"]);

  currentErrs[0] = Rcpp::as<double>(cpd0["lErr"]);
  currentErrs[1] = Rcpp::as<double>(cpd0["rErr"]);

  //1d indexing for matrices
  regimes[0] = 0;
  regimes[1] = Rcpp::as<int>(cpd0["cp"]);
  regimes[2] = Rcpp::as<int>(cpd0["cp"]);
  regimes[3] = nr;

  int nRegimes = 2;
  int bestIdx;

  while(nRegimes < maxNRegimes){

    double maxGain = -std::numeric_limits<double>::infinity();
    List cpdi;

    for(int i = 0; i < nRegimes; i++){

      if(not visited[i]){
        visited[i] = true;
        cpdi = miniOptCpp(Xnew, regimes[2*i], regimes[2*i+1]);

        if(!cpdi["valid"]){
          continue;
        }

        savedCps[i] = Rcpp::as<int>(cpdi["cp"]);
        leftErrs[i] = Rcpp::as<double>(cpdi["lErr"]);
        rightErrs[i] = Rcpp::as<double>(cpdi["rErr"]);
        gains[i] = currentErrs[i] - Rcpp::as<double>(cpdi["err"]);

      }

      if(gains[i] > maxGain){
        maxGain = gains[i];
        bestIdx = i;
      }
    }

    changePoints[nRegimes-1] = savedCps[bestIdx];
    IntegerVector aIdx = IntegerVector::create(bestIdx, nRegimes);
    visited[aIdx] = false;
    currentErrs[aIdx] = NumericVector::create(leftErrs[bestIdx], rightErrs[bestIdx]);
    int tEnd = regimes[2*bestIdx+1];
    regimes[2*bestIdx+1] = savedCps[bestIdx];
    regimes[2*nRegimes] = savedCps[bestIdx];
    regimes[2*nRegimes+1] = tEnd;
    nRegimes++;

  }

  return List::create(
    Named("Regimes") = changePoints
  );


}


