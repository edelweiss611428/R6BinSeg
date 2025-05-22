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
      Named("cp") = start + 1
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
      Named("cp") = cp
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

  if(nr == 0){
    stop("There can be no changepoint when there is only one observation.");
  } else if (nr == 1 || maxNRegimes == 2){
    return cpd0;
  }

  arma::Row<int> changePoints = arma::Row<int>{Rcpp::as<int>(cpd0["cp"])};
  arma::Row<int> regimes = arma::join_rows(arma::Row<int>{0}, changePoints, arma::Row<int>{nr});
  arma::Row<double> updatedErrs = arma::Row<double>{Rcpp::as<double>(cpd0["lErr"]), Rcpp::as<double>(cpd0["rErr"])};

  int nRegimes = 2;
  List bestCp;
  int bestSplit;
  arma::Row<double> currentErrs = updatedErrs;

  while(nRegimes < maxNRegimes){

    double maxGain = -std::numeric_limits<double>::infinity();

     for(int i = 0; i < nRegimes; i++){

      List cpdi = miniOptCpp(Xnew, regimes[i], regimes[i+1]);
      if(!cpdi["valid"]){
        continue;
      }

      double gain = currentErrs[i] - Rcpp::as<double>(cpdi["err"]);
      if(gain > maxGain){
        maxGain = gain;
        bestCp = cpdi;
        bestSplit = i;
      }
    }

    nRegimes++;
    updatedErrs.set_size(nRegimes);
    updatedErrs.subvec(bestSplit, bestSplit+1) = arma::Row<double>{Rcpp::as<double>(bestCp["lErr"]), Rcpp::as<double>(bestCp["rErr"])};

    for(int j = 0; j < bestSplit; j++){
      updatedErrs[j] = currentErrs[j];
    }

    for(int j = bestSplit+2; j < nRegimes; j++){
      updatedErrs[j] = currentErrs[j-1];
    }
    changePoints = arma::join_rows(changePoints, arma::Row<int>{Rcpp::as<int>(bestCp["cp"])});
    changePoints = arma::sort(changePoints);
    regimes = arma::join_rows(arma::Row<int>{0}, changePoints, arma::Row<int>{nr});
    currentErrs = updatedErrs;
  }

  return List::create(
    Named("Regimes") = regimes
  );


}


