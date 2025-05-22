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
      Named("cp") = start + 1,
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
List binSegCpp2(Rcpp::XPtr<Cost> Xptr, const int& maxNRegimes) {

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

  arma::Row<double> cost(maxNRegimes); //Final cpd results - will be updated
  arma::Row<int> changePoints = arma::Row<int>{Rcpp::as<int>(cpd0["cp"])}; //Final cpd results - will be updated

  cost[0] = Xnew.effEvalCpp(0, nr);
  cost[1] = cpd0["err"];
  int idx = 2;
  int nRegimes = 2;
  List gbestCp = cpd0;
  List gsbestCp; // global
  List gtbestCp; // global third best
  double gsmaxGain = -std::numeric_limits<double>::infinity(); // global
  double gtmaxGain = -std::numeric_limits<double>::infinity();

  while(nRegimes < maxNRegimes){

    List tbestCp; //temporary best Cp
    List tsbestCp; //temporary 2nd best Cp
    double tmaxGain = -std::numeric_limits<double>::infinity();
    double tsmaxGain = -std::numeric_limits<double>::infinity();
    List cpdi;
    double gain;


    for(int i = 0; i < 2; i++){
      // here we expect that the best split will never have size 1
      if(i == 0){

        cpdi = miniOptCpp(Xnew, Rcpp::as<int>(gbestCp["start"]), Rcpp::as<int>(gbestCp["cp"]));

        if(!cpdi["valid"]){
          continue;
        }

        gain = Rcpp::as<double>(gbestCp["lErr"]) - Rcpp::as<double>(cpdi["err"]);

      } else {

        cpdi = miniOptCpp(Xnew, Rcpp::as<int>(gbestCp["cp"]), Rcpp::as<int>(gbestCp["end"]));

        if(!cpdi["valid"]){
          continue;
        }

        gain = Rcpp::as<double>(gbestCp["rErr"]) - Rcpp::as<double>(cpdi["err"]);

      }


      if(gain > tmaxGain){

        tsmaxGain = tmaxGain;
        tmaxGain = gain;

        tsbestCp = tbestCp;
        tbestCp = cpdi;

      } else { //only 2 scenarios - will be optimised later

        tsmaxGain = gain;
        tsbestCp = cpdi;

      }

    }

    nRegimes++;
    // a > b & c > d
    if(tmaxGain > gsmaxGain){ // if (a > c)

      gbestCp = tbestCp; //max1 = a
      cost[idx] = cost[idx-1] - tmaxGain;

      if(tsmaxGain > gsmaxGain){ // if (b > c)

        gtbestCp = gsbestCp; //max3 = c
        gtmaxGain = gsmaxGain;

        gsbestCp = tsbestCp;//max2 = b
        gsmaxGain = tsmaxGain;

      } else{

        // gsbestCp = gsbestCp // no change - just a note
        // gsmaxGain = gsmaxGain // no change - just a note

        if(tsmaxGain > gtmaxGain){ //b > d?

          gtbestCp = tsbestCp; //max3 = b
          gtmaxGain = tsmaxGain;

        } else{
          // gtbestCp = gtbestCp // no change - just a note
          // gtmaxGain = gtmaxGain // no change - just a note
        }
      }
    } else{ //if (a < c)
      gbestCp = gsbestCp; //max1 = c
      cost[idx] = cost[idx-1] - gsmaxGain;

      if(tmaxGain > gtmaxGain){ //if (a > d)

        gsbestCp = tbestCp; // max2 = a
        gsmaxGain = tmaxGain;

        if(tsmaxGain > gtmaxGain){ //if (b >d)

          gtbestCp = tsbestCp; //max3 = b
          gtmaxGain = tsmaxGain;

        } else {

          // gtbestCp = gtbestCp // no change - just a note
          // gtmaxGain = gtmaxGain // no change - just a note

        }
      } else {

        gtbestCp = tbestCp; //max3 = a
        gtmaxGain = tmaxGain;

        gsbestCp = gtbestCp; //max2 = d
        gsmaxGain = gtmaxGain;

      }
    }
    idx++;

    changePoints = arma::join_rows(changePoints, arma::Row<int>{Rcpp::as<int>(gbestCp["cp"])});

  }



  return List::create(
    Named("second") = changePoints,
    Named("cost")  = cost
   );
}


