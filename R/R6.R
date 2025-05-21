# Load necessary libraries
library(R6)
library("Rcpp")
library("RcppArmadillo")
sourceCpp("./src/L2BinSeg.cpp")
sourceCpp("./src/L2BinSeg2.cpp")
miniOpt = function(Xnew, start, end){

  nr = end - start

  if(nr == 1){
    return(FALSE)
  } else if (nr == 2){
    return(list(cp = start + 1,
                err = 0,
                lErr = 0,
                rErr = 0))

  }

  Cv = numeric(nr-1)
  minErr = Inf

  for(i in 1:(nr-1)){

    lErr = Xnew$effEvalCpp(start,start+i)
    rErr = Xnew$effEvalCpp(start+i,end)
    Cv[i] = lErr + rErr
    # print(Cv[i])
    if(Cv[i] < minErr){
      cp = i
      minErr = Cv[i]
      minlErr = lErr
      minrErr = rErr

    }

  }

  return(list(cp = start + cp,
              err = minErr,
              lErr = minlErr,
              rErr = minrErr))
}



binSeg = function(Xnew, maxRegimes = 5, nr){

  cpd1 = miniOpt(Xnew, 0, nr)
  cp = cpd1$cp
  regimes = c(0, cp, nr)
  nregimes = 2
  newCost = c(cpd1$lErr, cpd1$rErr)

  repeat{

    currentCost = newCost
    maxG = -Inf

    for(i in 1:nregimes){

      cpdi = miniOpt(Xnew, regimes[i], regimes[i+1])
      if(is.logical(cpdi)){
        next
      }
      gain = currentCost[i] - cpdi$err
      if(gain > maxG){
        maxG = gain
        bestCpd = cpdi
        bestI = i
      }

    }

    cpnew = c(cp, bestCpd$cp)
    nregimes = nregimes  + 1
    newCost = numeric(nregimes)
    newCost[c(bestI, bestI+1)] = c(bestCpd$lErr, bestCpd$rErr)
    newCost[-c(bestI, bestI+1)] = currentCost[-bestI]

    cp = sort(cpnew)

    regimes = c(0, cp, nr)

    if(nregimes  == maxRegimes){
      break
    }

  }

  return(regimes)

}

binSegCpp(Xnew, k)

for(i in 1:1000){
  set.seed(i)
  N = 350
  k = 7
  p =  rep(1/7, 7)

  # Draw one sample
  counts =  as.vector(rmultinom(n = 1, size = N, prob = p))
  X = rnorm(N, rep(rnorm(k,0, 25), counts), 5)
  Xnew = createCostObj(as.matrix(X))
  cpd = as.vector(binSegCpp(Xnew, k)$Regimes)
  ts.plot(X, col = "red", main = "Binary Segmentation")
  lines(X[1:cpd[7]], col = "brown")
  lines(X[1:cpd[6]], col = "yellow")
  lines(X[1:cpd[5]], col = "blue")
  lines(X[1:cpd[4]], col = "green")
  lines(X[1:cpd[3]], col = "violet")
  lines(X[1:cpd[2]], col = "black")

  binseg = c(0, binsegRcpp::binseg_normal(X, 7)$splits[,4]$end)
  print(all.equal(cpd, sort(binseg)))
  gc()
  Sys.sleep(1)
}

N = 70000
k = 10
p =  rep(1/k, k)

# Draw one sample
counts =  as.vector(rmultinom(n = 1, size = N, prob = p))
X = rnorm(N, rep(rnorm(k,0, 25), counts), 5)
Xnew = createCostObj(as.matrix(X))
microbenchmark::microbenchmark(binSegCpp(Xnew, k),
                               binSegCpp2(Xnew, k),
                               binsegRcpp::binseg_normal(X, k))
