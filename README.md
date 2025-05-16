# R6BinSeg - An Object-Oriented R Package for Binary Segmentation via Modern C++

## This package is a demo intended for educational and presentation purposes only!!!

## Goals

✅ Provide an unoptimised version of BinSeg based on Rcpp/RcppArmadillo and implement the $\mathcal{O}(1)$ L2 cost function.  
Note: The current version only returns a user-specified K-partition of the time series without storing smaller partitions and their costs.  
⬜ Optimise BinSeg - not necessary to loop through all possible splits.  
⬜ Implement basic model selection criteria.  
⬜ Wrap these in an R6-based R package.  
⬜ Implement other popular cost functions to be used in BinSeg (L1, regression discontinuity, AR(k), kernel-based, etc).  
⬜ Extend BinSeg take various cost functions (e.g., based on a common class "Cost").  

## OOP Interface 

Users first create a C++ "Cost" object based on the time series. This object handles the computation of segment costs. 

```
costObj = createCostObj(tsMat, costFunc = "L2") #tsMat: a time series matrix
```
Then, initialise a (R6) BinSeg object, based on costObj. BinSeg does not rely on a specific distance method but instead uses an abstract "Cost" object.
```
BinSegObj = BinSeg$new(costObj = costObj)
```
The following methods are supported:

- $fit(): Perform binary segmentation.
```
BinSegObj$fit(minK = 2, maxK = 12, criterion = c("AIC", "BIC"), ...) 
```
- $predict(): Return the optimal segmentation based on the AIC (by default) or an user-provided K; must be called after BinSegObj$fit().
```
BinSegObj$predict(K = NULL, criterion = "AIC")
```
- $plot(): Plot the optimal segmentation based on the AIC (by default) or an user-provided K; must be called after BinSegObj$fit().
```
BinSegObj$plot(K = NULL, whichDim = NULL, criterion = "AIC") 
```
