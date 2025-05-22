# R6BinSeg - An Object-Oriented R Package for Binary Segmentation via Modern C++

## This package is a demo intended for educational and presentation purposes only!!!

## Goals

✅ Provide an unoptimised version of BinSeg based on Rcpp/RcppArmadillo and implement the $\mathcal{O}(1)$ L2 cost function.  
**NOTE: The current version only returns a user-specified K-partition of the time series without storing smaller partitions and their costs.**   
✅ Optimise BinSeg - not necessary to loop through all possible splits.    
**NOTE: Some refactoring is needed.**   
⬜ Implement basic model selection criteria.  
✅ Wrap these in an R6-based R package.  
⬜ Implement other popular cost functions to be used in BinSeg (L1, regression discontinuity, AR(k), kernel-based, etc).  
⬜ Extend BinSeg take various cost functions (e.g., based on a common class "Cost").  



## Cureent OOP interface 

Users first initialise a (R6) BinSeg object, based on a time series matrix, the number of regimes (k), and a cost function. 
```
BinSegObj = BinSeg$new(tsMat = tsMat, k, costFunc = "L2")
```
The following methods are supported:

- $fit(): Perform binary segmentation.
```
BinSegObj$fit() 
```
- $plot(): Plot the k-partition; must be called after BinSegObj$fit(); currently only supports one-dimensional data.
```
BinSegObj$plot() 
```



## Planned OOP interface 

Users first create a C++ "Cost" object based on a time series. This object handles the computation of segment costs. 

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
BinSegObj$fit(minK = 2, maxK = 12, ...) 
```
- $predict(): Return the optimal segmentation based on the AIC (by default) or a user-provided K; must be called after BinSegObj$fit().
```
BinSegObj$predict(K = NULL, criterion = "AIC")
```
- $plot(): Plot the optimal segmentation based on the AIC (by default) or a user-provided K; must be called after BinSegObj$fit().
```
BinSegObj$plot(K = NULL, whichDim = NULL, criterion = "AIC") 
```
