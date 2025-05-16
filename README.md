R6BinSeg - An Object-Oriented R Package for Binary Segmentation via Modern C++

Goals:

- Provide an unoptimised version of BinSeg based on Rcpp/RcppArmadillo and implement the O(1) $\mathcal{L}_2$ cost function.
- Implement other popular cost functions to be used in BinSeg (L1, regression discontinuity, AR(k), kernel-based, etc).
- Extend BinSeg take various cost functions (e.g., based on a common class "Cost").
- Wrap these in an R package based on R6.
