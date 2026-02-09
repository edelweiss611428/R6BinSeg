#' @useDynLib R6BinSeg, .registration = TRUE
#' @importFrom Rcpp loadModule
#' @import methods
#' @export



.onLoad <- function(libname, pkgname) {
  Rcpp::loadModule("cost_module", TRUE)
  Rcpp::loadModule("binseg_module", TRUE)
}
