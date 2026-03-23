#' Cost_L1 cost module
#'
#' @description Computes the L1 cost for a segment of a matrix.
#' @details The L1 cost is the sum of absolute deviations from the median.
#'
#' @name Cost_L1
#' @export
NULL

#' Cost_L2 cost module
#'
#' @description Computes the L2 cost (squared error) for a segment of a matrix.
#' @details The L2 cost is the sum of squared deviations from the mean.
#'
#' @name Cost_L2
#' @export
NULL

#' R-based custom cost module
#'
#' @description Provides a bridge for user-supplied cost functions implemented in R.
#' @details Allows flexible experimentation with any R-based cost function.
#'
#' @name RCostClass
#' @export
NULL

#' Binary segmentation with modular cost objects
#'
#' @docType class
#' @details Performs greedy binary segmentation using a cost object derived from
#' \code{CostBase}. Cost functions may be implemented either in R or in C++,
#' sharing a common interface. The full segmentation path is computed in C++
#' for efficiency, while changepoints are selected post hoc via a user-specified
#' penalty.
#'
#' @field nr Integer. Number of observations inferred from the cost object.
#' @field fitted Logical. Indicates whether \code{fit()} has been called.
#' @field cost_path Numeric vector. Cost values along the binary segmentation path.
#' @field bkps_path Integer vector. Breakpoints along the full segmentation path.
#' @field costObject CostBase. User-supplied cost object.
#' @field engine External pointer. Underlying C++ binary segmentation engine.
#' @field minSize Integer. Minimum allowed segment length.
#' @field jump Integer. Grid spacing for candidate split locations.
#'
#' @importFrom R6 R6Class
#' @export
BinSeg = R6Class(
  "BinSeg",

  public = list(
    nr = NULL,
    fitted = FALSE,
    cost_path = NULL,
    bkps_path = NULL,
    costObject = NULL,
    engine = NULL,
    minSize = NULL,
    jump = NULL,

    #' @description Initialise BinSeg
    #'
    #' @param costObject A CostBase object (e.g. \code{Cost_L1}, \code{Cost_L2} or \code{RCostClass}).
    #' @param minSize Minimum allowed segment length.
    #' @param jump Grid spacing for candidate split locations.
    #'
    #' @return Invisibly return `NULL`.
    #'
    #'
    initialize = function(costObject, minSize = 1, jump = 1) {
      if (missing(costObject)) stop("Provide a costObject (Cost_L1 / Cost_L2 / RCostClass / ...).")
      if (minSize < 1) stop("minSize must be >= 1")
      if (jump < 1) stop("jump must be >= 1")

      self$costObject = costObject
      self$nr = costObject$size()
      self$minSize = as.integer(minSize)
      self$jump = as.integer(jump)

      self$engine = new(binSegCpp, costObject, self$nr, self$minSize, self$jump)

      message("You have created a BinSeg object!")
      invisible(NULL)
    },

    #' @description Fit binary segmentation
    #'
    #' Computes the full binary segmentation solution path in C++.
    #' Stores the breakpoint path and cost path internally.
    #'
    #' @return Invisibly returns the \code{BinSeg} object.
    #'
    fit = function() {
      self$engine$fit()
      self$bkps_path = self$engine$bkpsVec
      self$cost_path = self$engine$costVec
      self$fitted = TRUE
      invisible(self)
    },

    #' @description Predict changepoints
    #'
    #' Selects changepoints from the fitted solution path using a linear penalty.
    #'
    #' @param penalty Non-negative penalty applied to the cost path.
    #'
    #' @return Integer vector of changepoint locations.
    predict = function(penalty = 0) {
      if (!self$fitted) stop("Call fit() first.")
      self$engine$predict(penalty)
    }
  )
)
