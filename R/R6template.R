#' @importFrom R6 R6Class
#' @export

BinSegL2 <- R6Class(
  "BinSegL2",

  public = list(
    nr = NULL,
    cp = NULL,
    tsMat = NULL,
    k = NULL,
    fitted = FALSE,

    initialize = function(tsMat, k) {
      self$tsMat = tsMat
      self$nr = nrow(tsMat)
      self$k = k
      print("You have created a BinSegL2 object!")
    },

    fit = function() {
      self$cp = fastBinSegCpp(self$tsMat, self$k)$changePoints
      self$fitted = TRUE
    },

    plot = function() {

      if(ncol(self$tsMat == 1)){

        if(!self$fitted){
          ts.plot(self$tsMat, xlab = "X")
          warning("Should run the fit() method first if want changepoints plotted!")
        } else{

          ts.plot(self$tsMat, xlab = "X",
                  main = "Binary Segmentation")

          sortedRegimes = c(sort(self$cp), self$nr)

          for(i in self$k:1){

            lines(self$tsMat[1:sortedRegimes[i]], col = i+1)

          }

        }

      } else {
        print("Currently does not support high-dimensional plots!")
      }

    }

  ),
)

