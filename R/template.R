miniOptR6 <- R6Class(
  "miniOpt",
  public = list(
    n = NULL,
    cp = NULL,
    tsMat = NULL, #time series matrix
    fitted = FALSE,

    initialize = function(tsMat) {
      self$tsMat = tsMat
      self$n = nrow(tsMat)
      print("Hi, you have created a miniOpt object!")
    },


    fit = function() {
      self$cp = miniOptCpp(self$tsMat, 0,self$n)
      self$fitted = TRUE
    },


    plot = function() {
      if(!self$fitted){
        ts.plot(self$tsMat)
        warning("Should run the fit() method first if want changepoints plotted!")
      } else{
        ts.plot(self$tsMat)
        lines(self$tsMat[1:self$cp,], col = "red")
      }
    }

  ),
)



X = matrix(c(rnorm(2000,10),rnorm(2000,0)))
model1 = miniOptR6$new(X)
model1$fit()
model1$plot()
model1$plot()
