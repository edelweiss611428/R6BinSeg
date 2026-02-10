library("R6BinSeg")

segSize = 50
p_outl = 0.3
x <- c(
  rnorm(segSize, -5, 0.3) + rt(segSize, 2)*sample(c(0,1), segSize, T, p = c(p_outl, 1- p_outl)),
  rnorm(segSize, -3, 0.3) + rt(segSize, 2)*sample(c(0,1), segSize, T, p = c(p_outl, 1- p_outl)),
  rnorm(segSize, -1, 0.3)  + rt(segSize, 2)*sample(c(0,1), segSize, T,p = c(p_outl, 1- p_outl)),
  rnorm(segSize, 1, 0.3)  + rt(segSize, 2)*sample(c(0,1), segSize, T, p = c(p_outl, 1- p_outl)),
  rnorm(segSize, 3, 0.3)  + rt(segSize, 2)*sample(c(0,1), segSize, T, p = c(p_outl, 1- p_outl)),
  rnorm(segSize, 5, 0.3) + rt(segSize, 2)*sample(c(0,1), segSize, T, p = c(p_outl, 1- p_outl))
)

f <- function(s,e){
  sum(abs(x[(s+1):e] - median(x[(s+1):e])))
}

cost <- RCostClass$new(f, length(x))

bs <- BinSeg$new(cost, minSize = 5)
bs$fit()
plot(x, type = "l", main = "Binary segmentation")
abline(v = bs$predict(5), col = "red", lwd = 2)


system.file("examples/custom_L1_costFunc.R", package = "R6BinSeg") |> source()
