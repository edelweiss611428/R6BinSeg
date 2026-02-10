library("R6BinSeg")

segSize = 50
p_outl = 0.4
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


costL2 = Cost_L2$new(matrix(x))
bs_l2 <- BinSeg$new(costL2, minSize = 5)
bs_l2$fit()


par(mfrow = c(1,2))
plot(x, type = "l", main = "binSeg (L1)")
abline(v = bs$predict(5), col = "red", lwd = 2)
plot(x, type = "l", main = "binSeg (L2)")
abline(v = bs_l2$predict(5), col = "red", lwd = 2)
par(mfrow = c(1,1))
