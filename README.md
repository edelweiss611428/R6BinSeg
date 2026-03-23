# R6BinSeg

**Object-Oriented Binary Segmentation with a Polymorphic Cost Interface**

*(Project currently in progress)*

## Overview

`R6BinSeg` is a demonstration R package illustrating how to implement
binary segmentation (BinSeg) using a clean object-oriented architecture
based on:

* Modern C++ cost objects (via Rcpp/RcppArmadillo)
* A high-level R6 interface for segmentation workflows
* A single abstract cost interface shared by both C++ and R-defined costs

The package is intended solely for educational and presentation purposes.
It prioritises clarity of design over performance, completeness, or API
stability.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("edelweiss611428/R6BinSeg")
```

## Example

### Simulated data with outliers

```r
library(R6BinSeg)

set.seed(256)

segSize <- 50
p_outl <- 0.4

x <- c(
  rnorm(segSize, -5, 0.3) + rt(segSize, 1) * sample(c(0,1), segSize, TRUE, prob = c(p_outl, 1 - p_outl)),
  rnorm(segSize, -3, 0.3) + rt(segSize, 1) * sample(c(0,1), segSize, TRUE, prob = c(p_outl, 1 - p_outl)),
  rnorm(segSize, -1, 0.3) + rt(segSize, 1) * sample(c(0,1), segSize, TRUE, prob = c(p_outl, 1 - p_outl)),
  rnorm(segSize,  1, 0.3) + rt(segSize, 1) * sample(c(0,1), segSize, TRUE, prob = c(p_outl, 1 - p_outl)),
  rnorm(segSize,  3, 0.3) + rt(segSize, 1) * sample(c(0,1), segSize, TRUE, prob = c(p_outl, 1 - p_outl)),
  rnorm(segSize,  5, 0.3) + rt(segSize, 1) * sample(c(0,1), segSize, TRUE, prob = c(p_outl, 1 - p_outl))
)
```
<img width="1736" height="822" alt="image" src="https://github.com/user-attachments/assets/0fcbc1da-db92-4513-92bb-63545d7dda51" />

---

### Custom cost function (L1, robust)

```r
# Define L1 cost
f <- function(s, e) {
  sum(abs(x[(s + 1):e] - median(x[(s + 1):e])))
}

# Wrap as cost object
cost <- RCostClass$new(f, length(x))

# Run BinSeg
bs <- BinSeg$new(cost, minSize = 5)
bs$fit()

# Extract change-points
cp_l1 <- bs$predict(10)
```

---

### Built-in cost (L2)

```r
costL2 <- Cost_L2$new(matrix(x))

bs_l2 <- BinSeg$new(costL2, minSize = 5)
bs_l2$fit()

cp_l2 <- bs_l2$predict(10)
```

---

### Visualisation

```r
par(mfrow = c(1, 2))

plot(x, type = "l", main = "BinSeg (L1)")
abline(v = cp_l1, col = "red", lwd = 2)

plot(x, type = "l", main = "BinSeg (L2)")
abline(v = cp_l2, col = "red", lwd = 2)

par(mfrow = c(1, 1))
```

<img width="1736" height="822" alt="image" src="https://github.com/user-attachments/assets/e67d6d02-71a7-4f69-8670-d30daf6ba8d6" />


---

## Key Idea

The core design principle of `R6BinSeg` is a **polymorphic cost interface**:

* C++ cost functions and R-defined cost functions share the same interface
* Binary Segmentation operates independently of the cost implementation
* Enables flexible experimentation with custom loss functions

---

## Status

This package is under active development and intended for demonstration and
educational purposes.
