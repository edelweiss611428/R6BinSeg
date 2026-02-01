# R6BinSeg - Object-Oriented Interface for Binary Segmentation via Modern C++

## Description

`R6BinSeg` is a **demonstration R package** illustrating how to implement
**binary segmentation (BinSeg)** using a clean **object-oriented design**
based on:

- **C++ cost objects** (via Rcpp / RcppArmadillo), and
- a **high-level R6 interface** for fitting, model selection, and plotting.

The package is intended **solely for educational and presentation purposes**.
It prioritises **clarity of design** over performance, completeness, or API
stability.

The central design principle is:

> **Binary segmentation should depend on an abstract cost interface, not on a specific cost function.**

---

## Goals

- Demonstrate BinSeg implemented in C++ with an $\mathcal{O}(1)$ L2 cost  
- Show how segmentation logic can be decoupled from cost computation  
- Provide an R6 interface suitable for teaching, prototyping, and extension  

---

## Design overview

The package separates responsibilities into two layers:

| Layer | Responsibility |
|-----|---------------|
| **Cost object (C++)** | Compute segment cost `cost(start, end)` |
| **BinSeg (R6)** | Run segmentation, select models, plot results |

This separation allows the same BinSeg logic to work with **any cost function**
that conforms to the expected interface.

---
