# R6BinSeg - Object-Oriented Interface for Binary Segmentation via Modern C++

## Description

`R6BinSeg` is a **demonstration R package** illustrating how to implement
**binary segmentation (BinSeg)** using a clean **object-oriented design**
based on:

- **C++ cost objects** (via Rcpp / RcppArmadillo), and
- a **high-level R6 interface** for data segmentation, model selection, and visualisation.

The package is intended **solely for educational and presentation purposes**.
It prioritises **clarity of design** over performance, completeness, or API
stability.

The central design principle is:

> **Binary segmentation should depend on an abstract cost interface, not on a specific cost function.**

This is an extension of [`AnimalCrossing`](https://github.com/edelweiss611428/AnimalCrossing).
