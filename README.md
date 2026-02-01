# R6BinSeg
Object-Oriented Binary Segmentation with a Polymorphic Cost Interface
*(Currently in progress)*

## Overview

`R6BinSeg` is a demonstration R package illustrating how to implement
Binary Segmentation (BinSeg) using a clean object-oriented architecture
based on:

- Modern C++ cost objects (via Rcpp/RcppArmadillo)
- A high-level R6 interface for segmentation workflows
- A single abstract cost interface shared by both C++ and R`defined costs

The package is intended solely for educational and presentation purposes.
It prioritises clarity of design over performance, completeness, or API
stability.

---

## Design principle

Binary segmentation should depend on an abstract cost interface, not on a
specific cost function.

The Binary Segmentation algorithm is implemented once, against an abstract
C++ base class `CostBase`. Any cost function—whether written in C++ or R—can be
used as long as it satisfies this interface.

---

## Core idea

At the C++ level, the BinSeg engine depends only on a reference to a `CostBase`
object and is completely agnostic to the concrete implementation.

This enables:

- Native C++ costs (fast, compiled, reusable)
- R-defined costs (flexible, user-extensible)

to be used interchangeably through one common C++ base class. From the perspective of the BinSeg algorithm, all costs are just `CostBase`.

This package extends the design pattern introduced in: [`AnimalCrossing`](https://github.com/edelweiss611428/AnimalCrossing)
