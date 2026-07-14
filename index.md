# Forest Informed Neural Networks (FINN)

FINN is an R package designed for modular, dynamic vegetation (forest)
models. Modularity is achieved by implementing all components and the
platform in R and by allowing users to pass their own processes to FINN
as R functions. Nevertheless, FINN is highly performant (has a low
runtime) because it is implemented in torch for R, which allows access
to highly optimised routines and enables FINN to run on the GPU.
Consequently, FINN scales extremely well with the number of sites.
Furthermore, FINN can seamlessly integrate neural networks into
processes, or even replace them entirely (hybrid modelling).
Additionally, FINN’s parameters can be optimised directly using
stochastic gradient descent.

## Internal workings of FINN

You can learn more about FINN in our preprint [Pichler & Käber, 2025 -
Inferring processes within dynamic forest models using hybrid
modeling](https://arxiv.org/abs/2508.01228)

## Installation

You can install the development version of FINN from
[GitHub](https://github.com/) with:

``` r

# install.packages("devtools")
devtools::install_github("FINNverse/FINN")
```

## Introduction

An introduction to FINN can be found in the package vignette, or
[online](https://github.com/FINNverse/FINN/blob/main/vignettes/Introduction-to-FINN.qmd)
