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

You can learn more about FINN in the accompanying paper [Pichler & Käber
(2026) - Inferring processes within dynamic forest models using hybrid
modelling](https://doi.org/10.1111/2041-210x.70347), published in
*Methods in Ecology and Evolution*.

## Installation

Install the released version of FINN from CRAN:

``` r

install.packages("FINN")
```

Or the development version from
[GitHub](https://github.com/FINNverse/FINN):

``` r

# install.packages("devtools")
devtools::install_github("FINNverse/FINN")
```

## Introduction

Start with [Introduction to
FINN](https://finnverse.github.io/FINN/articles/A-Introduction_to_FINN.html),
then browse the rest:

| Vignette | What it covers |
|----|----|
| [Introduction to FINN](https://finnverse.github.io/FINN/articles/A-Introduction_to_FINN.html) | the model, its processes, and a first simulation |
| [Parameterising FINN from Ellenberg indicator values](https://finnverse.github.io/FINN/articles/B-Succession_demo.html) | deriving species niches from Ellenberg values; succession along a climate gradient |
| [Preparing your data for FINN](https://finnverse.github.io/FINN/articles/C-Data_preparation.html) | going from a raw tree list to FINN’s input tables |
| [Fitting FINN to forest inventory data](https://finnverse.github.io/FINN/articles/D-Fit_to_FIA.html) | calibrating on US FIA data, with a held-out test |
| [Mortality: a binomial response and a NN process](https://finnverse.github.io/FINN/articles/E-Mortality.html) | choosing the right likelihood, and scoring it honestly |

The same pages ship with the package — `vignette(package = "FINN")`
lists them.
