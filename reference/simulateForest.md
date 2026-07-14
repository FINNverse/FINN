# Simulate

Simulate

## Usage

``` r
simulateForest(
  model,
  env,
  disturbance = NULL,
  patches = 100L,
  patch_size = 0.1,
  init_cohort = NULL,
  batchsize = NULL,
  device = c("cpu", "gpu"),
  debug = FALSE
)
```

## Arguments

- model:

  (`finn_class`)  
  Object of class `finn_class` created by
  [finn](https://finnverse.github.io/FINN/reference/finn.md).

- env:

  (`data.table|data.frame`)  
  Data with environmental covariates must be passed as `data.table` or
  `data.frame`.

- disturbance:

  (`data.table|data.frame`)  
  Data with disturbance rates must be passed as `data.table` or
  `data.frame`.

- patches:

  (`integer(1)`)  
  Number of patches.

- patch_size:

  (`numeric(1)`)  
  Patch size.

- init_cohort:

  (`CohortMat`)  
  Initial cohort matrix of class `CohortMat`, created by
  [CohortMat](https://finnverse.github.io/FINN/reference/CohortMat.md).

- batchsize:

  (`integer(1)`)  
  Batch size, model will be trained in random batch sizes of the data to
  preserve memory and improve convergence.

- device:

  (`character(1)`)  
  Should the model be fitted on the CPU or the GPU (Graphic card).
  Support is only for NVIDIA GPUs available.

- debug:

  (`logical(1)`)  
  Debug modus or not. If `TRUE`, individual tree states are stored.

## Details

Simulate from a fitted FINN model. This is a thin, backwards-compatible
alias for
[predict.finn_class](https://finnverse.github.io/FINN/reference/predict.finn_class.md);
new code should call `predict(model, ...)` directly.
