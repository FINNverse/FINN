# Predict from a FINN model

Predict from a FINN model

## Usage

``` r
# S3 method for class 'finn_class'
predict(
  object,
  env,
  disturbance = NULL,
  patches = 100L,
  patch_size = 0.1,
  init_cohort = NULL,
  device = c("cpu", "gpu"),
  return_cohorts = FALSE,
  debug = FALSE,
  ...
)
```

## Arguments

- object:

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

- device:

  (`character(1)`)  
  Should the simulation run on the CPU or the GPU (Graphics card).
  Support is only available for NVIDIA GPUs.

- return_cohorts:

  Controls whether the raw per-cohort state is returned in addition to
  the aggregated site output. Storing cohorts every timestep is
  expensive, so the default is `FALSE` (none). Use `TRUE` (or `"all"`)
  for every timestep, `"last"` for the final timestep only, or an
  integer vector of timesteps to store just those. Recorded cohorts
  appear as `$long$cohort` / `$wide$cohort`.

- debug:

  (`logical(1)`)  
  Debug modus or not. If `TRUE`, individual tree states are stored.

- ...:

  Advanced options forwarded to the internal simulator, chiefly
  `batchsize` (split the sites into batches of this many to cap memory
  for very large runs). The default processes all sites in one batch and
  is what almost all users want.

## Value

A named list of predictions. `$long$site` (and `$wide$site`) give the
site-level results (columns `siteID`, `year`, `species`, `variable`,
`value`). When `return_cohorts` is set, `$long$cohort` / `$wide$cohort`
add the per-cohort state (`dbh`, `trees`, `species`, growth `g`,
mortality `m`, ...) for the requested timesteps.

## Details

Simulate from a (fitted) FINN model. This is an S3 method for the
[stats::predict](https://rdrr.io/r/stats/predict.html) generic, so it is
dispatched as `predict(model, ...)`.
