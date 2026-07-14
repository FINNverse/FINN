# Learn z-standardization for environmental predictors

Computes the centering and scaling constants (mean and standard
deviation) of every numeric environmental predictor in `env` (all
columns except the keys `siteID` and `year`), so they can be re-applied
unchanged to new data at prediction time. A predictor with (near-)zero
standard deviation is given `scale = 1` (centred only) to avoid division
by zero, mirroring `recipes::step_normalize`.

## Usage

``` r
compute_env_scaling(env)
```

## Arguments

- env:

  A `data.frame`/`data.table` with `siteID`, `year` and the
  environmental predictor columns.

## Value

A `data.frame` with columns `variable`, `center`, `scale`; or `NULL` if
there are no numeric predictors. This is the object stored on a fitted
model as `model$env_scaling` when it is fit with `env_autoscale = TRUE`.

## See also

[`apply_env_scaling`](https://finnverse.github.io/FINN/reference/apply_env_scaling.md)
