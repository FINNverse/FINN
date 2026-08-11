# Apply stored z-standardization to environmental predictors

Re-applies the constants learned by
[`compute_env_scaling()`](https://finnverse.github.io/FINN/reference/compute_env_scaling.md)
to `env`. The transformation uses the *stored* mean/sd only and never
recomputes them from `env`, so calibration and prediction use an
identical transformation (the usual pitfall of
[`scale()`](https://rdrr.io/r/base/scale.html) inside a model formula is
avoided).

## Usage

``` r
apply_env_scaling(env, scaling)
```

## Arguments

- env:

  A `data.frame`/`data.table` of raw environmental data.

- scaling:

  The `data.frame` returned by
  [`compute_env_scaling()`](https://finnverse.github.io/FINN/reference/compute_env_scaling.md),
  or `NULL` (in which case `env` is returned unchanged). For a model fit
  with `env_autoscale = TRUE` this is `model$env_scaling`.

## Value

`env` with the scaled predictor columns, as a `data.table`.

## See also

[`compute_env_scaling`](https://finnverse.github.io/FINN/reference/compute_env_scaling.md),
[`fit`](https://finnverse.github.io/FINN/reference/fit.md)

## Examples

``` r
# reproduce the standardization a model fit with env_autoscale = TRUE used:
env <- data.frame(siteID = 1:3, year = 1L, temp = c(4, 6, 8), prec = c(700, 850, 1000))
scaling <- compute_env_scaling(env)
apply_env_scaling(env, scaling)
#>    siteID  year  temp  prec
#>     <int> <int> <num> <num>
#> 1:      1     1    -1    -1
#> 2:      2     1     0     0
#> 3:      3     1     1     1
```
