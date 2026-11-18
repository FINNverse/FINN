# Compute per-predictor environmental scaling

Computes the centering and scaling constants of every numeric
environmental predictor in `env` (all columns except the keys `siteID`
and `year`), so they can be re-applied unchanged to new data at
prediction time. Every supported mode is an affine map
`(x - center) / scale`, so the stored `(variable, center, scale)` table
is applied by
[`apply_env_scaling`](https://finnverse.github.io/FINN/reference/apply_env_scaling.md)
and inverted by ALE unchanged, whatever the mode. A predictor with
(near-)zero standard deviation is given `scale = 1` (centred only) to
avoid division by zero.

## Usage

``` r
compute_env_scaling(env, spec = TRUE)
```

## Arguments

- env:

  A `data.frame`/`data.table` with `siteID`, `year` and the
  environmental predictor columns.

- spec:

  how to scale each predictor. One of:

  - `TRUE` – z-standardise every predictor (`center = mean`,
    `scale = sd`); the default and the historical behaviour.

  - `FALSE` – handled by the caller (no scaling table computed).

  - a length-1 string applied to all predictors, or a vector/list with
    one entry per predictor. Entries may be:

    - `"auto"` z-standardise (mean / sd),

    - `"identity"` (alias `"none"`) leave unchanged (center 0, scale 1),

    - `"0to1"` min-max to `[0, 1]` (center min, scale range),

    - a **function** `f(x)` returning `list(center=, scale=)` – a custom
      affine scaler (e.g. robust:
      `function(x) list(center=median(x), scale=IQR(x))`). Non-affine
      transforms are intentionally not supported because ALE needs an
      invertible mapping. Unnamed vectors/lists are matched to
      predictors by position; named ones by variable name (predictors
      not named default to `"auto"`).

## Value

A `data.frame` with columns `variable`, `center`, `scale`; or `NULL` if
there are no numeric predictors. This is the object stored on a fitted
model as `model$env_scaling` when it is fit with `env_autoscale = TRUE`.

## See also

[`apply_env_scaling`](https://finnverse.github.io/FINN/reference/apply_env_scaling.md),
[`fit`](https://finnverse.github.io/FINN/reference/fit.md)
(`env_autoscale`),
[`finn_data`](https://finnverse.github.io/FINN/reference/finn_data.md)

## Examples

``` r
# the constants fit(env_autoscale = TRUE) stores on a model for this table
compute_env_scaling(fia_env_dt)
#>    variable      center      scale
#> 1      temp    6.866604   2.388455
#> 2   tempmax   25.362500   2.223931
#> 3   tempmin   -5.205000   4.158342
#> 4      prec 1103.250000 685.523918
#> 5  precseas   54.545718  15.704577
#> 6 precwarmq   90.955000  35.459819
# keep one predictor on its own scale
compute_env_scaling(fia_env_dt, spec = c(temp = "auto", prec = "identity"))
#>    variable    center     scale
#> 1      temp  6.866604  2.388455
#> 2   tempmax 25.362500  2.223931
#> 3   tempmin -5.205000  4.158342
#> 4      prec  0.000000  1.000000
#> 5  precseas 54.545718 15.704577
#> 6 precwarmq 90.955000 35.459819
# min-max scaling for every predictor, or a custom (robust) scaler for one
head(compute_env_scaling(fia_env_dt, spec = "0to1"), 3)
#>   variable      center  scale
#> 1     temp   0.4041666 11.725
#> 2  tempmax  20.0000000 11.400
#> 3  tempmin -12.5000000 15.800
robust <- function(x) list(center = median(x), scale = IQR(x))
head(compute_env_scaling(fia_env_dt, spec = list(temp = robust)), 3)
#>   variable    center    scale
#> 1     temp  6.239583 3.930208
#> 2  tempmax 25.362500 2.223931
#> 3  tempmin -5.205000 4.158342
```
