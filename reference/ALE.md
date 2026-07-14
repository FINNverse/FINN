# Accumulated local effect plots

Calculates accumulated local effects (ALE) for the three processes

## Usage

``` r
ALE(
  model,
  env = NULL,
  init_cohort = NULL,
  env_autoscale = TRUE,
  sim_seed = 42L,
  plot = TRUE,
  process = NULL,
  scale = FALSE,
  ...
)
```

## Arguments

- model:

  (`finn_class`)  
  Model object of class `finn_class`.

- env:

  (`data.table|data.frame`)  
  Environmental covariates for which the ALE should be calculated. If
  `NULL` (default) the training `env` cached by
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) is used.

- init_cohort:

  (`CohortMat`)  
  Initial cohort of class `CohortMat`. If `NULL` (default) the training
  init_cohort cached by
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) is used
  (or bare ground if none was cached).

- env_autoscale:

  (`logical(1)`)  
  If `TRUE` (default) `env` is assumed to be on the raw (unscaled) scale
  and the model's stored `env_scaling` is applied internally before the
  effects are computed, mirroring how the model was fitted (see
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)'s
  `env_autoscale`). Set `FALSE` if `env` is already on the model scale,
  or for a model fitted without autoscaling.

- sim_seed:

  (`integer(1)`)  
  Seed applied via
  [`FINN.seed()`](https://finnverse.github.io/FINN/reference/FINN.seed.md)
  before the state-harvesting simulation, so the (stochastic)
  conditional effects / ALE are reproducible and cacheable. `NULL`
  disables seeding.

- plot:

  (`logical(1)`)  
  If `TRUE` (default) an ALE plot is drawn via
  [`plot.FINNale()`](https://finnverse.github.io/FINN/reference/plot.FINNale.md)
  (rows = processes, columns = environmental predictors, one coloured
  line per species).

- process:

  (`character`\|`NULL`)  
  If given (one of `"growth"`, `"mortality"`, `"regeneration"`) only
  that process is plotted. `NULL` (default) plots all three.

- scale:

  (`logical(1)`)  
  If `TRUE` each curve is divided by the SD of its process x species
  rate, yielding dimensionless, comparable effects (the curve's variance
  then equals the Sobol-style normalised importance). Default `FALSE`.

- ...:

    
  Not supported yet.
