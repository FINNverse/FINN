# Average conditional effects of a FINN model

Summarises the conditional effects into a per process x species x
variable average marginal effect: the mean local derivative
(`mean_effect`, an approximate linear effect). Derived cheaply from the
cached conditional effects.

## Usage

``` r
averageConditionalEffects(
  model,
  env = NULL,
  init_cohort = NULL,
  env_autoscale = TRUE,
  sim_seed = 42L,
  env_only = TRUE
)
```

## Arguments

- model:

  (`finn_class`)  
  fitted model.

- env:

  (`data.table|data.frame`\|`NULL`)  
  env covariates; `NULL` uses the cached training env.

- init_cohort:

  (`CohortMat`\|`NULL`)  
  initial cohort; `NULL` uses the cached training init_cohort.

- env_autoscale:

  (`logical(1)`)  
  see [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md).

- sim_seed:

  (`integer(1)`)  
  seed for the state-harvesting simulation (see
  [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md)).

- env_only:

  (`logical(1)`)  
  if `TRUE` (default) report only the environmental predictors of each
  process; if `FALSE` also include the stand-state inputs (dbh, light,
  growth) each process depends on.

## Value

a named list (one entry per process) of data.frames with columns
`species, variable, mean_effect`.
