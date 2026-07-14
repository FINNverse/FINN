# Conditional effects of a FINN model

Computes (and caches on the model) the conditional effects — the local
derivatives of each demographic process (growth, mortality,
regeneration) with respect to its inputs, per process and per species.
This is the shared primitive behind
[`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md),
[`averageConditionalEffects()`](https://finnverse.github.io/FINN/reference/averageConditionalEffects.md)
and the analytical variable importance in
[`summary.finn_class()`](https://finnverse.github.io/FINN/reference/summary.finn_class.md).

## Usage

``` r
conditionalEffects(
  model,
  env = NULL,
  init_cohort = NULL,
  env_autoscale = TRUE,
  sim_seed = 42L
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

## Value

an object of class `FINNconditionalEffects` (also cached on
`model$conditional_effects`).
