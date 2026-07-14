# Summarise a fitted FINN model

Prints, per process and per species, the environmental variable
importance and the average conditional effects. Both are derived from
the model's conditional effects, which are computed once and cached — so
if [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) was
already run (or a previous
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)),
the default path needs no further simulation.

## Usage

``` r
# S3 method for class 'finn_class'
summary(
  object,
  env = NULL,
  init_cohort = NULL,
  importance = c("ale", "permutation"),
  env_autoscale = TRUE,
  sim_seed = 42L,
  nperm = 20L,
  scale = TRUE,
  ...
)
```

## Arguments

- object:

  (`finn_class`)  
  fitted model.

- env, init_cohort:

  (`NULL` or data)  
  interpretation data; `NULL` uses the cached training data (see
  [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md)).

- importance:

  (`character(1)`)  
  `"ale"` (default) = analytical ALE-variance importance (cheap, from
  the cache); `"permutation"` = re-simulating permutation importance
  (see
  [`feature_importance()`](https://finnverse.github.io/FINN/reference/feature_importance.md)),
  cached on the model under the same input key.

- env_autoscale:

  (`logical(1)`)  
  see [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md).

- sim_seed:

  (`integer(1)`)  
  seed for the simulation (see
  [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md)).

- nperm:

  (`integer(1)`)  
  replicates for `importance = "permutation"`.

- scale:

  (`logical(1)`)  
  for `importance = "ale"`, divide `Var(ALE)` by the process x species
  rate variance so the importances are dimensionless and comparable
  across processes and species (Sobol-style). Default `TRUE`.

- ...:

  passed through (e.g. to
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  for the permutation option).

## Value

invisibly, a list with `importance`, `average_conditional_effects`, and
`method`.
