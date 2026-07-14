# Permutation feature importance for FINN demographic rates

Permutes each environmental predictor of a process and measures how
strongly the permutation shifts that process's predicted rate, **per
process and per species**. Unlike the analytical ALE-variance
importance, this re-simulates the model, so it captures the full
dynamical response (feedback through the stand state) rather than the
process response function alone. It does NOT use the conditional-effects
cache.

## Usage

``` r
feature_importance(
  model,
  env = NULL,
  init_cohort = NULL,
  nperm = 20L,
  method = c("rmse", "sobol"),
  seed = NULL,
  sim_seed = 42L,
  env_autoscale = TRUE,
  ...
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
  init cohort; `NULL` uses the cached training init_cohort.

- nperm:

  (`integer(1)`)  
  number of permutation replicates (default 20).

- method:

  (`character(1)`)  
  `"rmse"` or `"sobol"`.

- seed:

  (`integer(1)`\|`NULL`)  
  R RNG seed controlling which permutations are drawn.

- sim_seed:

  (`integer(1)`\|`NULL`)  
  torch seed for common random numbers across runs.

- env_autoscale:

  (`logical(1)`)  
  see [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md).
  `TRUE` (default) leaves
  [`predict()`](https://rdrr.io/r/stats/predict.html) to rescale raw
  `env` internally; `FALSE` treats `env` as already on the model scale.

- ...:

  passed to [`predict()`](https://rdrr.io/r/stats/predict.html) (e.g.
  `patches`, `patch_size`).

## Value

a named list (one per process) of data.frames with columns
`species, variable, importance`, sorted within species.

## Details

Two scorings via `method`:

- `"rmse"` — RMSE between the unpermuted and permuted rate, in units of
  that species' rate SD. Unbounded; larger = more important.

- `"sobol"` — total-effect estimator
  `0.5 * mean(MSE_shift) / Var(rate)`; dimensionless (only bounded in
  `[0, 1]` under independent predictors — FINN's climate predictors are
  usually correlated, so treat it as relative).

Predictors are read per-process from the model formulas, so processes
with different formulas get different variable sets. Common random
numbers (`sim_seed`, applied to the torch simulation RNG) make the
stochastic mortality/regeneration draws shared across the reference and
permuted runs, so a driver with no effect returns ~0.
