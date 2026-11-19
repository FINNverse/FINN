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
  [`predict()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  to rescale raw `env` internally; `FALSE` treats `env` as already on
  the model scale.

- ...:

  passed to
  [`predict()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  (e.g. `patches`, `patch_size`).

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

## See also

[`summary.finn_class()`](https://finnverse.github.io/FINN/reference/summary.finn_class.md)
(`importance = "permutation"`),
[`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md)

## Examples

``` r
library(data.table)
# a small model fitted on three Oregon FIA sites (?fia_obs_dt, ?fit); the NN
# variant is the same call with `hidden = c(10, 10)` on the process
sites <- 1:3
Nsp   <- max(fia_obs_dt$species)
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
FINN.seed(1)
m <- finn(
  N_species            = Nsp, recruits_dbh = 12.9,
  growth_process       = createProcess(~ temp + prec, FINN::growth,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~ temp + prec, FINN::mortality,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~ temp + prec, FINN::regeneration,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE))
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06, plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.768)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%

# permutation importance: how much the simulated rate of a process shifts
# when a driver is shuffled between sites and years. `nperm` is small here
# to keep the example quick; use more for a real reading.
fi <- feature_importance(m, nperm = 2L, seed = 1L)
head(fi$growth)
#>    species variable importance
#> 1        1     temp  1.2552243
#> 12       1     prec  1.0058492
#> 13       2     prec  1.2301807
#> 2        2     temp  0.4590797
#> 3        3     temp  1.2019306
#> 14       3     prec  1.0795486
# `method = "sobol"` gives the dimensionless Sobol-style score instead
```
