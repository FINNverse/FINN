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
`model$conditional_effects`): a list with one element per process
(`growth`, `mortality`, `regeneration`), each a list with one element
per species holding the inputs the process was evaluated at (`X`), the
local derivatives with respect to them (`grads`) and the predicted rate
(`pred`).

## See also

[`averageConditionalEffects()`](https://finnverse.github.io/FINN/reference/averageConditionalEffects.md)
for one number per driver,
[`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) for the
curves

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

# the local slope of each process with respect to each of its inputs, at the
# states the fitted model passes through: a list per process and species
ce <- conditionalEffects(m)
names(ce)
#> [1] "processes"     "env"           "env_autoscale" "sim_seed"     
#> [5] "key"          
# species 1, growth: the inputs, the slopes and the predicted rate
str(ce$growth[[1]], max.level = 1)
#>  NULL
```
