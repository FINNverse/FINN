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

## See also

[`conditionalEffects()`](https://finnverse.github.io/FINN/reference/conditionalEffects.md),
[`summary.finn_class()`](https://finnverse.github.io/FINN/reference/summary.finn_class.md),
which prints these effects

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

# one number per driver and process: the average slope over the observed range,
# which is the linear-model coefficient a mechanistic process would report and
# the closest thing a neural-network process has to one
ace <- averageConditionalEffects(m)
head(ace$growth)
#>   species variable  mean_effect
#> 1       1     temp -0.018280701
#> 2       1     prec -0.009924146
#> 3       5     temp -0.015855194
#> 4       5     prec  0.008324687
#> 5      11     temp -0.028679522
#> 6      11     prec -0.006263988
# with the stand-state inputs of each process (dbh, light, ...) as well
head(averageConditionalEffects(m, env_only = FALSE)$mortality)
#>   species    variable   mean_effect
#> 1       1         dbh -5.155978e-06
#> 2       1      growth -1.288468e-01
#> 3       1       light -2.999704e-03
#> 4       1 (Intercept) -5.143399e-02
#> 5       1        temp  2.621655e-04
#> 6       1        prec -2.621650e-04
```
