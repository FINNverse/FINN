# Summarise a fitted FINN model

Prints, per process and per species, the environmental variable
importance and the average conditional effects. Both are derived from
the model's conditional effects, which are computed once and cached — so
if [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) was
already run (or a previous
[`summary()`](https://rdrr.io/r/base/summary.html)), the default path
needs no further simulation.

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
  [`predict()`](https://rdrr.io/r/stats/predict.html) for the
  permutation option).

## Value

invisibly, a list with `importance`, `average_conditional_effects`, and
`method`.

## See also

[`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md),
[`averageConditionalEffects()`](https://finnverse.github.io/FINN/reference/averageConditionalEffects.md),
[`feature_importance()`](https://finnverse.github.io/FINN/reference/feature_importance.md)

## Examples

``` r
library(data.table)
# a small model fitted on two Oregon FIA sites (see ?fit)
sites <- 1:2
Nsp   <- max(fia_obs_dt$species)
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
FINN.seed(1)
m <- finn(N_species = Nsp, recruits_dbh = 12.9,
          growth_process       = createProcess(~ temp + prec, growth),
          mortality_process    = createProcess(~ temp + prec, mortality),
          regeneration_process = createProcess(~ temp + prec, regeneration_capped))
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 1L, patches = 4, patch_size = 0.06,
    plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00301, ba=8.77, trees=0.852, growth=1.58e+04, mortality=32.4, regeneration=0.677)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%

# which predictor matters for which process and species, and in which direction
s <- summary(m)
#> FINN model summary
#> ==================
#> 
#> ### GROWTH
#> 
#> Analytical ALE importance (rate-normalised) (species in columns):
#>  variable    sp1     sp2    sp3 sp4    sp5 sp6 sp7      sp8 sp9  sp10    sp11
#>      temp 0.2789 0.00418 0.5204   0 0.0673   0   0 0.048926   0 0.291 0.03872
#>      prec 0.0643 0.03167 0.0968   0 0.0217   0   0 0.000299   0 0.228 0.00434
#> 
#> Average conditional effects (mean; species in columns):
#>  variable      sp1      sp2     sp3      sp4      sp5     sp6     sp7      sp8
#>      temp -0.01504 -0.00548  0.0477  0.00101 -0.01324 -0.0063  0.0529  0.08218
#>      prec -0.00722  0.01509 -0.0206 -0.04514  0.00752 -0.0552 -0.0246 -0.00643
#>     sp9   sp10     sp11
#>  0.0262 0.0366 -0.01694
#>  0.0123 0.0324 -0.00567
#> 
#> ### MORTALITY
#> 
#> Analytical ALE importance (rate-normalised) (species in columns):
#>  variable   sp1    sp2 sp3 sp4    sp5 sp6 sp7    sp8 sp9 sp10  sp11
#>      temp 0.241 0.0475   0   0 0.0518   0   0 0.0589   0  0.2 0.263
#>      prec 0.241 0.0475   0   0 0.0518   0   0 0.0589   0  0.2 0.263
#> 
#> Average conditional effects (mean; species in columns):
#>  variable       sp1       sp2 sp3 sp4       sp5 sp6       sp7       sp8
#>      prec -0.000188 -0.000183   0   0 -0.000171   0  0.000222  0.000182
#>      temp  0.000188  0.000183   0   0  0.000171   0 -0.000222 -0.000182
#>        sp9      sp10      sp11
#>   0.000211  0.000216 -0.000186
#>  -0.000211 -0.000216  0.000186
#> 
#> ### REGENERATION
#> 
#> Analytical ALE importance (rate-normalised) (species in columns):
#>  variable     sp1    sp2      sp3    sp4    sp5      sp6   sp7    sp8    sp9
#>      prec 0.17037 0.0166 0.065128 0.1151 0.1499 0.010799 0.307 0.1571 0.0209
#>      temp 0.00513 0.0975 0.000859 0.0524 0.0753 0.000872 0.144 0.0128 0.2588
#>    sp10    sp11
#>  0.0918 0.05192
#>  0.2072 0.00816
#> 
#> Average conditional effects (mean; species in columns):
#>  variable     sp1     sp2     sp3    sp4     sp5      sp6     sp7     sp8
#>      prec -0.0857 -0.0196 0.04632 0.0543 -0.0671 -0.01241 -0.1151  0.0760
#>      temp  0.0149  0.0476 0.00532 0.0366  0.0476  0.00353 -0.0787 -0.0217
#>      sp9    sp10    sp11
#>   0.0246 -0.0482 -0.0462
#>  -0.0865 -0.0724  0.0183
# the printed tables are also returned
head(s$importance$growth)
#>   species variable  importance
#> 1       1     temp 0.278906198
#> 2       1     prec 0.064345682
#> 3       2     prec 0.031673597
#> 4       2     temp 0.004184051
#> 5       3     temp 0.520382398
#> 6       3     prec 0.096755242
```
