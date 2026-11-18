# Predict growth and mortality for individual observed trees

Predict growth and mortality for individual observed trees

## Usage

``` r
predictTrees(
  model,
  tree_data,
  env,
  patches = 100L,
  patch_size = NULL,
  device = c("cpu", "gpu")
)
```

## Arguments

- model:

  (`finn_class`)  
  A (fitted) model.

- tree_data:

  (`data.table|data.frame`)  
  Observed trees: `siteID`, `year`, `patchID`, `species`, `dbh` (at the
  interval start) and `period_length` (the interval, in timesteps).

- env:

  (`data.table|data.frame`)  
  Environment, as for
  [fit](https://finnverse.github.io/FINN/reference/fit.md).

- patches:

  (`integer(1)`)  
  Not used: the patch dimension of the observed-tree tensors is the
  number of distinct `patchID`s in `tree_data` (it used to be `patches`,
  and the extra patches were empty). Kept so existing calls run
  unchanged.

- patch_size:

  (`numeric(1)`)  
  Area of one patch in ha. Defaults to the one the model was last fitted
  or simulated with (`model$patch_size_ha`); pass a value only to
  override it. A model that was never fitted or simulated has none, so
  give it here. The value is stored on the model.

- device:

  (`character(1)`)  
  `"cpu"` or `"gpu"`; `"cuda"` is accepted as a synonym for `"gpu"`.

## Value

One row per row of `tree_data`, in the same order (`.row` records the
original position), with `growth_pred` (the relative growth rate per
year over the interval, compounded as described above) and `mort_pred`
(probability of dying during the interval). Bind these onto `tree_data`
directly rather than joining on siteID/year/patchID, which would be
many-to-many.

## Details

The tree-level counterpart of
[predict.finn_class](https://finnverse.github.io/FINN/reference/predict.finn_class.md).
Where [`predict()`](https://rdrr.io/r/stats/predict.html) simulates a
stand forward and reports site-level states, this evaluates the growth
and mortality kernels on trees you observed: light is computed from the
observed stand, the kernels run on those trees at their observed
diameters, and the annual rates are accumulated over each site's
remeasurement interval — growth as the per-year rate of the compounded
annual rates (the geometric mean of `1 + g`, as
[`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
derives the observed rate), mortality as the complement of the product
of annual survivals.

It is the same code path
`fit(..., loss_aggregation = c(growth = "tree", mortality = "tree"))`
scores, so a validation report cannot drift from the objective it is
validating.

## See also

[`predict.finn_class()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
/
[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
for the simulated stand,
[finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md) for
the per-tree losses that share this code,
[fia_tree_dt](https://finnverse.github.io/FINN/reference/fia_tree_dt.md)
for a table of observed trees with their fates

## Examples

``` r
library(data.table)
# the growth and mortality kernels applied to the trees as observed, so the
# error of the two functions alone, with no simulated stand in between
sites <- 1:3
Nsp   <- max(fia_obs_dt$species)
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
FINN.seed(1)
m <- finn(
  N_species            = Nsp, recruits_dbh = 12.9,
  growth_process       = createProcess(~ temp + prec, FINN::growth),
  mortality_process    = createProcess(~ temp + prec, FINN::mortality),
  regeneration_process = createProcess(~ temp + prec, FINN::regeneration))
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
    plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.768)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%

# `fia_init_trees` is a pure initial state, so name the census the trees are
# scored against and the interval, in timesteps (years), that reaches it
trees <- fia_init_trees[siteID %in% sites,
                        .(siteID, patchID, species, dbh, year = 10L,
                          period_length = 10)]
pred  <- predictTrees(m, tree_data = trees,
                      env = fia_env_dt[siteID %in% sites], patch_size = 0.06)
#> [tree_data] 87 trees, up to 16 per patch and observation year
head(pred)
#>     .row siteID  year patchID growth_pred mort_pred
#>    <int>  <int> <int>   <int>       <num>     <num>
#> 1:     1      1    10       1 0.006177030 0.1652927
#> 2:     2      1    10       1 0.000000000 0.1834491
#> 3:     3      1    10       1 0.001696481 0.1784589
#> 4:     4      1    10       1 0.007360663 0.1722142
#> 5:     5      1    10       1 0.000000000 0.1835063
#> 6:     6      1    10       2 0.016580803 0.1630244

# the bundled tree table carries the observed fates, so prediction and
# observation can be compared tree by tree: here the predicted probability
# of dying within the interval for the trees that died and those that lived
obs_trees <- fia_tree_dt[siteID %in% sites]
p <- predictTrees(m, tree_data = obs_trees, env = fia_env_dt[siteID %in% sites],
                  patch_size = 0.06)
#> [tree_data] 113 trees, up to 14 per patch and observation year
tapply(p$mort_pred, obs_trees$died, mean)
#>         0         1 
#> 0.1686863 0.1718107 
```
