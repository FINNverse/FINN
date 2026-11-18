# Make initial cohorts for FINN

Builds the initial state of a fit or a simulation, the `init_cohort` of
[`fit`](https://finnverse.github.io/FINN/reference/fit.md) and
[`simulateForest`](https://finnverse.github.io/FINN/reference/predict.finn_class.md),
from the trees of the first inventory: every tree becomes a cohort of
its own, or, with `dbh_binsize`, the trees of a patch and species are
pooled into diameter classes. The function prepares `obs_df` in the
format expected by
[`CohortMat`](https://finnverse.github.io/FINN/reference/CohortMat.md)
and calls it. The long-table schema of `obs_df` is:

- **siteID** integer index of sites \\1..S\\.

- **patchID** integer index of patches within site \\1..P_s\\.

- **species** integer species code \\1..sp\\.

- **dbh** numeric DBH in cm (either exact or binned midpoints).

- **trees** integer count of trees in the cohort.

## Usage

``` r
makeInitCohorts(
  init_trees,
  dbh_binsize = NULL,
  min_dbh = NULL,
  Nspecies,
  treeID_table = FALSE,
  singleCohortTreeNames = NULL
)
```

## Arguments

- init_trees:

  data.table of initial trees with columns `siteID`, `patchID`,
  `species`, `dbh` (cm), `treeName` and `trees` (the number of trees the
  row stands for, usually 1). Without `dbh_binsize` a column `cohortID`
  is added to the table by reference.

- dbh_binsize:

  numeric or NULL. Bin width in cm. If NULL, keep exact DBH.

- min_dbh:

  numeric or NULL. Lower bound for binning. If NULL, uses min DBH.

- Nspecies:

  integer. Number of species levels passed to `sp`.

- treeID_table:

  logical. If TRUE, also return the cohort table used.

- singleCohortTreeNames:

  character vector or NULL. Tree names excluded from binning and kept as
  single-tree cohorts.

## Value

If `treeID_table = FALSE`, a
[`FINN::CohortMat`](https://finnverse.github.io/FINN/reference/CohortMat.md)
object. If `TRUE`, a list with:

- `initCohort`: the `CohortMat` object.

- `init_trees`: the `obs_df` passed to
  [`CohortMat()`](https://finnverse.github.io/FINN/reference/CohortMat.md),
  with the `cohortID` of every row.

## Details

Aggregates initial trees into cohorts by DBH bins or exact DBH and
constructs a
[`FINN::CohortMat`](https://finnverse.github.io/FINN/reference/CohortMat.md)
object.

Internally calls:


    FINN::CohortMat(
      obs_df = <data.frame with columns siteID, patchID, cohortID, species, dbh, trees>,
      dbh    = NULL,
      trees  = NULL,
      species= NULL,
      dims   = c(S, P, K),   # inferred from obs_df
      sp     = Nspecies,     # passed from argument
      device = "cpu"
    )

Key fields of the resulting `CohortMat` (a `torch` `nn_module`):

- `dbh`, `trees`, `species`: tensors per cohort.

- `dims`: integer vector `c(S, P, K)` for sites, patches, cohorts.

- `sp`: integer number of species.

- `device_r`: the device given at construction ("cpu").

- `dbh_r`, `trees_r`, `species_r`: R arrays.

- `obsDF2arrays`: method converting `obs_df` into arrays.

## See also

[`CohortMat`](https://finnverse.github.io/FINN/reference/CohortMat.md),
[`resolveSiteIDs`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md),
which calls this function for the first inventory,
[`fit`](https://finnverse.github.io/FINN/reference/fit.md) and
[`simulateForest`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
(`init_cohort`),
[`finn_data`](https://finnverse.github.io/FINN/reference/finn_data.md)

## Examples

``` r
library(data.table)
# one cohort per observed tree, which is what `loss_resolution = "tree"` needs
init <- makeInitCohorts(fia_init_trees[siteID %in% 1:3],
                        Nspecies = max(fia_obs_dt$species))
dim(init$dbh)   # sites x patches x cohorts
#> [1]  3  4 16

sum(init$trees_r)   # the living trees of the three sites
#> [1] 61

# 5 cm diameter bins instead: fewer cohorts, so a much cheaper simulation.
# Binning adds up `trees`; rows without a count (the records of dead trees
# in this table, trees = NA) carry no stems, so leave them out first, or a
# bin holding only such rows becomes an empty cohort
living <- fia_init_trees[siteID %in% 1:3 & !is.na(trees)]
binned <- makeInitCohorts(living, dbh_binsize = 5,
                          Nspecies = max(fia_obs_dt$species))
dim(binned$dbh)
#> [1]  3  4 13
sum(binned$trees_r)   # the same trees in fewer cohorts
#> [1] 61

# the cohort table behind the tensors, with the cohort of every tree
tab <- makeInitCohorts(living, Nspecies = max(fia_obs_dt$species),
                       treeID_table = TRUE)
head(tab$init_trees[, .(siteID, patchID, cohortID, species, dbh, trees, treeName)])
#>    siteID patchID cohortID species    dbh trees      treeName
#>     <int>   <int>    <int>   <int>  <num> <int>        <char>
#> 1:      1       2        1       1 80.264     1 1_50638_2_121
#> 2:      1       2        2       5 13.716     1 1_50638_2_111
#> 3:      1       2        3       5 14.986     1 1_50638_2_112
#> 4:      1       2        4       5 20.574     1 1_50638_2_114
#> 5:      1       2        5      11 81.026     1 1_50638_2_115
#> 6:      1       2        6      11 68.326     1 1_50638_2_120
```
