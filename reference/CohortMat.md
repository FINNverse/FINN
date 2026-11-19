# Cohort Matrix Class

A `torch` `nn_module` holding the cohorts of every site and patch as
three aligned tensors – `dbh`, `trees` and `species` – each shaped
`[sites, patches, cohorts]`. It is FINN's initial state: build it from
an inventory with
[`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md)
or from a cohort table with `obs_df`, and pass it to
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md),
[`predict.finn_class()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
or
[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
as `init_cohort`. The three tensors are registered as buffers, so a
`CohortMat` moves between devices with the model and carries no
gradient. `asDF()` converts it back to a cohort table (a `data.frame`
with `siteID`, `patchID`, `cohortID`, `species`, `dbh`, `trees`).

Build the cohort tensors, either from a cohort table (`obs_df`) or from
arrays of `dbh`, `trees` and `species`

Store a plain R copy of the cohort arrays (used by `restore()`)

Rebuild the tensors from the stored R copy when they are no longer
valid, as after [`readRDS()`](https://rdrr.io/r/base/readRDS.html) of a
saved `CohortMat` (or of the list
[`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
returns) in a new session

Convert the cohort tensors back to a cohort table

## Usage

``` r
CohortMat(
  obs_df = NULL,
  dbh = NULL,
  trees = NULL,
  species = NULL,
  dims = c(50, 30, 10),
  sp = 10,
  device = "cpu"
)
```

## Arguments

- obs_df:

  A data frame containing columns "siteID", "patchID", "species", "dbh",
  and "trees". If provided, it will be used to initialize the tensors.

- dbh:

  A tensor or array representing the diameter at breast height. Defaults
  to `self$dbh`.

- trees:

  A tensor or array representing the number of trees. Defaults to
  `self$trees`.

- species:

  A tensor or array representing the species. Defaults to
  `self$species`.

- dims:

  A numeric vector representing the dimensions of the arrays (sites,
  patches, cohorts). Defaults to `self$dims`.

- sp:

  An integer representing the number of species. Defaults to `self$sp`.

- device:

  A character string specifying the device to use ('cpu' or 'cuda').
  Defaults to `self$device`.

## Value

A `CohortMat`, a `torch` `nn_module` whose `dbh`, `trees` and `species`
buffers hold the cohorts as `[sites, patches, cohorts]` tensors.

The `CohortMat`, invisibly.

A `data.frame` with one row per cohort and the columns `siteID`,
`patchID`, `cohortID`, `species`, `dbh` and `trees`.

## Fields

- `dbh`:

  A tensor representing the diameter at breast height for each cohort.

- `trees`:

  A tensor representing the number of trees in each cohort.

- `species`:

  A tensor representing the species of each cohort.

- `dims`:

  A vector representing the dimensions of the arrays (sites, patches,
  cohorts).

- `sp`:

  An integer representing the number of species.

- `dbh_r`:

  A numeric array representing the diameter at breast height in R array
  format.

- `trees_r`:

  A numeric array representing the number of trees in R array format.

- `species_r`:

  An integer array representing the species in R array format.

- `device_r`:

  A character string specifying the device in R format ('cpu' or
  'cuda').

- `obsDF2arrays`:

  Transforms data.table into array

- `site_ids`:

  The user `siteID`s in the order of the site dimension when the cohorts
  were built from `obs_df` (the lookup table from array position to user
  id); `NULL` when built from arrays, where the site dimension is `1:N`.

## See also

[`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md)
to build one from a tree list,
[`obsDF2arrays()`](https://finnverse.github.io/FINN/reference/obsDF2arrays.md)
and
[`array2obsDF()`](https://finnverse.github.io/FINN/reference/array2obsDF.md)
for the conversion between table and arrays,
[finn_data](https://finnverse.github.io/FINN/reference/finn_data.md) for
the role of `init_cohort`

## Examples

``` r
# from a cohort table: one site, two patches, three cohorts
cohorts <- data.frame(siteID = 1, patchID = c(1, 1, 2), cohortID = c(1, 2, 1),
                      species = c(1, 2, 1), dbh = c(35, 12, 20), trees = c(3, 10, 6))
init <- CohortMat(obs_df = cohorts, sp = 2)
init$dbh$shape   # sites x patches x cohorts
#> [1] 1 2 2
init$trees_r     # the same content as an R array; empty slots hold 0
#> , , 1
#> 
#>      [,1] [,2]
#> [1,]    3    6
#> 
#> , , 2
#> 
#>      [,1] [,2]
#> [1,]   10    0
#> 
init$asDF()      # and back to a table
#>   siteID patchID cohortID species dbh trees
#> 1      1       1        1       1  35     3
#> 2      1       1        2       2  12    10
#> 3      1       2        1       1  20     6

# from arrays: every patch of two sites starts with one cohort of twenty
# 10-cm trees of species 1
dims <- c(2, 4, 1)
init2 <- CohortMat(dbh = array(10, dims), trees = array(20, dims),
                   species = array(1L, dims), dims = dims, sp = 2)
head(init2$asDF(), 3)
#>   siteID patchID cohortID species dbh trees
#> 1      1       1        1       1  10    20
#> 2      1       2        1       1  10    20
#> 3      1       3        1       1  10    20

# trees of the bundled Oregon FIA data as an initial state (?fia_init_trees)
first3 <- fia_init_trees[fia_init_trees$siteID %in% 1:3, ]
init3 <- makeInitCohorts(first3, Nspecies = max(fia_init_trees$species))
init3$dims
#> [1]  3  4 16
```
