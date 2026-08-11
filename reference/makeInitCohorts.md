# Make initial cohorts for FINN

This function prepares `obs_df` in the exact format expected by
`FINN::CohortMat$new()` and calls it. The required long-table schema is:

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
  `species`, `dbh`, `treeName`; optional `trees` for pre-counted
  individuals.

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

- `init_trees`: the `obs_df` passed to `CohortMat$new()`.

## Details

Aggregates initial trees into cohorts by DBH bins or exact DBH and
constructs a
[`FINN::CohortMat`](https://finnverse.github.io/FINN/reference/CohortMat.md)
object.

Internally calls:


    FINN::CohortMat$new(
      obs_df = <data.frame with columns siteID, patchID, species, dbh, trees>,
      dbh    = NULL,
      trees  = NULL,
      species= NULL,
      dims   = c(S, P, K),   # inferred from obs_df
      sp     = Nspecies,     # passed from argument
      device = "cpu"
    )

Key fields of the resulting R6 object:

- `dbh`, `trees`, `species`: tensors per cohort.

- `dims`: integer vector `c(S, P, K)` for sites, patches, cohorts.

- `sp`: integer number of species.

- `device`: "cpu" or "cuda".

- `dbh_r`, `trees_r`, `species_r`: R arrays.

- `obsDF2arrays`: method converting `obs_df` into arrays.

## See also

[`CohortMat`](https://finnverse.github.io/FINN/reference/CohortMat.md)
