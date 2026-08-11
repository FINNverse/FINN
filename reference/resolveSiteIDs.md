# Resolve site, patch, and year indices for FINN inputs

Joins tree, environment, and observation tables; assigns integer indices
for site, patch, period; standardizes species coding; and returns
aligned data plus optional initial cohorts.

## Usage

``` r
resolveSiteIDs(tree_dt, env_dt, obs_dt, createInitCohorts = TRUE)
```

## Arguments

- tree_dt:

  data.table with tree-level data including siteName, patchName, year,
  species_name, dbh, status, living, and optional trees.

- env_dt:

  data.table with environment data including siteName and year.

- obs_dt:

  data.table with observations including siteName, patchName, year,
  species_name, and stand metrics (ba, dbh, trees, growth, mort,
  n_at_risk, n_died, reg), as returned by `makeObsData`.

- createInitCohorts:

  logical. If TRUE, build FINN initial cohorts from trees with
  `living == TRUE` and `year == 1`. Default TRUE.

## Value

A list with:

- `siteID_dt`: site/patch/year index map.

- `tree_dt`: tree table with indices and standardized species.

- `env_dt`: environment table with indices.

- `obs_dt`: site-aggregated observations by species.

- `obs_dt_patches`: patch-level observations by species.

- `species_dt`: species lookup (`species`, `species_name`).

- `initCohorts` (optional): FINN `CohortMat` object.

- `initCohort_dt` (optional): trees used to build cohorts.
