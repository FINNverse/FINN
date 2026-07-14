# Create observation data from trees

Derives stand metrics by site/patch/species/year from tree measurements,
filters sites, harmonizes years, optionally aggregates by site.

## Usage

``` r
makeObsData(
  tree_dt,
  plotsize,
  aggregate_by_site = T,
  minNyears = 2,
  fix_period_length = NULL,
  dbh_growth_thresh = c(-10, 50),
  Npatches = NULL,
  Nspecies = NULL,
  NspeciesQuantile = NULL
)
```

## Arguments

- tree_dt:

  data.table of tree records with siteName, patchName, year, treeName,
  species_name, dbh, status, living, and optional `mort`.

- plotsize:

  numeric plot area used to scale recruitment.

- aggregate_by_site:

  logical. Aggregate patches to site level. Default TRUE.

- minNyears:

  integer or NULL. Keep sites with at least this many years and equal
  counts across patches. Default 2.

- fix_period_length:

  integer or NULL. If set, drop sites whose inventory interval differs.
  Default NULL.

- dbh_growth_thresh:

  length-2 numeric or NULL. Drop sites where any tree `dbh` change falls
  outside this range. Default `c(-10, 50)`.

- Npatches:

  integer or NULL. If set, keep only sites with exactly this many
  patches.

- Nspecies:

  integer or NULL. Cap species to the top N (others merged to "other").

- NspeciesQuantile:

  numeric in (0,1\] or NULL. Choose the smallest N covering this
  fraction of individuals; overrides `Nspecies` if supplied.

## Value

A list with:

- `obs_dt`: observations at site or patch level.

- `tree_dt`: input trees with added growth fields and species recode.
