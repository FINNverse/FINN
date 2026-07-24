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
  species_name, dbh, `status` (one of "alive", "new", "dead") and
  `living`. Mortality is derived from `status`; a `mort` column, if
  present, is ignored.

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

- `obs_dt`: observations at site or patch level. `growth` is the mean
  relative diameter increment (`dbh/dbh_before - 1`), and `growth_n` is
  how many trees that mean rests on. Both matter: means are aggregated
  across patches TREE-WEIGHTED (by `growth_n`), because that is the
  quantity the model predicts - it simulates `sum(g*trees)/sum(trees)`
  over the whole site, so an unweighted mean of patch means would be a
  different number. `growth_n` also travels with the response so `fit`
  can weight the likelihood by it: 34% of patch-level growth
  observations rest on a single tree while others average 30+, and the
  variance of a mean is \\\sigma^2/n\\. Mortality comes back as a
  closed-cohort pair of counts, `n_at_risk` (trees alive at the start of
  the interval) and `n_died` (how many of them were dead at the end),
  plus the derived rate `mort = n_died / n_at_risk` (`NA` where no
  cohort was at risk). The counts are the binomial response; pass them
  to `fit` with `mortality = "binomial"`.

- `tree_dt`: input trees with added growth fields and species recode.
