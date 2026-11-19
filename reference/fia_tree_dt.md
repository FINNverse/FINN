# Oregon FIA calibration data: the trees behind the mortality counts

One row per tree and remeasurement interval for the
[fia_obs_dt](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
training sites (`fia_tree_dt`) and the
[fia_obs_test](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
holdout sites (`fia_tree_test`): every tree alive at the start of an
interval, with its diameter then and its fate over the interval. These
are the closed cohorts whose counts (`n_at_risk`, `n_died`) the site
tables carry, and the `tree_data` of
[`fit(loss_aggregation = c(mortality = "tree"))`](https://finnverse.github.io/FINN/reference/fit.md)
and
[`predictTrees()`](https://finnverse.github.io/FINN/reference/predictTrees.md),
which score the mortality function on every observed tree instead of on
the site's pooled rate (see
[finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md)).
The ids, time axis and species coding are those of the site tables, so
the two join on `siteID`, `year` and `species`. All trees of the FIA
macroplot are kept: the per-tree loss scores fates and needs no area
expansion.

## Usage

``` r
fia_tree_dt

fia_tree_test
```

## Format

A `data.table` with one row per tree and interval (9979 rows and 661
deaths / 10109 rows and 734 deaths) and 8 columns:

- `siteID`, `year`:

  site index and the census that closes the interval, on the axis of
  [fia_obs_dt](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
  (`10`, `20`; `21` or `22`).

- `patchID`:

  patch (FIA subplot) index, as in
  [fia_init_trees](https://finnverse.github.io/FINN/reference/fia_init_trees.md).

- `species`:

  species index `1..11`
  ([fia_species_dt](https://finnverse.github.io/FINN/reference/fia_species_dt.md)),
  the species the tree was recorded as at the start of the interval.

- `dbh`:

  diameter at breast height at the start of the interval, cm.

- `growth`:

  the tree's yearly relative growth rate over the interval,
  `(dbh_end / dbh)^(1 / period_length) - 1`; `NA` for a tree that died.

- `died`:

  `1` if the tree was dead at the end of the interval, `0` if alive.

- `period_length`:

  years from the start to the end of the interval, as in
  [fia_obs_dt](https://finnverse.github.io/FINN/reference/fia_obs_dt.md).

An object of class `data.table` (inherits from `data.frame`) with 10109
rows and 8 columns.

## Source

A seeded subsample of Oregon plots of the US Forest Inventory and
Analysis (FIA) database, with WorldClim bioclimatic variables attached
per plot. The preparation script and the provenance of the source tables
are in the package repository (`dev/make_extdata.R`,
`data-raw/README.md`).

## Details

In a few cells (21 of the training and 22 of the holdout site x year x
species cells) the tree table holds trees that the site table's counts
do not: the site table has no row for a species whose trees had all died
by a census. The tree table keeps them (5 and 7 deaths); everywhere else
the counts agree exactly. Two holdout sites whose first inventory found
only dead trees have no tree at risk and so no row.

## See also

[fia_obs_dt](https://finnverse.github.io/FINN/reference/fia_obs_dt.md),
[`predictTrees()`](https://finnverse.github.io/FINN/reference/predictTrees.md),
[finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md)

Other fia datasets:
[`fia_env_dt`](https://finnverse.github.io/FINN/reference/fia_env_dt.md),
[`fia_init_trees`](https://finnverse.github.io/FINN/reference/fia_init_trees.md),
[`fia_obs_dt`](https://finnverse.github.io/FINN/reference/fia_obs_dt.md),
[`fia_species_dt`](https://finnverse.github.io/FINN/reference/fia_species_dt.md)

## Examples

``` r
library(data.table)
# deaths and trees at risk per interval, and the yearly mortality by size class
fia_tree_dt[, .(trees = .N, deaths = sum(died)), by = year]
#>     year trees deaths
#>    <int> <int>  <num>
#> 1:    10  4935    279
#> 2:    20  5044    382
fia_tree_dt[, .(trees = .N, annual_mortality = round(1 - (1 - mean(died))^(1 / 10), 4)),
            by = .(dbh_class = cut(dbh, c(0, 20, 30, 45, 61, 80, 100, Inf)))][order(dbh_class)]
#>    dbh_class trees annual_mortality
#>       <fctr> <int>            <num>
#> 1:    (0,20]  3010           0.0084
#> 2:   (20,30]  2502           0.0087
#> 3:   (30,45]  1624           0.0058
#> 4:   (45,61]   734           0.0047
#> 5:   (61,80]   982           0.0049
#> 6:  (80,100]   702           0.0030
#> 7: (100,Inf]   425           0.0033
```
