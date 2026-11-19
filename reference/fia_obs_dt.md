# Oregon FIA calibration data: training sites

A ready-to-fit sample of the US Forest Inventory and Analysis (FIA)
network in Oregon: 200 training sites (`fia_*_dt` / `fia_init_trees`)
and 200 disjoint holdout sites (`fia_*_test`) that are never seen during
fitting, so a model fit on the training tables can be scored genuinely
out of sample. The tables are already in FINN's input format – IDs
resolved, species re-indexed `1..K` – so they go straight into
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) /
[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
without
[`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
or
[`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md).
The `D-Fit_to_FIA` and `E-Mortality` vignettes are the worked examples.

## Usage

``` r
fia_obs_dt

fia_obs_test
```

## Format

`fia_obs_dt` and `fia_obs_test`: one row per site x year x species (4400
rows, 14 columns). The year-0 initial state is in `fia_init_trees` /
`fia_init_test`.

- `siteID`, `year`:

  integer site index (`1..200`) and the years since the site's first
  inventory (`10`, `20`; `21` or `22` for one site each).

- `period_length`:

  years since the previous inventory (`10`; `11` or `12` for the late
  inventory): the window over which the simulation is aggregated before
  it meets the observation.

- `species`, `species_name`:

  species index `1..11` and its name; `"other"` pools the species
  outside the 98 % most abundant stems.

- `ba`:

  living basal area, m2 – the mean over the site's patches (0.06 ha
  each;
  [`makeObsData(aggregate_by_site = TRUE)`](https://finnverse.github.io/FINN/reference/makeObsData.md)).

- `trees`:

  living stems, mean per patch (hence fractional).

- `dbh`:

  mean diameter at breast height of living trees, cm (`NA` when there
  are none).

- `growth`, `growth_n`:

  mean annual relative dbh growth rate of the survivors,
  `(dbh / dbh_before)^(1 / period_length) - 1` per tree (tree-weighted
  over patches), and the number of trees behind it.

- `n_at_risk`, `n_died`, `mort`:

  trees alive at the start of the interval and how many of them had died
  by its end, summed over patches, and the annual mortality rate
  `1 - (1 - n_died / n_at_risk)^(1 / period_length)`. The binomial
  mortality loss scores `mort` weighted by `n_at_risk`.

- `reg`:

  recruits over the interval, trees per hectare (mean over patches).

An object of class `data.table` (inherits from `data.frame`) with 4400
rows and 14 columns.

## Source

A seeded subsample of Oregon plots of the US Forest Inventory and
Analysis (FIA) database, with WorldClim bioclimatic variables attached
per plot. The preparation script and the provenance of the source tables
are in the package repository (`dev/make_extdata.R`,
`data-raw/README.md`).

## Details

All nine objects are `data.table`s. Attach `data.table` to use its `[`
semantics on them, and
[`data.table::copy()`](https://rdrr.io/pkg/data.table/man/copy.html) one
before modifying it by reference (`:=`), as with any package dataset.

The tables are on an **annual** time axis: one simulation step is one
year. The first inventory (2001-02) is the initial state, year 0; the
two remeasurements are observations at the number of years since it, 10
and 20 (21 for one training site and 22 for one holdout site, whose
third inventory came late), each with `period_length`, the years since
the previous inventory. Growth and mortality are yearly rates;
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) compounds
the simulated yearly rates of the `period_length` years before each
inventory into the same per-year rate before comparing, and compares
`reg` with the recruits simulated over those years. The trees behind the
mortality counts are in
[fia_tree_dt](https://finnverse.github.io/FINN/reference/fia_tree_dt.md).
Before FINN 0.2.0 the tables had one step per ten-year interval (`year`
1 and 2, rates over the interval);
[`makeObsData(rates = "interval")`](https://finnverse.github.io/FINN/reference/makeObsData.md)
and
[`resolveSiteIDs(time_axis = "census")`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
still build such tables.

## See also

[fia_env_dt](https://finnverse.github.io/FINN/reference/fia_env_dt.md),
[fia_init_trees](https://finnverse.github.io/FINN/reference/fia_init_trees.md),
[fia_tree_dt](https://finnverse.github.io/FINN/reference/fia_tree_dt.md),
[fia_species_dt](https://finnverse.github.io/FINN/reference/fia_species_dt.md)

Other fia datasets:
[`fia_env_dt`](https://finnverse.github.io/FINN/reference/fia_env_dt.md),
[`fia_init_trees`](https://finnverse.github.io/FINN/reference/fia_init_trees.md),
[`fia_species_dt`](https://finnverse.github.io/FINN/reference/fia_species_dt.md),
[`fia_tree_dt`](https://finnverse.github.io/FINN/reference/fia_tree_dt.md)

## Examples

``` r
library(data.table)
fia_obs_dt[year == 10, .(sites = uniqueN(siteID), ba = mean(ba)), by = species_name]
#>                species_name sites          ba
#>                      <char> <int>       <num>
#>  1:   Pseudotsuga menziesii   200 0.726985711
#>  2:         Pinus ponderosa   200 0.184966585
#>  3:          Pinus contorta   200 0.024908897
#>  4:          Abies concolor   200 0.055690344
#>  5:           Abies grandis   200 0.054632263
#>  6:      Tsuga heterophylla   200 0.056597997
#>  7:       Tsuga mertensiana   200 0.015395718
#>  8: Lithocarpus densiflorus   200 0.010621666
#>  9:             Alnus rubra   200 0.013027215
#> 10:          Abies amabilis   200 0.006718681
#> 11:                   other   200 0.123505996
# the observation windows: inventory year x interval length, sites
unique(fia_obs_dt[, .(siteID, year, period_length)])[, .N, by = .(year, period_length)]
#>     year period_length     N
#>    <int>         <int> <int>
#> 1:    10            10   200
#> 2:    20            10   199
#> 3:    21            11     1
```
