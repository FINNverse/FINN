# Oregon FIA calibration data: site climate

Climate per site and year for the
[fia_obs_dt](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
training sites (`fia_env_dt`) and the
[fia_obs_test](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
holdout sites (`fia_env_test`), in **natural units** (not standardised –
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) z-scales
internally with `env_autoscale`). One row per simulated year `1..22`:
the row of year `k` drives the year that ends `k` years after the first
inventory, so every site x year in the observation tables is covered,
and the years after a site's last inventory are simulated but not
scored. The climate is the plot's WorldClim normal, the same in every
year.

## Usage

``` r
fia_env_dt

fia_env_test
```

## Format

A `data.table` with 4400 rows and 8 columns: 200 sites x 22 years, one
row per site x year.

- `siteID`, `year`:

  site index and simulated year (`1..22`).

- `temp`, `tempmax`, `tempmin`:

  annual mean, warmest-month maximum and coldest-month minimum
  temperature, degrees C (bioclim 1, 5, 6).

- `prec`, `precwarmq`:

  annual and warmest-quarter precipitation, mm (bioclim 12, 18).

- `precseas`:

  precipitation seasonality, coefficient of variation (bioclim 15).

An object of class `data.table` (inherits from `data.frame`) with 4400
rows and 8 columns.

## Source

A seeded subsample of Oregon plots of the US Forest Inventory and
Analysis (FIA) database, with WorldClim bioclimatic variables attached
per plot. The preparation script and the provenance of the source tables
are in the package repository (`dev/make_extdata.R`,
`data-raw/README.md`).

## See also

Other fia datasets:
[`fia_init_trees`](https://finnverse.github.io/FINN/reference/fia_init_trees.md),
[`fia_obs_dt`](https://finnverse.github.io/FINN/reference/fia_obs_dt.md),
[`fia_species_dt`](https://finnverse.github.io/FINN/reference/fia_species_dt.md),
[`fia_tree_dt`](https://finnverse.github.io/FINN/reference/fia_tree_dt.md)

## Examples

``` r
library(data.table)
head(fia_env_dt)
#>    siteID  year     temp tempmax tempmin  prec precseas precwarmq
#>     <int> <int>    <num>   <num>   <num> <int>    <num>     <int>
#> 1:      1     1 3.958333    24.7   -10.4   547 30.04612       101
#> 2:      1     2 3.958333    24.7   -10.4   547 30.04612       101
#> 3:      1     3 3.958333    24.7   -10.4   547 30.04612       101
#> 4:      1     4 3.958333    24.7   -10.4   547 30.04612       101
#> 5:      1     5 3.958333    24.7   -10.4   547 30.04612       101
#> 6:      1     6 3.958333    24.7   -10.4   547 30.04612       101
# one row per site and simulated year; the climate of a site is its normal,
# repeated for every year
fia_env_dt[, .(sites = uniqueN(siteID), first_year = min(year), last_year = max(year))]
#>    sites first_year last_year
#>    <int>      <int>     <int>
#> 1:   200          1        22
# the range of the two predictors the package examples use, over the sites
fia_env_dt[year == 1, .(temp_min = min(temp), temp_max = max(temp),
                        prec_min = min(prec), prec_max = max(prec))]
#>     temp_min temp_max prec_min prec_max
#>        <num>    <num>    <int>    <int>
#> 1: 0.4041666 12.12917      249     2786
```
