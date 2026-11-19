# Create observation data from trees

Derives stand metrics by site/patch/species/year from tree measurements,
filters sites, harmonizes years, optionally aggregates by site.

## Usage

``` r
makeObsData(
  tree_dt,
  plotsize,
  aggregate_by_site = TRUE,
  minNyears = 2,
  fix_period_length = NULL,
  dbh_growth_thresh = c(-10, 50),
  Npatches = NULL,
  Nspecies = NULL,
  NspeciesQuantile = NULL,
  rates = c("annual", "interval")
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
  fraction of individuals. Used only when `Nspecies` is `NULL`: if both
  are given, `Nspecies` takes precedence.

- rates:

  `"annual"` (default) or `"interval"`: the unit of `growth` and `mort`.
  `"annual"` returns them per year, for a model simulated in yearly
  steps: a tree's growth over an interval of \\p\\ years is
  \\(dbh/dbh\_{before})^{1/p} - 1\\, and the mortality of a cohort is
  \\1 - (1 - n\_{died}/n\_{at\\risk})^{1/p}\\, the constant yearly rates
  that compound to what was observed over the interval. `"interval"`
  returns them over the whole interval (`dbh/dbh_before - 1` and
  `n_died / n_at_risk`), the output of FINN 0.1.0, for a model in which
  one simulation step is one census interval. \\p\\ is
  `year - year_before` in the units of `year`, so `year` must be the
  calendar year of the inventory.

## Value

A list with:

- `obs_dt`: observations at site or patch level, one row per site (and
  patch), inventory `year` and species. `period_length` is the number of
  years since the site's previous inventory (`NA` at its first).
  `growth` is the mean relative diameter growth rate of the surviving
  trees (per year, or over the interval, see `rates`), and `growth_n` is
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
  plus the derived rate `mort` (`NA` where no cohort was at risk): the
  yearly rate `1 - (1 - n_died / n_at_risk)^(1 / period_length)`, or
  `n_died / n_at_risk` with `rates = "interval"`.
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) scores
  `mort` with the binomial loss weighted by `n_at_risk`. `reg` is the
  number of recruits per hectare over the interval.

- `tree_dt`: input trees with added growth fields and species recode:
  `year_before`, `dbh_before`, `period_length` (years since the tree's
  previous record), `rel_growth` (the increment over the interval,
  `dbh/dbh_before - 1`) and `growth` (the tree's growth rate in the unit
  of `rates`, the quantity behind `obs_dt$growth` and the `growth`
  column of
  [`fit(tree_data = )`](https://finnverse.github.io/FINN/reference/fit.md));
  an input column `growth` is replaced.

## Details

FINN simulates the forest in steps of one year. With `rates = "annual"`
the observation table is in the units that
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) compares
with such a simulation:
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) compounds
the simulated yearly growth and mortality of the `period_length` years
before each census into one rate per year with the same functional (the
geometric mean of \\1 + g\\ minus one, one minus the geometric mean
survival), and sums the simulated recruits over them (see
[`finn_data`](https://finnverse.github.io/FINN/reference/finn_data.md)).
`rates = "interval"` belongs with
[`resolveSiteIDs(time_axis = "census")`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md),
where one simulation step is one census interval and nothing is
compounded. The counts (`n_at_risk`, `n_died`, `growth_n`) and `reg`
stay quantities of the whole interval in both modes: the counts weight
the likelihood by the number of trees behind an observation, and `reg`
is compared with the recruits summed over the interval. To use the
table, express `year` as years since the first inventory, give `env` one
row per site and year from 1 to the last census, and build the initial
cohorts from the first inventory;
[`resolveSiteIDs`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
does all three.

## See also

[`resolveSiteIDs`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
for the next step (integer ids, the annual time axis, the initial
cohorts),
[`finn_data`](https://finnverse.github.io/FINN/reference/finn_data.md)
for the tables
[`fit`](https://finnverse.github.io/FINN/reference/fit.md) takes,
[`fia_obs_dt`](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
for a finished example, and
[`vignette("C-Data_preparation", package = "FINN")`](https://finnverse.github.io/FINN/articles/C-Data_preparation.md)

## Examples

``` r
library(data.table)
# A tree list of eight Oregon FIA plots with four patches and three inventories
# each, one row per tree and inventory: siteName, patchName, treeName, year,
# species_name, dbh (cm), status ("new", "alive", "dead") and living.
tree_dt <- fread(system.file("extdata", "example_tree_dt.csv", package = "FINN"))
head(tree_dt[, .(siteName, patchName, treeName, year, species_name, dbh, status)])
#>    siteName  patchName       treeName  year          species_name    dbh status
#>      <char>     <char>         <char> <int>                <char>  <num> <char>
#> 1: 11_66410 11_66410_1 11_66410_1_100  2001 Pseudotsuga menziesii 47.244    new
#> 2: 11_66410 11_66410_1 11_66410_1_100  2011 Pseudotsuga menziesii 57.404  alive
#> 3: 11_66410 11_66410_1 11_66410_1_100  2021 Pseudotsuga menziesii 65.024  alive
#> 4: 11_66410 11_66410_2 11_66410_2_103  2001 Pseudotsuga menziesii 26.924    new
#> 5: 11_66410 11_66410_2 11_66410_2_103  2011 Pseudotsuga menziesii 30.734  alive
#> 6: 11_66410 11_66410_2 11_66410_2_103  2021 Pseudotsuga menziesii 34.036  alive

obs <- makeObsData(
  tree_dt[complete == TRUE],    # fully remeasured plots only
  plotsize          = 0.06,     # ha per patch; recruits are scaled to per hectare
  aggregate_by_site = FALSE,    # keep the patches: resolveSiteIDs() needs them
  fix_period_length = 10,       # drop plots whose remeasurement interval is not 10 years
  NspeciesQuantile  = 0.98)     # species covering 98 % of the stems; the rest is "other"
#> excluded 0 sites due to period length not beeing 10
#> excluded 0 sites due to dbh_growth being not within the range of -10 and 50
#> Nspecies set to 11 to cover 98% of individuals
#> Selected species: Pseudotsuga menziesii, Pinus ponderosa, Tsuga mertensiana, Pinus contorta, Abies grandis, Pinus monticola, Abies concolor, Abies shastensis, other, Juniperus occidentalis, Larix occidentalis
#> Combined 4 unselected species into 'other'
#> actuall number of species including 'others' is 12
#> Filtered to sites with at least 2 years and all patches having the same number of years:
#> Sites with 3 Inventories: 1056
#> Sites have between 4 and 4 patches
#> Kept obs_dt unaggregated by siteName. For aggregated data, set aggregate_by_site = TRUE
head(obs$obs_dt)
#>    siteName  patchName  year   species_name    ba   dbh trees growth growth_n
#>      <char>     <char> <int>         <char> <num> <num> <int>  <num>    <int>
#> 1: 11_66410 11_66410_1  2001 Abies concolor     0    NA     0     NA        0
#> 2: 11_66410 11_66410_1  2011 Abies concolor     0    NA     0     NA        0
#> 3: 11_66410 11_66410_1  2021 Abies concolor     0    NA     0     NA        0
#> 4: 11_66410 11_66410_1  2001  Abies grandis     0    NA     0     NA        0
#> 5: 11_66410 11_66410_1  2011  Abies grandis     0    NA     0     NA        0
#> 6: 11_66410 11_66410_1  2021  Abies grandis     0    NA     0     NA        0
#>     mort n_at_risk n_died   reg period_length
#>    <num>     <int>  <int> <num>         <int>
#> 1:    NA         0      0     0            NA
#> 2:    NA         0      0     0            10
#> 3:    NA         0      0     0            10
#> 4:    NA         0      0     0            NA
#> 5:    NA         0      0     0            10
#> 6:    NA         0      0     0            10

# mortality is a closed cohort: the trees alive at the start of the interval
# and how many of them had died by its end; `mort` is the yearly rate over the
# period_length years of the interval
obs$obs_dt[n_at_risk > 0, .(siteName, patchName, year, period_length, species_name,
                            n_at_risk, n_died, mort)][1:5]
#>    siteName  patchName  year period_length          species_name n_at_risk
#>      <char>     <char> <int>         <int>                <char>     <int>
#> 1: 11_66410 11_66410_1  2011            10 Pseudotsuga menziesii         1
#> 2: 11_66410 11_66410_1  2021            10 Pseudotsuga menziesii         1
#> 3: 11_66410 11_66410_2  2011            10 Pseudotsuga menziesii         8
#> 4: 11_66410 11_66410_2  2021            10 Pseudotsuga menziesii         8
#> 5: 11_66410 11_66410_3  2011            10 Pseudotsuga menziesii         4
#>    n_died       mort
#>     <int>      <num>
#> 1:      0 0.00000000
#> 2:      0 0.00000000
#> 3:      0 0.00000000
#> 4:      1 0.01326438
#> 5:      0 0.00000000
# growth is a yearly rate too, and travels with the number of trees behind it
obs$obs_dt[growth_n > 0, .(siteName, year, species_name, growth, growth_n)][1:5]
#>    siteName  year          species_name      growth growth_n
#>      <char> <int>                <char>       <num>    <int>
#> 1: 11_66410  2011 Pseudotsuga menziesii 0.019669783        1
#> 2: 11_66410  2021 Pseudotsuga menziesii 0.012542247        1
#> 3: 11_66410  2011 Pseudotsuga menziesii 0.011816466        8
#> 4: 11_66410  2021 Pseudotsuga menziesii 0.007917063        7
#> 5: 11_66410  2011 Pseudotsuga menziesii 0.018626114        4

# the same table with the rates over the whole ten-year interval (FINN 0.1.0)
obs10 <- makeObsData(tree_dt[complete == TRUE], plotsize = 0.06, aggregate_by_site = FALSE,
                     fix_period_length = 10, NspeciesQuantile = 0.98, rates = "interval")
#> excluded 0 sites due to period length not beeing 10
#> excluded 0 sites due to dbh_growth being not within the range of -10 and 50
#> Nspecies set to 11 to cover 98% of individuals
#> Selected species: Pseudotsuga menziesii, Pinus ponderosa, Tsuga mertensiana, Pinus contorta, Abies grandis, Pinus monticola, Abies concolor, Abies shastensis, other, Juniperus occidentalis, Larix occidentalis
#> Combined 4 unselected species into 'other'
#> actuall number of species including 'others' is 12
#> Filtered to sites with at least 2 years and all patches having the same number of years:
#> Sites with 3 Inventories: 1056
#> Sites have between 4 and 4 patches
#> Kept obs_dt unaggregated by siteName. For aggregated data, set aggregate_by_site = TRUE
obs10$obs_dt[growth_n > 0, .(siteName, year, species_name, growth, growth_n)][1:5]
#>    siteName  year          species_name     growth growth_n
#>      <char> <int>                <char>      <num>    <int>
#> 1: 11_66410  2011 Pseudotsuga menziesii 0.21505376        1
#> 2: 11_66410  2021 Pseudotsuga menziesii 0.13274336        1
#> 3: 11_66410  2011 Pseudotsuga menziesii 0.12644955        8
#> 4: 11_66410  2021 Pseudotsuga menziesii 0.08319808        7
#> 5: 11_66410  2011 Pseudotsuga menziesii 0.20321636        4
```
