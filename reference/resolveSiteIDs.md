# Resolve site, patch, and year indices for FINN inputs

Joins tree, environment, and observation tables; assigns integer indices
for site, patch and year; standardizes species coding; and returns
aligned data plus optional initial cohorts.

## Usage

``` r
resolveSiteIDs(
  tree_dt,
  env_dt,
  obs_dt,
  createInitCohorts = TRUE,
  time_axis = c("annual", "census")
)
```

## Arguments

- tree_dt:

  data.table with tree-level data including siteName, patchName, year,
  species_name, dbh, status, living, and optional trees.

- env_dt:

  data.table with environment data including siteName and year (calendar
  year; see Details for the years needed).

- obs_dt:

  data.table with observations including siteName, patchName, year,
  species_name, and stand metrics (ba, dbh, trees, growth, mort,
  n_at_risk, n_died, reg, period_length), as returned by `makeObsData`.

- createInitCohorts:

  logical. If TRUE, build FINN initial cohorts from the living trees of
  the first inventory (`year == 0`) on either time axis. Default TRUE.

- time_axis:

  `"annual"` (default) or `"census"`. `"census"` is the axis of FINN
  0.1.0: `year` is the index of the inventory (0, 1, 2, ...), one
  simulation step is one census interval (`period_length = 1`), `env_dt`
  keeps the inventory years, and the observations of the first inventory
  stay in `obs_dt`. Use it with
  [`makeObsData(rates = "interval")`](https://finnverse.github.io/FINN/reference/makeObsData.md).
  Unlike 0.1.0, which built the initial cohorts from the second
  inventory, the census axis now starts from the first one, so 0.1.0
  results are not reproduced exactly.

## Value

A list with:

- `siteID_dt`: site/patch/year index map (`OrigYear` is the calendar
  year).

- `tree_dt`: tree table with indices and standardized species.

- `env_dt`: environment table with indices, the years of the simulation:
  `siteID`, `year` and the predictors, what
  [`fit(env = )`](https://finnverse.github.io/FINN/reference/fit.md)
  takes (`siteID_dt` maps the ids back to the site names).

- `obs_dt`: site-aggregated observations by species, with
  `period_length`.

- `obs_dt_patches`: patch-level observations by species.

- `species_dt`: species lookup (`species`, `species_name`).

- `initCohorts` (optional): FINN `CohortMat` object.

- `initCohort_dt` (optional): trees used to build cohorts.

The list can be stored with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and read back with
[`readRDS()`](https://rdrr.io/r/base/readRDS.html): `initCohorts` keeps
a plain R copy of its arrays and rebuilds its torch tensors, which
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) cannot store, when
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md),
[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
or
[`predict.finn_class()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
use it (or by `x$initCohorts$restore()`, see
[`CohortMat`](https://finnverse.github.io/FINN/reference/CohortMat.md)).
Cohorts saved with FINN before 0.2.0 carry no such copy; rebuild them
from the table with
`x$initCohorts <- makeInitCohorts(x$initCohort_dt, Nspecies = nrow(x$species_dt))`.

## Details

With `time_axis = "annual"` (the default) one simulation step is one
year: `year` becomes the number of years since the first inventory of
the site's patch (0 at the first inventory, e.g. 10 and 20 at the next
two), the first inventory is the initial state (the initial cohorts, not
an observation), the observations carry `period_length` in years, and
`env_dt` is returned with one row per site and year `1..Tmax`, `Tmax`
being the longest time from a first to a last inventory. The environment
row of year \\k\\ drives the simulated year that ends at year \\k\\, so
an inventory \\k\\ years after the first is compared with the stand
after \\k\\ simulated years. `env_dt` therefore needs one row per site
and calendar year from the year after the site's first inventory to
`Tmax` years after it; a climate that does not change is repeated for
every year.

## See also

[`makeObsData`](https://finnverse.github.io/FINN/reference/makeObsData.md)
for the step before,
[`makeInitCohorts`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md)
for the initial cohorts on their own,
[`fit`](https://finnverse.github.io/FINN/reference/fit.md) and
[`finn_data`](https://finnverse.github.io/FINN/reference/finn_data.md)
for where the tables go

## Examples

``` r
library(data.table)
# eight Oregon FIA plots, three inventories each; see the data-preparation vignette
tree_dt <- fread(system.file("extdata", "example_tree_dt.csv", package = "FINN"))
env_dt  <- fread(system.file("extdata", "example_env_dt.csv",  package = "FINN"))
obs <- makeObsData(tree_dt[complete == TRUE], plotsize = 0.06, aggregate_by_site = FALSE,
                   fix_period_length = 10, NspeciesQuantile = 0.98)
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

# Names become integer ids shared by every table, calendar years become years
# since the first inventory: 0 (the initial state), 10 and 20. Pass the tree
# list makeObsData() returned: it carries the same species lumping as obs_dt.
# createInitCohorts = FALSE leaves out the starting cohorts, which need torch.
inputs <- resolveSiteIDs(tree_dt = obs$tree_dt, obs_dt = obs$obs_dt, env_dt = env_dt,
                         createInitCohorts = FALSE)
inputs$species_dt
#>     species           species_name
#>       <int>                 <char>
#>  1:       9  Pseudotsuga menziesii
#>  2:      11                  other
#>  3:       8        Pinus ponderosa
#>  4:       5     Larix occidentalis
#>  5:       2          Abies grandis
#>  6:       4 Juniperus occidentalis
#>  7:       7        Pinus monticola
#>  8:       6         Pinus contorta
#>  9:       1         Abies concolor
#> 10:      10      Tsuga mertensiana
#> 11:       3       Abies shastensis
head(inputs$obs_dt)      # per site, year and species: what fit(data = ) takes
#>    siteID  year    ba   dbh trees growth growth_n  mort n_at_risk n_died   reg
#>     <int> <int> <num> <num> <num>  <num>    <int> <num>     <int>  <int> <num>
#> 1:      1    10     0   NaN     0    NaN        0    NA         0      0     0
#> 2:      1    20     0   NaN     0    NaN        0    NA         0      0     0
#> 3:      2    10     0   NaN     0    NaN        0    NA         0      0     0
#> 4:      2    20     0   NaN     0    NaN        0    NA         0      0     0
#> 5:      3    10     0   NaN     0    NaN        0    NA         0      0     0
#> 6:      3    20     0   NaN     0    NaN        0    NA         0      0     0
#>    species   species_name period_length
#>      <int>         <char>         <int>
#> 1:       1 Abies concolor            10
#> 2:       1 Abies concolor            10
#> 3:       1 Abies concolor            10
#> 4:       1 Abies concolor            10
#> 5:       1 Abies concolor            10
#> 6:       1 Abies concolor            10
range(inputs$env_dt$year) # one row per site and simulated year: what fit(env = ) takes
#> [1]  1 20
head(inputs$siteID_dt)   # the map back to the inventory's own names
#>    siteID patchID siteName  patchName  year OrigYear
#>     <int>   <int>   <char>     <char> <int>    <int>
#> 1:      1       1 11_66410 11_66410_1     0     2001
#> 2:      1       2 11_66410 11_66410_2     0     2001
#> 3:      1       3 11_66410 11_66410_3     0     2001
#> 4:      1       4 11_66410 11_66410_4     0     2001
#> 5:      1       1 11_66410 11_66410_1    10     2011
#> 6:      1       2 11_66410 11_66410_2    10     2011
# with torch, the same call also builds the starting cohorts (a CohortMat)
inputs <- resolveSiteIDs(tree_dt = obs$tree_dt, obs_dt = obs$obs_dt, env_dt = env_dt)
inputs$initCohorts$dims  # sites, patches, cohorts
#> [1]  8  4 15
```
