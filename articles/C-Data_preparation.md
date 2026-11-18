# Preparing your data for FINN

This vignette is a standalone reference for preparing forest-inventory
data into FINN. It documents the inputs of FINN: columns names and
units, and the helper functions for data preparation
[`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
→
[`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
→
[`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md).

``` r

library(FINN)
library(data.table)
```

## What FINN needs

FINN is calibrated against clearly defined inputs:

- **`obs_dt`**: demographic responses per site × year × species:
  `siteID, year, ba, dbh, trees, growth, mort, n_at_risk, n_died, reg, species, species_name, period_length`.
  FINN simulates the forest in **yearly steps**, so `year` counts years
  since the first inventory, which is the initial state (year 0) and
  *not* included; with inventories ten years apart the observations sit
  at years 10 and 20, each with `period_length = 10`, the years since
  the previous inventory. NAs are allowed in any response. Units: `ba`
  m², `dbh` cm, `trees` count, `growth` relative growth rate per year,
  `reg` trees ha⁻¹ over the interval. Mortality comes back as a **pair
  of counts**, `n_at_risk` (trees alive at the start of the interval)
  and `n_died` (how many of them were dead at the end), with the yearly
  rate `mort = 1 - (1 - n_died / n_at_risk)^(1 / period_length)` derived
  from them. The pair is the `cbind(died, survived)` response of a
  binomial GLM: the default
  `loss_family = c(..., mortality = "binomial", ...)` weights each
  observation by `n_at_risk`, so a 60-tree observation outweighs a
  1-tree one.
- **`env_dt`**: `siteID, year, <env vars…>` in **natural units**, one
  row per site and simulated year, from year 1 to the last inventory. No
  manual scaling is needed: FINN z-standardizes the predictors
  internally when fitting (`env_autoscale = TRUE`, the default) and
  reuses the same constants at prediction time.
- **`init_trees`**: the initial state used to build the starting
  cohorts: `siteID, patchID, year, species, dbh, trees, living`.
- **`dist_dt`** (optional): `siteID, year, intensity` (fraction of
  patches disturbed that year).

You assemble these from two raw ingredients: a **tree list** (one row
per tree × inventory) and a **site environment table**. The rest of this
vignette builds the FINN tables from those two.

## Start from a raw tree list

``` r

# An example from raw data: Oregon FIA sites, 4 patches each, 3 inventories.
# example_tree_dt.csv and example_env_dt.csv (the climate, loaded later) are both
# built by dev/make_extdata.R from the data-raw/ source; see data-raw/README.md.
tree_dt <- fread(system.file("extdata", "example_tree_dt.csv", package = "FINN"))
str(tree_dt)
#> Classes 'data.table' and 'data.frame':   500 obs. of  13 variables:
#>  $ siteName     : chr  "11_66410" "11_66410" "11_66410" "11_66410" ...
#>  $ patchName    : chr  "11_66410_1" "11_66410_1" "11_66410_1" "11_66410_2" ...
#>  $ treeName     : chr  "11_66410_1_100" "11_66410_1_100" "11_66410_1_100" "11_66410_2_103" ...
#>  $ year         : int  2001 2011 2021 2001 2011 2021 2001 2011 2021 2001 ...
#>  $ species_name : chr  "Pseudotsuga menziesii" "Pseudotsuga menziesii" "Pseudotsuga menziesii" "Pseudotsuga menziesii" ...
#>  $ dbh          : num  47.2 57.4 65 26.9 30.7 ...
#>  $ status       : chr  "new" "alive" "alive" "new" ...
#>  $ status_before: chr  "" "new" "alive" "" ...
#>  $ mort_cause   : chr  "" "" "" "" ...
#>  $ reg          : logi  TRUE NA NA TRUE NA NA ...
#>  $ mort         : logi  NA NA NA NA NA NA ...
#>  $ living       : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ complete     : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  - attr(*, ".internal.selfref")=<externalptr>
```

A raw tree record must include: `siteName, patchName, treeName, year`
(calendar year of the inventory), `species_name`, `dbh` (cm), and a
`status` of `new` / `alive` / `dead`. That is enough:
[`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
reconstructs each tree’s previous state itself, by sorting on `treeName`
and `year`. Our example also contains columns for `status_before` and a
`complete` flag marking trees from fully re-measured plots.

## Acquiring raw data (example: US FIA)

Turning a public inventory into the raw tree list above needs external
downloads, so it is shown but not run. For the US FIA program:

``` r

library(rFIA)
fia <- readFIA(dir = "data/FIA", common = TRUE, states = "OR")
# build tree_dt: dbh = DIA * 2.54 (inch -> cm); status from STATUSCD / PREV_STATUS_CD;
# mort_cause from AGENTCD. See playground/calibration/FIA/FIA_dataprep_vignette.R
```

## Define status flags

FINN’s builders key off two logical flags derived from `status`:

``` r

# reg    = newly recruited this inventory
# living = currently alive
tree_dt[, reg    := status == "new"]
tree_dt[, living := status %chin% c("alive", "new")]
tree_dt[, .(n = .N), by = .(status, reg, living)]
#>    status    reg living     n
#>    <char> <lgcl> <lgcl> <int>
#> 1:    new   TRUE   TRUE   193
#> 2:  alive  FALSE   TRUE   239
#> 3:   dead  FALSE  FALSE    68
```

There is deliberately **no tree-level `mort` flag** here.
[`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
derives mortality itself, from `status` and `living`, as a closed
cohort: for each interval it takes the trees that were alive at the
*start*, and counts how many of them were dead at the *end*. Both
columns are pinned to the tree’s state at the start of the interval,
which matters because FIA re-identifies a tree’s species between visits;
a hand-rolled rate can otherwise book a death against the species a tree
*ended* as while the denominator counted the species it *started* as,
and report more deaths than trees.

## Build `obs_dt` with `makeObsData()`

[`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
aggregates the tree list to per-patch demographic responses and lumps
rare species into a single `"other"` class. Each argument controls a
data-quality filter:

``` r

obs_list <- makeObsData(
  tree_dt           = tree_dt[complete == TRUE],  # only fully re-measured plots
  plotsize          = 0.06,            # ha per patch (scales recruitment to per-ha)
  aggregate_by_site = FALSE,           # keep patch-level resolution
  fix_period_length = 10,              # drop sites whose remeasure gap != 10 yr
  dbh_growth_thresh = c(-10, 50),      # drop sites with implausible dbh change
  NspeciesQuantile  = 0.98,            # keep species covering 98% of stems; rest -> "other"
  minNyears         = 2                # need >= 2 inventories per patch
)
obs_dt <- obs_list$obs_dt
head(obs_dt)
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
```

Use `Nspecies` (a hard cap) or `NspeciesQuantile` (smallest set covering
a stem fraction) to control how many species you model; everything else
becomes `"other"`, which FINN still tracks as a competitor.

The rates come back **per year**, because FINN simulates the forest in
yearly steps. A tree that grew from 20 to 24 cm over the ten years of an
interval has the yearly growth rate `(24 / 20)^(1 / 10) - 1 = 0.0184`,
the constant rate that compounds to the observed increment, and a cohort
of which 1 in 2 trees died over the interval has the yearly mortality
rate `1 - 0.5^(1 / 10) = 0.067`. `period_length` is the length of the
interval in years. When fitting, FINN averages the simulated yearly
growth and mortality over these ten years and sums the simulated
recruits, so `reg` stays the recruitment of the whole interval, and the
counts `n_at_risk`, `n_died` and `growth_n` are the trees behind each
observation. The growth and mortality of every observation, per year and
as compounded back over its interval:

``` r

obs_dt[n_at_risk > 0, .(siteName, year, species_name, period_length, growth, mort,
                        growth_interval = (1 + growth)^period_length - 1,
                        mort_interval   = n_died / n_at_risk)][1:6]
#>    siteName  year          species_name period_length      growth       mort
#>      <char> <int>                <char>         <int>       <num>      <num>
#> 1: 11_66410  2011 Pseudotsuga menziesii            10 0.019669783 0.00000000
#> 2: 11_66410  2021 Pseudotsuga menziesii            10 0.012542247 0.00000000
#> 3: 11_66410  2011 Pseudotsuga menziesii            10 0.011816466 0.00000000
#> 4: 11_66410  2021 Pseudotsuga menziesii            10 0.007917063 0.01326438
#> 5: 11_66410  2011 Pseudotsuga menziesii            10 0.018626114 0.00000000
#> 6: 11_66410  2021 Pseudotsuga menziesii            10 0.010948501 0.00000000
#>    growth_interval mort_interval
#>              <num>         <num>
#> 1:       0.2150538         0.000
#> 2:       0.1327434         0.000
#> 3:       0.1246501         0.000
#> 4:       0.0820516         0.125
#> 5:       0.2026744         0.000
#> 6:       0.1150397         0.000
```

`makeObsData(rates = "interval")` returns the rates over the whole
interval instead, for a model in which one simulation step is one census
interval (the set-up of FINN 0.1.0). One step per interval distorts the
dynamics, though: recruitment, competition and growth then act once per
interval, and the mortality limits of a step become limits per decade.
Keep the default.

## Environmental data

Acquisition from rasters (e.g. WorldClim) needs external files, so it is
shown but not run. The result is a table of environmental predictors in
**natural units**, one row per site and **calendar year**, because every
simulated year needs its environment:

``` r

library(terra)
pts  <- vect(site_dt, geom = c("x", "y"), crs = "EPSG:4326")
bio1 <- rast("wc2.1_30s_bio_1.tif")
env_dt[, temp := terra::extract(bio1, project(pts, crs(bio1)))[, 2]]
```

``` r

env_dt <- fread(system.file("extdata", "example_env_dt.csv", package = "FINN"))
head(env_dt)
#>    siteName  year     temp tempmax tempmin  prec precseas precwarmq
#>      <char> <int>    <num>   <num>   <num> <int>    <num>     <int>
#> 1: 11_66410  2001 10.69167    22.8     2.7  1880 70.15234       103
#> 2: 11_66410  2002 10.69167    22.8     2.7  1880 70.15234       103
#> 3: 11_66410  2003 10.69167    22.8     2.7  1880 70.15234       103
#> 4: 11_66410  2004 10.69167    22.8     2.7  1880 70.15234       103
#> 5: 11_66410  2005 10.69167    22.8     2.7  1880 70.15234       103
#> 6: 11_66410  2006 10.69167    22.8     2.7  1880 70.15234       103
env_dt[, .(first = min(year), last = max(year), years = .N), by = siteName]
#>    siteName first  last years
#>      <char> <int> <int> <int>
#> 1: 11_66410  2001  2021    21
#> 2: 17_65658  2002  2022    21
#> 3:  1_50638  2001  2021    21
#> 4: 23_98340  2002  2022    21
#> 5: 25_81184  2002  2022    21
#> 6: 33_95570  2001  2021    21
#> 7: 35_70911  2001  2021    21
#> 8: 35_90216  2001  2021    21
```

The climate here is the long-term WorldClim normal of each plot, the
same in every year, so it is simply repeated; a yearly climate series
would give every year its own values. `env_dt` must cover every calendar
year from the first inventory of a site to its last (precisely: to the
longest first-to-last span of any site after the site’s first inventory,
since all sites share one time axis).

There is **no manual standardization step**. FINN z-scales the
predictors internally when you fit (`env_autoscale = TRUE`, the
default), stores the per-variable mean/sd on the model, and re-applies
them at prediction time. So you supply (and keep) `env_dt` in natural
units throughout.

## Resolve IDs across tables

Raw tables use names (`siteName`, `patchName`, `species_name`); FINN
works with integer indices aligned across all tables.
[`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
does that alignment and can build the initial cohorts in one pass. Pass
the **processed** tree list `obs_list$tree_dt` (not the raw one):
[`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
returns it with the same species lumping as `obs_dt`, so the two stay
consistent.

``` r

inputs <- resolveSiteIDs(
  tree_dt = obs_list$tree_dt,   # processed tree list (species lumped like obs_dt)
  obs_dt  = obs_dt,
  env_dt  = env_dt,
  createInitCohorts = TRUE
)
# obs_dt / env_dt now carry integer siteID/patchID/species, aligned across tables
head(inputs$obs_dt)
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
```

[`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
also sets the time axis. A calendar year becomes the number of years
since the site’s first inventory: the first inventory is year 0, the
initial state, and leaves `obs_dt`; the next two are observations at
years 10 and 20. `env_dt` comes back with one row per site and year 1 to
20: the row of year *k* drives the simulated year that ends *k* years
after the first inventory, so the inventory at year 10 is compared with
the stand after ten simulated years.

``` r

unique(inputs$siteID_dt[, .(siteID, year, OrigYear)])[siteID == 1]
#>    siteID  year OrigYear
#>     <int> <int>    <int>
#> 1:      1     0     2001
#> 2:      1    10     2011
#> 3:      1    20     2021
inputs$obs_dt[, .(rows = .N), by = .(year, period_length)]
#>     year period_length  rows
#>    <int>         <int> <int>
#> 1:    10            10    88
#> 2:    20            10    88
inputs$env_dt[siteID == 1][c(1:2, .N)]
#>    siteID  year     temp tempmax tempmin  prec precseas precwarmq
#>     <int> <int>    <num>   <num>   <num> <int>    <num>     <int>
#> 1:      1     1 10.69167    22.8     2.7  1880 70.15234       103
#> 2:      1     2 10.69167    22.8     2.7  1880 70.15234       103
#> 3:      1    20 10.69167    22.8     2.7  1880 70.15234       103
```

## Initial cohorts

With `createInitCohorts = TRUE`,
[`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
returns a ready `CohortMat` describing the starting state, the living
trees of the first inventory; its `dims` are
`c(sites, patches, cohorts)`:

``` r

init_cohorts <- inputs$initCohorts
init_cohorts
#> An `nn_module` containing 0 parameters.
#> 
#> ── Buffers ─────────────────────────────────────────────────────────────────────
#> • dbh: Float [1:8, 1:4, 1:15]
#> • trees: Float [1:8, 1:4, 1:15]
#> • species: Long [1:8, 1:4, 1:15]
```

To build cohorts yourself from a year-0 tree slice, call
[`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md)
directly; it expects `siteID, patchID, species, dbh, treeName` (+
optional `trees`):

``` r

init_trees <- inputs$tree_dt[living == TRUE & year == 0]
makeInitCohorts(init_trees, Nspecies = uniqueN(inputs$obs_dt$species))
```

## Checklist before fitting

- **One step is one year**: `year` counts years since the first
  inventory, `growth` and `mort` are yearly rates, and `period_length`
  gives each observation’s interval in years. The intervals may differ
  between sites (FINN then scores every site over its own interval);
  `fix_period_length` keeps only sites with one interval length if you
  prefer a uniform design.
- **The initial state** (the first inventory) is used for `init_trees`,
  and the responses in `obs_dt` start *after* it.
- **`env_dt` covers every simulated year** from 1 to the last inventory,
  for every site; missing rows break the fit.
- **Units**: `dbh` in cm, `ba` in m², `reg` in trees ha⁻¹; env
  predictors in natural units (FINN standardizes them internally).
- **Species IDs** are contiguous integers `1..N`, aligned across tables;
  always go through
  [`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
  rather than hand-coding them.

You now have `obs_dt`, `env_dt`, and `init_cohorts`, the inputs the
**Fitting FINN to forest inventory data** vignette consumes.

``` r

sessionInfo()
#> R version 4.5.0 (2025-04-11)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS 26.5.1
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Europe/Berlin
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] patchwork_1.3.0     ggplot2_3.5.2       data.table_1.18.6.1
#> [4] FINN_0.2.0         
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.6.5        cli_3.6.6          knitr_1.50         rlang_1.2.0       
#>  [5] xfun_0.57          processx_3.8.6     generics_0.1.4     torch_0.15.1      
#>  [9] coro_1.1.0         labeling_0.4.3     glue_1.8.0         bit_4.6.0         
#> [13] ps_1.9.1           scales_1.4.0       grid_4.5.0         abind_1.4-8       
#> [17] evaluate_1.0.5     tibble_3.3.0       lifecycle_1.0.5    compiler_4.5.0    
#> [21] dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.1.2         pkgconfig_2.0.3   
#> [25] farver_2.1.2       R6_2.6.1           tidyselect_1.2.1   pillar_1.11.0     
#> [29] callr_3.7.6        magrittr_2.0.3     withr_3.0.2        tools_4.5.0       
#> [33] bit64_4.6.0-1      gtable_0.3.6
```
