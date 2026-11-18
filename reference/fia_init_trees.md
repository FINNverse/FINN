# Oregon FIA calibration data: initial tree list

The year-0 tree list of the
[fia_obs_dt](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
training sites (`fia_init_trees`) and the
[fia_obs_test](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
holdout sites (`fia_init_test`): what a site looked like at the first
census, from which
[`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md)
builds the initial cohorts a fit or simulation starts from
(`makeInitCohorts(fia_init_trees, Nspecies = 11)`).

## Usage

``` r
fia_init_trees

fia_init_test
```

## Format

A `data.table` with one row per tree (6039 / 6083 rows) and 28 columns.
The ones FINN uses:

- `siteID`, `patchID`, `year`:

  site and patch index; `year` is 0.

- `species`, `species_name`:

  species index `1..11` and name, coded as in
  [fia_species_dt](https://finnverse.github.io/FINN/reference/fia_species_dt.md).

- `dbh`:

  diameter at breast height, cm.

- `trees`:

  trees represented by the record (`1`; `NA` for records that are not
  living stems).

- `treeName`:

  unique tree identifier.

The remaining columns (`siteName`, `patchName`, `OrigYear`, FIA's `SPCD`
and `STATUSCD`, `status`, `living`, and the per-tree change fields such
as `dbh_before`, `dbh_growth`, `mort`, `reg`) are carried over from the
remeasurement tree list produced by
[`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md);
at year 0 most of them are `NA` by construction.

An object of class `data.table` (inherits from `data.frame`) with 6083
rows and 28 columns.

## Source

A seeded subsample of Oregon plots of the US Forest Inventory and
Analysis (FIA) database, with WorldClim bioclimatic variables attached
per plot. The preparation script and the provenance of the source tables
are in the package repository (`dev/make_extdata.R`,
`data-raw/README.md`).

## See also

Other fia datasets:
[`fia_env_dt`](https://finnverse.github.io/FINN/reference/fia_env_dt.md),
[`fia_obs_dt`](https://finnverse.github.io/FINN/reference/fia_obs_dt.md),
[`fia_species_dt`](https://finnverse.github.io/FINN/reference/fia_species_dt.md),
[`fia_tree_dt`](https://finnverse.github.io/FINN/reference/fia_tree_dt.md)

## Examples

``` r
library(data.table)
head(fia_init_trees[, .(siteID, patchID, species, species_name, dbh, trees, treeName)])
#>    siteID patchID species          species_name     dbh trees      treeName
#>     <int>   <int>   <int>                <char>   <num> <int>        <char>
#> 1:      1       1      11                 other 120.650    NA 1_50638_1_101
#> 2:      1       1      11                 other  83.566    NA 1_50638_1_102
#> 3:      1       1      11                 other 102.870    NA 1_50638_1_104
#> 4:      1       1      11                 other 112.776    NA 1_50638_1_105
#> 5:      1       1      11                 other  70.104    NA 1_50638_1_106
#> 6:      1       2       1 Pseudotsuga menziesii  80.264     1 1_50638_2_121
# living stems and their mean diameter per species at the first inventory
fia_init_trees[!is.na(trees), .(stems = sum(trees), mean_dbh = round(mean(dbh), 1)),
               by = species_name][order(-stems)]
#>                species_name stems mean_dbh
#>                      <char> <int>    <num>
#>  1:   Pseudotsuga menziesii  1862     48.9
#>  2:                   other   849     30.6
#>  3:         Pinus ponderosa   792     39.9
#>  4:          Pinus contorta   516     20.0
#>  5:          Abies concolor   374     30.5
#>  6:           Abies grandis   306     33.4
#>  7:      Tsuga heterophylla   197     44.1
#>  8:             Alnus rubra   142     28.2
#>  9:       Tsuga mertensiana   140     30.1
#> 10: Lithocarpus densiflorus   100     21.8
#> 11:          Abies amabilis    72     24.9
# the initial state of a fit or a simulation, here for three sites: one
# cohort per tree
init <- makeInitCohorts(fia_init_trees[siteID %in% 1:3], Nspecies = 11)
init$dims   # sites, patches, cohorts
#> [1]  3  4 16
```
