# Oregon FIA calibration data: species coding

The species index shared by all `fia_*` tables. Derived from the
training sites, so `species` is `1..11` in every table.

## Usage

``` r
fia_species_dt
```

## Format

A `data.table` with 11 rows and 2 columns, `species` (integer index) and
`species_name` (Latin name, or `"other"` for the pooled rare species).

## Source

A seeded subsample of Oregon plots of the US Forest Inventory and
Analysis (FIA) database, with WorldClim bioclimatic variables attached
per plot. The preparation script and the provenance of the source tables
are in the package repository (`dev/make_extdata.R`,
`data-raw/README.md`).

## See also

Other fia datasets:
[`fia_env_dt`](https://finnverse.github.io/FINN/reference/fia_env_dt.md),
[`fia_init_trees`](https://finnverse.github.io/FINN/reference/fia_init_trees.md),
[`fia_obs_dt`](https://finnverse.github.io/FINN/reference/fia_obs_dt.md),
[`fia_tree_dt`](https://finnverse.github.io/FINN/reference/fia_tree_dt.md)

## Examples

``` r
library(data.table)
fia_species_dt
#>     species            species_name
#>       <int>                  <char>
#>  1:       1   Pseudotsuga menziesii
#>  2:       2         Pinus ponderosa
#>  3:       3          Pinus contorta
#>  4:       4          Abies concolor
#>  5:       5           Abies grandis
#>  6:       6      Tsuga heterophylla
#>  7:       7       Tsuga mertensiana
#>  8:       8 Lithocarpus densiflorus
#>  9:       9             Alnus rubra
#> 10:      10          Abies amabilis
#> 11:      11                   other
# put names on a table that carries only the species index
deaths <- fia_tree_dt[, .(trees = .N, died = sum(died)), by = species]
merge(fia_species_dt, deaths, by = "species")
#> Key: <species>
#>     species            species_name trees  died
#>       <int>                  <char> <int> <num>
#>  1:       1   Pseudotsuga menziesii  3408   153
#>  2:       2         Pinus ponderosa  1485    83
#>  3:       3          Pinus contorta  1029   116
#>  4:       4          Abies concolor   688    57
#>  5:       5           Abies grandis   613    48
#>  6:       6      Tsuga heterophylla   382    15
#>  7:       7       Tsuga mertensiana   277    23
#>  8:       8 Lithocarpus densiflorus   226     7
#>  9:       9             Alnus rubra   225    22
#> 10:      10          Abies amabilis   148    13
#> 11:      11                   other  1498   124
```
