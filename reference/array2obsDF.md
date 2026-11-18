# Transform Arrays to Observation Data Table

This function transforms arrays of species, dbh, and trees back into a
cohort table, the inverse of
[`obsDF2arrays()`](https://finnverse.github.io/FINN/reference/obsDF2arrays.md).

## Usage

``` r
array2obsDF(obs_array)
```

## Arguments

- obs_array:

  A list containing three arrays `species`, `dbh` and `trees`, each of
  dimension sites x patches x cohorts.

## Value

A data.frame with one row per cohort and the columns `siteID`,
`patchID`, `cohortID` (the positions in the three array dimensions),
`species`, `dbh` and `trees`. Slots with `dbh` of 0 or `NA` (no cohort)
are left out.

## See also

[`obsDF2arrays()`](https://finnverse.github.io/FINN/reference/obsDF2arrays.md),
[CohortMat](https://finnverse.github.io/FINN/reference/CohortMat.md)

## Examples

``` r
# 2 sites x 2 patches x 2 cohort slots; a dbh of 0 marks an empty slot
obs_array <- list(
  species = array(c(1, 2), dim = c(2, 2, 2)),
  dbh = array(c(10, 20, 30, 0), dim = c(2, 2, 2)),
  trees = array(c(100, 200, 150, 0), dim = c(2, 2, 2)))
result <- array2obsDF(obs_array)
result   # six cohorts: the two empty slots are dropped
#>   siteID patchID cohortID species dbh trees
#> 1      1       1        1       1  10   100
#> 2      1       1        2       1  10   100
#> 3      1       2        1       1  30   150
#> 4      1       2        2       1  30   150
#> 5      2       1        1       2  20   200
#> 6      2       1        2       2  20   200
```
