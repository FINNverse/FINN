# Transform Arrays to Observation Data Table

This function transforms arrays of species, dbh, and trees back into an
observation data table.

## Usage

``` r
array2obsDF(obs_array)
```

## Arguments

- obs_array:

  A list containing three arrays: species, dbh, and trees.

## Value

A data.frame with columns siteID, patchID, cohortID, species, dbh, and
trees.

## Examples

``` r
obs_array <- list(
  species = array(c("A", "B"), dim = c(2, 2, 2)),
  dbh = array(c(10, 20, 30, 40), dim = c(2, 2, 2)),
  trees = array(c(100, 200, 150, 250), dim = c(2, 2, 2)))
result <- array2obsDF(obs_array)
```
