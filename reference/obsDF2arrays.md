# Convert observation data frame to arrays

Convert observation data frame to arrays

## Usage

``` r
obsDF2arrays(obs_dt, additional_cols = character(0))
```

## Arguments

- obs_dt:

  data.frame The observation data table containing siteID, patchID,
  cohortID, species, dbh, and trees columns.

- additional_cols:

  character vector Optional. Additional columns to be included as
  arrays.

## Value

A list of arrays for species, dbh, trees, and additional columns.

## Examples

``` r
obs_dt <- data.frame(siteID = c(1, 1, 2), patchID = c(1, 2, 1), cohortID = c(1, 1, 2), species = c(1, 2, 1), dbh = c(10, 20, 30), trees = c(100, 200, 150), height = c(5, 10, 15))
result <- obsDF2arrays(obs_dt, additional_cols = c("height"))
```
