# Convert a climate data frame to a FINN environment array

Reshapes a long-format climate table into the site x year x variable
array FINN uses internally, auto-detecting the time resolution from the
columns.

## Usage

``` r
climateDF2array(climate_dt, env_vars)
```

## Arguments

- climate_dt:

  A data.frame/data.table of climate values in long format, with
  `siteID`, a time column (e.g. `year`), and one column per variable.

- env_vars:

  character Names of the environmental variable columns to extract.

## Value

A numeric array of environmental values indexed by site, time and
variable.
