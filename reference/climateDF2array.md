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
variable: `[site, year, variable]` for a yearly table,
`[site, year, month, variable]` with a `month` column and
`[site, year, month, day, variable]` with a `day` column. Sites and
years are in increasing order; the array has no dimnames, and the year
axis holds only the years present in the table. A combination of
`siteID`, `year` (and `month`, `day` where present) that appears more
than once is an error, whether or not the rows are identical. Missing
(`siteID`, `year`) combinations become `NA` in the array.

## See also

[finn_data](https://finnverse.github.io/FINN/reference/finn_data.md) for
the `env` table of
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md), which is
converted with this function

## Examples

``` r
library(data.table)
# yearly climate of three sites over 22 years -> array [site, year, variable]
env <- fia_env_dt[siteID %in% 1:3]
a <- climateDF2array(env, env_vars = c("temp", "prec"))
dim(a)
#> [1]  3 22  2
head(a[1, , ])                # site 1; the rows are the years 1, 2, ...
#>          [,1] [,2]
#> [1,] 3.958333  547
#> [2,] 3.958333  547
#> [3,] 3.958333  547
#> [4,] 3.958333  547
#> [5,] 3.958333  547
#> [6,] 3.958333  547

# a monthly table (siteID, year, month) gives a fourth dimension of 12 months
monthly <- CJ(siteID = 1:2, year = 2001:2002, month = 1:12)
monthly[, temp := 10 + 8 * sin((month - 4) / 12 * 2 * pi)]
#> Key: <siteID, year, month>
#>     siteID  year month      temp
#>      <int> <int> <int>     <num>
#>  1:      1  2001     1  2.000000
#>  2:      1  2001     2  3.071797
#>  3:      1  2001     3  6.000000
#>  4:      1  2001     4 10.000000
#>  5:      1  2001     5 14.000000
#>  6:      1  2001     6 16.928203
#>  7:      1  2001     7 18.000000
#>  8:      1  2001     8 16.928203
#>  9:      1  2001     9 14.000000
#> 10:      1  2001    10 10.000000
#> 11:      1  2001    11  6.000000
#> 12:      1  2001    12  3.071797
#> 13:      1  2002     1  2.000000
#> 14:      1  2002     2  3.071797
#> 15:      1  2002     3  6.000000
#> 16:      1  2002     4 10.000000
#> 17:      1  2002     5 14.000000
#> 18:      1  2002     6 16.928203
#> 19:      1  2002     7 18.000000
#> 20:      1  2002     8 16.928203
#> 21:      1  2002     9 14.000000
#> 22:      1  2002    10 10.000000
#> 23:      1  2002    11  6.000000
#> 24:      1  2002    12  3.071797
#> 25:      2  2001     1  2.000000
#> 26:      2  2001     2  3.071797
#> 27:      2  2001     3  6.000000
#> 28:      2  2001     4 10.000000
#> 29:      2  2001     5 14.000000
#> 30:      2  2001     6 16.928203
#> 31:      2  2001     7 18.000000
#> 32:      2  2001     8 16.928203
#> 33:      2  2001     9 14.000000
#> 34:      2  2001    10 10.000000
#> 35:      2  2001    11  6.000000
#> 36:      2  2001    12  3.071797
#> 37:      2  2002     1  2.000000
#> 38:      2  2002     2  3.071797
#> 39:      2  2002     3  6.000000
#> 40:      2  2002     4 10.000000
#> 41:      2  2002     5 14.000000
#> 42:      2  2002     6 16.928203
#> 43:      2  2002     7 18.000000
#> 44:      2  2002     8 16.928203
#> 45:      2  2002     9 14.000000
#> 46:      2  2002    10 10.000000
#> 47:      2  2002    11  6.000000
#> 48:      2  2002    12  3.071797
#>     siteID  year month      temp
#>      <int> <int> <int>     <num>
dim(climateDF2array(monthly, env_vars = "temp"))
#> [1]  2  2 12  1
```
