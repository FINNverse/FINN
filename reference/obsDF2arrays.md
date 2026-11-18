# Convert a cohort table to arrays

Reshapes a long cohort table (one row per cohort) into the
`[sites, patches, cohorts]` arrays FINN works on. It is the step behind
`CohortMat(obs_df = )` and
[`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md);
[`array2obsDF()`](https://finnverse.github.io/FINN/reference/array2obsDF.md)
is the way back.

## Usage

``` r
obsDF2arrays(obs_dt, additional_cols = character(0))
```

## Arguments

- obs_dt:

  data.frame The cohort table, with the columns `siteID`, `patchID`,
  `cohortID` (integers; `cohortID` numbers the cohorts within a patch
  from 1), `species` (integer code), `dbh` (cm) and `trees` (number of
  trees in the cohort; rounded to whole trees).

- additional_cols:

  character vector Optional. Further numeric columns of `obs_dt` to
  return as arrays of the same shape.

## Value

A named list of arrays `species`, `dbh`, `trees` (and one per
`additional_cols`), each of dimension sites x patches x largest number
of cohorts in a patch. Slots without a cohort hold `dbh = 0`,
`trees = 0` and the placeholder `species = 1`. The site dimension
follows the sorted user `siteID`s (sorted with `method = "radix"`, so
independent of the locale); the user ids in that order are attached as
attribute `"site_ids"`, so a site set that is not `1:N` keeps its
identity. `CohortMat` stores them as `$site_ids`, and
[fit](https://finnverse.github.io/FINN/reference/fit.md) /
[predict.finn_class](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
check them against `data` and `env`.

## See also

[`array2obsDF()`](https://finnverse.github.io/FINN/reference/array2obsDF.md),
[CohortMat](https://finnverse.github.io/FINN/reference/CohortMat.md),
[`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md)

## Examples

``` r
obs_dt <- data.frame(
  siteID = c(1, 1, 2), patchID = c(1, 2, 1), cohortID = c(1, 1, 2),
  species = c(1, 2, 1), dbh = c(10, 20, 30), trees = c(100, 200, 150),
  height = c(5, 10, 15))
result <- obsDF2arrays(obs_dt, additional_cols = c("height"))
dim(result$dbh)     # 2 sites x 2 patches x 2 cohort slots
#> [1] 2 2 2
result$dbh[1, , ]   # site 1: one cohort in each patch
#>      [,1] [,2]
#> [1,]   10    0
#> [2,]   20    0
result$trees[2, , ] # site 2: patch 1 holds a cohort in slot 2, patch 2 is empty
#>      [,1] [,2]
#> [1,]    0  150
#> [2,]    0    0
```
