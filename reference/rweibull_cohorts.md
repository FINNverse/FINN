# Generate Cohorts Using Weibull Distribution

This function generates cohort data for tree populations using the
Weibull distribution based on specified tree counts, diameter at breast
height (DBH) shape, and scale parameters.

## Usage

``` r
rweibull_cohorts(
  trees = NULL,
  dbh_shape = NULL,
  dbh_scale = NULL,
  dbh_class_range = 1,
  siteID = 1,
  patchID = 1,
  species = 1
)
```

## Arguments

- trees:

  Integer vector specifying the number of trees for each cohort. If a
  single integer is provided, it will be replicated for each draw.

- dbh_shape:

  Numeric vector specifying the shape parameters of the Weibull
  distribution for DBH for each cohort. If a single numeric value is
  provided, it will be replicated for each draw.

- dbh_scale:

  Numeric vector specifying the scale parameters of the Weibull
  distribution for DBH for each cohort. If a single numeric value is
  provided, it will be replicated for each draw.

- dbh_class_range:

  Numeric value specifying the range of DBH classes. Default is 1.

- siteID:

  Integer vector specifying the site ID for each cohort. If a single
  integer is provided, it will be replicated for each draw. Default is
  1.

- patchID:

  Integer vector specifying the patch ID for each cohort. If a single
  integer is provided, it will be replicated for each draw. Default is
  1.

- species:

  Integer vector specifying the species ID for each cohort. If a single
  integer is provided, it will be replicated for each draw. Default is
  1.

## Value

A data frame containing the cohort data with columns for site ID, patch
ID, cohort ID, species, number of trees, and DBH.

## Details

The function generates cohort data by drawing samples from the Weibull
distribution for each cohort based on the specified shape and scale
parameters. The resulting DBH values are binned into classes, and the
cohort data is generated accordingly.

## Examples

``` r
obs_df <- rweibull_cohorts(
  trees = c(300, 10),
  dbh_shape = c(3, 3),
  dbh_scale = c(1, 50),
  dbh_class_range = 0.1,
  siteID = c(1, 1),
  patchID = c(1, 1),
  species = c(3, 4)
)
head(obs_df)
#>   siteID patchID cohortID species dbh trees
#> 1      1       1        1       3 0.4     5
#> 2      1       1        2       3 0.5    10
#> 3      1       1        3       3 0.6    14
#> 4      1       1        4       3 0.7    21
#> 5      1       1        5       3 0.8    30
#> 6      1       1        6       3 0.9    36
```
