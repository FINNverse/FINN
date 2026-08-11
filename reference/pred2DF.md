# Convert Prediction Arrays to Data Frames

This function takes prediction arrays from a model output and converts
them into data frames. The data frames can be returned in either 'wide'
or 'long' format. The function processes site-level, patch-level, and
cohort-level predictions.

## Usage

``` r
pred2DF(pred, format = "wide")
```

## Arguments

- pred:

  A list containing prediction arrays for site-level, patch-level, and
  cohort-level data. Each element of `pred` should have a structure
  similar to what is outlined in the details section.

- format:

  A character string indicating the desired format of the output data
  frames. Must be either "wide" (default) or "long".

## Value

A list of data frames. The list may contain up to three elements:
`site`, `patch`, and `cohort`, corresponding to the processed
site-level, patch-level, and cohort-level predictions, respectively.

## Details

The `pred` argument should be a list containing at least a `Predictions`
element, which itself is a list of arrays. The arrays represent
predictions for different metrics such as dbh/ba, tree counts, AL
(aboveground live biomass), growth rates, mortality rates, and
regeneration rates for sites, patches, or cohorts. The dimensionality of
the arrays should correspond to different factors, such as siteID, year,
species, and optionally patch or cohortID.

The function first converts each prediction array into a data frame,
properly naming and converting the relevant dimensions. It then merges
these data frames by common identifiers such as siteID, year, species,
patch, and cohortID. Depending on the `format` parameter, the data
frames are returned in either a wide format (one row per
site/patch/cohort per year with multiple columns for different metrics)
or a long format (one row per site/patch/cohort per year per metric).

## Examples

``` r
# a minimal prediction object of the shape pred2DF() expects
# (normally produced by predict()/simulateForest()):
metrics <- c("dbh", "ba", "trees", "growth", "mort", "reg", "r_mean_ha")
arr <- array(1, dim = c(1, 2, 3), dimnames = list(1, 1:2, 1:3))  # [site, year, species]
pred <- list(Predictions = list(Site = stats::setNames(
  lapply(metrics, function(m) arr), metrics)))
result <- pred2DF(pred, format = "long")
head(result$site)
#>    siteID  year species variable value
#>     <int> <int>   <int>   <fctr> <num>
#> 1:      1     1       1      dbh     1
#> 2:      1     2       1      dbh     1
#> 3:      1     1       2      dbh     1
#> 4:      1     2       2      dbh     1
#> 5:      1     1       3      dbh     1
#> 6:      1     2       3      dbh     1
```
