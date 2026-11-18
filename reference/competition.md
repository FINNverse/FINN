# Compute the fraction of available light (light) for each cohort based on the given parameters

This function calculates the fraction of available light for each cohort
of trees based on their diameter at breast height (dbh), species, number
of trees, and global parameters.

## Usage

``` r
competition(
  dbh,
  species,
  trees,
  parComp,
  h = NULL,
  patch_size_ha,
  ba = NULL,
  cohortHeights = NULL,
  n_quantiles = 10,
  continuous = FALSE,
  light_mask = TRUE
)
```

## Arguments

- dbh:

  torch.Tensor Diameter at breast height for each cohort, in cm,
  `[sites, patches, cohorts]`.

- species:

  torch.Tensor (int64) species index for each cohort,
  `[sites, patches, cohorts]`.

- trees:

  torch.Tensor Number of trees in each cohort,
  `[sites, patches, cohorts]`.

- parComp:

  torch.Tensor Competition parameters, one row per species and two
  columns: the height parameter of
  [`height()`](https://finnverse.github.io/FINN/reference/height.md) and
  the competition strength.

- h:

  (Optional) Defaults to NULL: the light each cohort receives at its own
  height. Any non-`NULL` value returns the light at the forest floor
  instead, one value per patch (`[sites, patches, 1]`), which is what
  the regeneration process is given; the value itself is not used.

- patch_size_ha:

  numeric Patch size in hectares.

- ba:

  torch.Tensor (Optional) Pre-computed basal area. Defaults to NULL.

- cohortHeights:

  torch.Tensor (Optional) Pre-computed cohort heights. Defaults to NULL.

- n_quantiles:

  integer Number of height quantiles used when `continuous = FALSE`.
  Defaults to 10.

- continuous:

  logical Use the continuous competition formulation. Defaults to FALSE.

- light_mask:

  logical Take the height quantiles over the live cohorts only
  (`trees > 0.5`), leaving every cohort's real height as the bin query.
  Defaults to TRUE since 0.2.0: unmasked, the padded cohort width is set
  by the widest site in the batch, so the same stand's light depends on
  which other sites happen to share its batch. `FALSE` reproduces
  earlier results. See
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  for the measurement.

## Value

torch.Tensor Fraction of available light (between 0 and 1) for each
cohort, `[sites, patches, cohorts]`; with `h` given, the light at the
forest floor, `[sites, patches, 1]`.

## Details

The light of a cohort is one minus the weighted basal area of the
cohorts of its patch that are taller than it (by more than 0.1 m),
`light = max(0, 1 - sum(BA_taller * parComp[, 2] * 0.1))`, where `BA` is
the basal area in m^2/ha
([`BA_stand()`](https://finnverse.github.io/FINN/reference/BA_stand.md))
and the height comes from
[`height()`](https://finnverse.github.io/FINN/reference/height.md) with
the species' height parameter `parComp[, 1]`. With a competition
strength of 0.2 (about the start value) a cohort under 50 m^2/ha of
taller trees receives no light. In a model the function is the
competition process, `createProcess(func = competition)`; it takes no
environment.

## See also

[`height()`](https://finnverse.github.io/FINN/reference/height.md),
[`BA_stand()`](https://finnverse.github.io/FINN/reference/BA_stand.md),
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
(`n_quantiles`, `continuous`, `light_mask`),
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md)

## Examples

``` r
# the light of three cohorts on one 0.1 ha patch: two 50-cm trees, ten 20-cm
# trees and forty 5-cm trees of one species (height parameter 0.5,
# competition strength 0.2)
shape   <- c(1, 1, 3)
dbh     <- torch::torch_tensor(array(c(50, 20, 5), shape))
trees   <- torch::torch_tensor(array(c(2, 10, 40), shape))
species <- torch::torch_tensor(array(1L, shape), dtype = torch::torch_long())
parComp <- torch::torch_tensor(matrix(c(0.5, 0.2), nrow = 1))
# the tallest cohort is in full light, each smaller one is shaded by the
# basal area above it
competition(dbh, species, trees, parComp, patch_size_ha = 0.1, continuous = TRUE)
#> torch_tensor
#> (1,.,.) = 
#>   1.0000  0.9215  0.8586
#> [ CPUFloatType{1,1,3} ]
# the light that reaches the forest floor, where the recruits start
competition(dbh, species, trees, parComp, h = 0, patch_size_ha = 0.1)
#> torch_tensor
#> (1,.,.) = 
#>   0.8429
#> [ CPUFloatType{1,1,1} ]

# the light model as the competition process of a model. It takes no
# environment; `continuous = TRUE` compares every pair of cohorts instead of
# binning them into height quantiles
m <- finn(N_species = 2,
          competition_process = createProcess(func = competition, continuous = TRUE))

# a 30-year simulation from bare ground; an unfitted model uses `env` as
# given (no scaling constants yet), so the predictor is passed standardised
env <- data.frame(siteID = 1L, year = 1:30, temp = 0)
sim <- simulateForest(m, env = env, patches = 5, patch_size = 0.1)
tail(sim$wide$site)
#>    siteID species  year      dbh        ba trees    growth       mort   reg
#>     <int>   <int> <int>    <num>     <num> <num>     <num>      <num> <num>
#> 1:      1       2    25 41.34618 0.3323392   1.6 0.2630191 0.01891340   0.2
#> 2:      1       2    26 39.42628 0.3603613   1.8 0.3062723 0.01889439   0.2
#> 3:      1       2    27 42.14650 0.3896927   1.8 0.3313038 0.01886775   0.0
#> 4:      1       2    28 40.76972 0.3218157   1.6 0.2931730 0.01885694   0.0
#> 5:      1       2    29 44.11386 0.3501407   1.6 0.2866735 0.01888959   0.0
#> 6:      1       2    30 37.63353 0.3489603   1.8 0.2344551 0.01888839   0.4
#>    r_mean_ha
#>        <num>
#> 1: 0.7281690
#> 2: 0.7196229
#> 3: 0.7095281
#> 4: 0.7148335
#> 5: 0.7074409
#> 6: 0.7030461
```
