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
  continuous = FALSE
)
```

## Arguments

- dbh:

  torch.Tensor Diameter at breast height for each cohort.

- species:

  torch.Tensor species index for each cohort.

- trees:

  torch.Tensor Number of trees in each cohort.

- parComp:

  torch.Tensor Competition / height-allometry parameters per species.

- h:

  torch.Tensor (Optional) Height of each cohort. Defaults to NULL.

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

## Value

torch.Tensor Fraction of available light (light) for each cohort.
