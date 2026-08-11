# Calculate growth

This function calculates growth based on specified parameters.

## Usage

``` r
growth(
  dbh,
  species,
  parGrowth,
  pred,
  light,
  light_steepness = 10,
  debug = FALSE,
  trees = NULL
)
```

## Arguments

- dbh:

  torch.Tensor Diameter at breast height.

- species:

  torch.Tensor species of tree.

- parGrowth:

  torch.Tensor Growth parameters.

- pred:

  torch.Tensor Predicted values.

- light:

  torch.Tensor Accumulated Light.

- light_steepness:

  numeric Steepness of the light-response sigmoid. Defaults to 10.

- debug:

  logical If TRUE, return the intermediate components as a list.
  Defaults to FALSE.

- trees:

  torch.Tensor (Optional) Number of trees per cohort. Defaults to NULL.

## Value

torch.Tensor A tensor representing the forest plot growth.
