# Calculate the regeneration of forest patches based on the input parameters.

This function calculates the regeneration of forest patches based on
species information, regeneration parameters, prediction values, and
available light.

## Usage

``` r
regeneration(species, parReg, pred, light, debug = F)
```

## Arguments

- species:

  torch.Tensor species information.

- parReg:

  torch.Tensor Regeneration parameters. 0 \<= parReg \<= 1 This
  parameter denotes the fraction of light needed for a species to
  regenerate. In general low values for high regeneration and high
  values for low regeneration.

- pred:

  torch.Tensor Prediction values.

- light:

  torch.Tensor Available light variable for calculation.

- debug:

  logical If TRUE, return the intermediate components as a list.
  Defaults to FALSE.

## Value

torch.Tensor Regeneration values for forest patches.
