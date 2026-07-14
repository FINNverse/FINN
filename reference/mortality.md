# Mortality

Mortality

## Usage

``` r
mortality(
  dbh,
  species,
  trees,
  parMort,
  pred,
  light,
  base_steepness = 5,
  debug = F,
  growth = NULL
)
```

## Arguments

- dbh:

  dbh

- species:

  species

- trees:

  trees

- parMort:

  parMort

- pred:

  predictions

- light:

  available light

- base_steepness:

  numeric Steepness of the shade-response sigmoid. Defaults to 5.

- debug:

  logical If TRUE, return the intermediate components as a list.
  Defaults to FALSE.

- growth:

  torch.Tensor (Optional) Growth entering the mortality response;
  defaults to the model's current growth.
