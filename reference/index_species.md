# Index species

Index species

## Usage

``` r
index_species(pred, species)
```

## Arguments

- pred:

  predictions

- species:

  species index vector, must be int64

## Value

A `torch` tensor of `pred` gathered along its species dimension by
`species`.
