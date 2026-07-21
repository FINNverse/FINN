# group by mean

group by mean

## Usage

``` r
groupby_mean(values, labels)
```

## Arguments

- values:

  list of value tensors

- labels:

  labels

## Value

A `torch` tensor of the mean of `values` within each unique group in
`labels`.
