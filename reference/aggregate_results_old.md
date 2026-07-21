# Aggregate function

Aggregate function

## Usage

``` r
aggregate_results_old(
  labels,
  samples,
  Results,
  drop_rows = TRUE,
  sp_max = NULL
)
```

## Arguments

- labels:

  labels

- samples:

  samples

- Results:

  Results

- drop_rows:

  logical Drop empty rows from the aggregated output. Defaults to TRUE.

- sp_max:

  integer (Optional) Maximum species index to aggregate over. Defaults
  to NULL.

## Value

A list of `torch` tensors, one per element of `samples`, each holding
the per-species aggregated result.
