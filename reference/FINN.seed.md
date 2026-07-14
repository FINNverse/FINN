# Set Seed for Reproducibility in R and Torch

This function sets the seed for both R's random number generator and
Torch's random number generator, ensuring reproducibility across
operations that involve both R and Torch.

## Usage

``` r
FINN.seed(seed)
```

## Arguments

- seed:

  An integer value to set as the seed for both R and Torch. This ensures
  that random operations in both environments produce consistent
  results.

## Value

This function does not return a value. It sets the seed internally for
both R and Torch.

## Details

The function calls [`set.seed()`](https://rdrr.io/r/base/Random.html) to
set the seed for R's random number generator and
[`torch::torch_manual_seed()`](https://torch.mlverse.org/docs/reference/torch_manual_seed.html)
to set the seed for Torch's random number generator. This is useful for
ensuring reproducibility in scripts that rely on both R and Torch for
random operations.

## Examples

``` r
if (FALSE) { # \dontrun{
FINN.seed(123)
# Now both R and Torch are seeded with 123, ensuring reproducible results
} # }
```
