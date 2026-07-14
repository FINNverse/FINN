# Generate random numbers from a uniform distribution

This function generates random numbers from a uniform distribution with
specified low and high values and size similar to np.random.uniform in
Python.

## Usage

``` r
np_runif(low, high, size)
```

## Arguments

- low:

  numeric Lower bound of the uniform distribution.

- high:

  numeric Upper bound of the uniform distribution.

- size:

  numeric Size of the output array.

## Value

array A numeric array of random numbers.

## Examples

``` r
np_runif(0, 1, c(2, 3))
#>           [,1]      [,2]       [,3]
#> [1,] 0.7353196 0.9805397 0.05144628
#> [2,] 0.1959567 0.7415215 0.53021246
```
