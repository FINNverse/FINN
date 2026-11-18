# Draw deaths from a Gamma matched to the binomial

The mortality draw of FINN before 0.2.0, selected with
`createProcess(mortality_draw = "gamma")`: a Gamma distribution with the
mean `n * p` and the variance `n * p * (1 - p)` of the binomial, sampled
with `rsample()` so that the draw carries a pathwise gradient to `n` and
`p`. The draw is continuous, not an integer count; the simulation rounds
it afterwards. See
[`binomial_from_poisson()`](https://finnverse.github.io/FINN/reference/binomial_from_poisson.md)
for what that rounding costs and
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
for the other draws.

## Usage

``` r
binomial_from_gamma(n, p, sample_size = 1)
```

## Arguments

- n:

  torch.Tensor Number of trees at risk per cohort.

- p:

  torch.Tensor Mortality probability per cohort, strictly between 0 and
  1 (the Gamma's shape `n * p / (1 - p)` degenerates at both ends).

- sample_size:

  integer Number of draws; only the default `1` is used by FINN (the
  leading sample dimension is dropped).

## Value

A `torch` tensor of the shape of `n` with continuous, positive draws
whose mean is `n * p` and whose variance is `n * p * (1 - p)`.

## See also

[`binomial_from_poisson()`](https://finnverse.github.io/FINN/reference/binomial_from_poisson.md),
[`binomial_from_bernoulli()`](https://finnverse.github.io/FINN/reference/binomial_from_bernoulli.md),
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
(`mortality_draw`)

## Examples

``` r
FINN.seed(1)
n <- torch::torch_full(c(10000), 40)   # 10 000 cohorts of 40 trees
p <- torch::torch_full(c(10000), 0.3)
d <- as.numeric(binomial_from_gamma(n, p))
# the binomial has mean 12 and variance 8.4; the draws are not integers
c(mean = mean(d), var = var(d))
#>      mean       var 
#> 11.957824  8.343056 
head(d, 3)
#> [1] 13.76861 12.54947 10.51678
```
