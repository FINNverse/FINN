# Draw binomial counts from per-trial Bernoulli probabilities

The mortality draw selected with
`createProcess(mortality_draw = "bernoulli")`: every element of `n` gets
`max(n)` reparameterised Bernoulli trials with probability `p`, of which
the first `n` are kept and summed. The result is integer-valued and
matches the binomial in mean and variance, and it keeps a gradient path
to `p` (a straight-through sigmoid) and to `n` (a soft mask). The price
is memory: the trials are laid out along an extra dimension of length
`max(n)`, so the cost grows with the largest cohort.

## Usage

``` r
binomial_from_bernoulli(n, p)
```

## Arguments

- n:

  torch.Tensor Number of Bernoulli trials (trees at risk) per element,
  with at least two dimensions (FINN passes `[sites, patches, cohorts]`)
  and a largest value of at least 1.

- p:

  torch.Tensor Success (death) probability per element, of the shape of
  `n`.

## Value

A `torch` tensor of the shape of `n`: the number of successes per
element, a whole number between 0 and `n`.

## See also

[`binomial_from_gamma()`](https://finnverse.github.io/FINN/reference/binomial_from_gamma.md),
[`binomial_from_poisson()`](https://finnverse.github.io/FINN/reference/binomial_from_poisson.md),
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
(`mortality_draw`)

## Examples

``` r
FINN.seed(1)
# one site, 5000 patches, two cohorts of 5 and 40 trees, 30 % mortality
n <- torch::torch_tensor(array(rep(c(5, 40), each = 5000), dim = c(1, 5000, 2)))
p <- torch::torch_full(c(1, 5000, 2), 0.3)
deaths <- torch::as_array(binomial_from_bernoulli(n, p))
dim(deaths)
#> [1]    1 5000    2
# binomial mean n * p = 1.5 and 12, variance n * p * (1 - p) = 1.05 and 8.4
sapply(c("n = 5" = 1, "n = 40" = 2),
       function(k) c(mean = mean(deaths[1, , k]), var = var(deaths[1, , k])))
#>         n = 5    n = 40
#> mean 1.493600 11.940400
#> var  1.043768  8.316511
```
