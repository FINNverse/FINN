# Draw binomial counts from per-trial Bernoulli probabilities

Draw binomial counts from per-trial Bernoulli probabilities

## Usage

``` r
binomial_from_bernoulli(n, p)
```

## Arguments

- n:

  torch.Tensor Number of Bernoulli trials per element.

- p:

  torch.Tensor Success probability per element.

## Value

torch.Tensor The number of successes per element.
