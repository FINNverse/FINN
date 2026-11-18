# Draw deaths as a Bernoulli/Poisson count

The alternative mortality draw, selected with
`createProcess(mortality_draw = "poisson")`. The default since 0.2.0 is
`mortality_draw = "binomial"`, the exact binomial;
[`binomial_from_gamma()`](https://finnverse.github.io/FINN/reference/binomial_from_gamma.md)
(`mortality_draw = "gamma"`) reproduces results published before 0.2.0.
A cohort of one tree or fewer (`n <= 1.5`) is a Bernoulli trial and is
drawn as one; every larger cohort is drawn from a Poisson capped at `n`.

## Usage

``` r
binomial_from_poisson(n, p)
```

## Arguments

- n:

  torch.Tensor Number of trees at risk per cohort.

- p:

  torch.Tensor Annual mortality probability per cohort.

## Value

A `torch` tensor of integer death counts, never exceeding `n`.

## Details

[`binomial_from_gamma()`](https://finnverse.github.io/FINN/reference/binomial_from_gamma.md)
matches a Gamma to the binomial's first two moments and `rsample()`s it,
which buys a pathwise gradient – the hazard it guards against is a
discrete draw having no gradient at all. The forward pass then rounds
that draw straight-through (`finn-class.R`), so the gradient flows
through a quantity the simulation discards, and the rounded count
delivers a measured 0.68-0.72 of the calibrated deaths for one-tree
cohorts at annual probabilities of 0.01-0.10 (5e5 draws per probability;
0.76 at `p = 0.26`). This draw is integer-valued, exact in expectation
for small `n*p`, and carries no pathwise gradient, which is how Lin et
al. (2026) treat discrete events in a differentiable crop model: resolve
them before differentiation and hold them fixed during gradient
evaluation. The Bernoulli branch is not cosmetic – capping a Poisson at
a one-tree cohort delivers only 88\\ suppressed stem, which is the stem
that decides thinning.

The choice shows up in the simulated stand rather than in a holdout
score: a long simulation to equilibrium separates the two draws by a
factor of several in basal area and stem number, while the fitted loss
moves far less, because the rounded deficit is a bias on the trajectory
and not on the parameters.

## References

Lin et al. (2026).
[doi:10.48550/arXiv.2610.36827](https://doi.org/10.48550/arXiv.2610.36827)

## See also

[`binomial_from_gamma()`](https://finnverse.github.io/FINN/reference/binomial_from_gamma.md),
[`binomial_from_bernoulli()`](https://finnverse.github.io/FINN/reference/binomial_from_bernoulli.md),
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
(`mortality_draw`)

## Examples

``` r
FINN.seed(1)
# cohorts of one tree are Bernoulli trials, larger cohorts a capped Poisson
n <- torch::torch_tensor(rep(c(1, 40), each = 10000))
p <- torch::torch_full(c(20000), 0.1)
d <- as.numeric(binomial_from_poisson(n, p))
# mean deaths per cohort: expected 0.1 (one tree) and 4 (forty trees)
tapply(d, rep(c("n = 1", "n = 40"), each = 10000), mean)
#>  n = 1 n = 40 
#> 0.1009 4.0245 
range(d)   # whole trees, and never more than the cohort holds
#> [1]  0 13
```
