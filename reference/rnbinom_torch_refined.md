# Negative-binomial recruit draw without the recruitment floor

The default regeneration sampler since 0.2.0, selected with
`createProcess(regeneration_sampler = "refined")`; `rnbinom_torch()`
(`regeneration_sampler = "gaussian"`) reproduces results published
before 0.2.0. Three differences: the probability clamp is `1e-10`
instead of `1e-4`, the Gamma draw is clamped at that same `eps` instead
of taking `+0.001` and an absolute value, and the draw goes through
`sample_poisson_refined()` (an exact `torch_poisson` below `lambda = 2`,
the Gaussian approximation above) instead of
`sample_poisson_gaussian()`.

## Usage

``` r
rnbinom_torch_refined(pred, theta)
```

## Arguments

- pred:

  torch.Tensor Expected recruits per patch.

- theta:

  torch.Tensor Negative-binomial dispersion.

## Value

A `torch` tensor of recruit counts, of the shape of `pred`.

## Details

The `+0.001` offset and the `1e-4` clamp exist for the Gaussian branch:
they keep the Gamma draw off zero and stop `log(0)` sending a NaN
gradient back through `distr_gamma()$rsample()`. Neither hazard applies
to `torch_poisson`, and together they impose a recruitment FLOOR of
0.0011 per cell-year, which at a FINN patch of 0.0168 ha is 0.066
recruits ha^-1 yr^-1 *per species* – over a pool of a hundred-odd
species, as large as the whole ingrowth a temperate inventory records.
`sample_poisson_gaussian()` (`lambda + eps*sqrt(lambda)`, rounded
straight-through by the forward pass) rounds the floor away but also
rounds the signal away: at a FINN patch of 0.0168 ha the per-cell mean
is 0.001-0.008, and over 2e6 draws the delivered fraction of the
calibrated mean was 0.000-0.048 at `theta = 1`.

One bias is left and is not repaired in 0.2.0: above `lambda = 2` the
draw is the Gaussian approximation FLOORED, so it delivers about half a
recruit per draw too few (measured over 1e6 draws: -0.43 at
`lambda = 2.5`, -0.50 from `lambda = 20`, against about 0 with
rounding). It bites only where a cell's mean recruitment is above two,
which the per-patch means above are not. Rounding instead of flooring is
planned for 0.2.1; it changes results, so it is not in this release.

For the measured consequence of both repairs on the simulated stand, and
why the holdout cannot separate them, see
[`binomial_from_poisson()`](https://finnverse.github.io/FINN/reference/binomial_from_poisson.md).

## See also

[recruitment](https://finnverse.github.io/FINN/reference/recruitment.md),
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
(`regeneration_sampler`, `theta_floor`)

## Examples

``` r
FINN.seed(1)
# sparse recruitment: an expected 0.05 recruits per patch and year
mu    <- torch::torch_full(c(100000), 0.05)
theta <- torch::torch_tensor(1)
r <- as.numeric(rnbinom_torch_refined(mu, theta))
mean(r)    # the draw delivers the expected count ...
#> [1] 0.05034
table(r)   # ... as whole recruits, most patches getting none
#> r
#>     0     1     2     3 
#> 95213  4548   231     8 
```
