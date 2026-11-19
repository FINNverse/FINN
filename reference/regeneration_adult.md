# Regeneration limited by conspecific adult density

A mechanistic extension of
[regeneration_capped](https://finnverse.github.io/FINN/reference/regeneration_capped.md)
in which recruitment needs a seed source: the mean recruitment of a
species scales with the density of its own ADULTS on the site, the stems
whose diameter exceeds an adult threshold `x`.

## Usage

``` r
regeneration_adult(
  species,
  parReg,
  pred,
  light,
  dbh = NULL,
  trees = NULL,
  debug = FALSE
)
```

## Arguments

- species:

  torch.Tensor species information.

- parReg:

  torch.Tensor Regeneration parameters. 0 \<= parReg \<= 1 This
  parameter denotes the fraction of light needed for a species to
  regenerate. In general low values for high regeneration and high
  values for low regeneration.

- pred:

  torch.Tensor Prediction values.

- light:

  torch.Tensor Available light variable for calculation.

- dbh:

  torch.Tensor cohort diameters, `[sites, patches, cohorts]`.

- trees:

  torch.Tensor stems per cohort, `[sites, patches, cohorts]`.

- debug:

  logical If TRUE, return the intermediate components as a list.
  Defaults to FALSE.

## Value

torch.Tensor mean recruitment, stems ha^-1 step^-1,
`[sites, patches, species]`. With `debug = TRUE`, a list that also holds
`adults` (A, stems ha^-1, `[sites, 1, species]`) and `seed` (f(A)).

## Details

For species `k` on a site, \$\$A_k = \frac{1}{P a}
\sum\_{\mathrm{cohorts\\ of\\ } k} n \\ \sigma\\\left(\frac{d -
x_k}{\tau}\right)\$\$ is the conspecific adult density (stems ha^-1),
pooled over the `P` patches of the site (area `a` each) because seeds
disperse beyond one patch. The soft threshold (`tau` = 2 dbh units)
keeps `x_k` differentiable; a hard `d > x` would give it no gradient.
The seed-limitation term \$\$f(A_k) = \frac{A_k}{A_k + A\_{1/2,k}}\$\$
is 0 without adults and approaches 1 once adults are abundant, so the
mean recruitment \$\$m = r(\mathrm{light}) \\ e^{\mathrm{env}} f(A) +
\mathrm{reg\\floor}\$\$ falls to the floor (read it as immigration from
outside the site) where the species has no adult, and rises to the
light- and environment-driven potential of
[regeneration](https://finnverse.github.io/FINN/reference/regeneration.md)
where adults are many. The result is capped with the same Beverton-Holt
term as
[regeneration_capped](https://finnverse.github.io/FINN/reference/regeneration_capped.md),
`K m / (K + m)`.

The cap `reg_logK` is created by
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md) as for
[regeneration_capped](https://finnverse.github.io/FINN/reference/regeneration_capped.md),
with its resolution from
[`createProcess(cap = )`](https://finnverse.github.io/FINN/reference/createProcess.md)
and its prior from `prior`/`prior_lambda`, so it need not be declared. A
`reg_logK` declared in `custom_parameters` (as in the example) is kept
and sets the resolution by its length, but its start value is replaced
by the resolved prior at the first
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md). Declare
the two adult parameters with
[`createProcess(custom_parameters = )`](https://finnverse.github.io/FINN/reference/createProcess.md);
a length-1 init gives one shared value and a length-`N_species` init one
value per species:

- `reg_logK` (optional, see above): the log Beverton-Holt cap, stems
  ha^-1 step^-1;

- `reg_adult_logx`: the log adult threshold, `x = exp(reg_adult_logx)`
  in the units of `dbh`. Once `x` falls below the smallest diameter
  present, every stem counts as an adult and `x` stops receiving
  gradient;

- `reg_adult_logA`: the log half-saturation adult density,
  `A_half = exp(reg_adult_logA)` stems ha^-1.

`forward()` passes `dbh` and `trees` to a regeneration function only
when its signature has them, so the other regeneration functions are
unaffected. Called without `dbh`/`trees` (as the ALE code in `xAI.R`
does), the function returns recruitment with the seed source present
(`f = 1`).

## See also

[regeneration_capped](https://finnverse.github.io/FINN/reference/regeneration_capped.md),
[`regeneration_adult_effect()`](https://finnverse.github.io/FINN/reference/regeneration_adult_effect.md)
(adult density as an effect, not a gate),
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md),
[recruitment](https://finnverse.github.io/FINN/reference/recruitment.md)
for the recruitment options and how to combine them

## Examples

``` r
library(data.table)
# three sites of the bundled Oregon FIA sample (?fia_obs_dt)
sites <- 1:3
Nsp   <- max(fia_obs_dt$species)
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)

FINN.seed(1)
m <- finn(
  N_species = Nsp, recruits_dbh = 12.9,
  regeneration_process = createProcess(
    ~ temp + prec, FINN::regeneration_adult,
    custom_parameters = list(reg_logK       = rep(log(50), Nsp),
                             reg_adult_logx = rep(log(20), Nsp),   # adult from 20 cm
                             reg_adult_logA = rep(log(20), Nsp)),  # f = 0.5 at 20 adults/ha
    optimizeSpecies = TRUE, optimizeEnv = TRUE))
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06, plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.768)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%
# the fitted cap, adult threshold (cm) and half-saturation adult density (per ha)
exp(as.numeric(m$reg_logK)); exp(as.numeric(m$reg_adult_logx)); exp(as.numeric(m$reg_adult_logA))
#>  [1] 1.687334 1.633717 1.662593 1.662593 1.700353 1.662593 1.662405 1.662593
#>  [9] 1.662593 1.662593 1.634049
#>  [1] 19.83587 20.39247 20.00000 20.00000 19.60848 20.00000 20.33682 20.00000
#>  [9] 20.00000 20.00000 20.40406
#>  [1] 20.40364 19.60620 20.00000 20.00000 19.60513 20.00000 20.33681 20.00000
#>  [9] 20.00000 20.00000 20.40423
```
