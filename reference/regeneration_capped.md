# Regeneration with a Beverton-Holt recruitment cap

The default regeneration process since 0.2.0, and a drop-in alternative
to the uncapped
[regeneration](https://finnverse.github.io/FINN/reference/regeneration.md).
It caps the otherwise unbounded mean recruitment with a Beverton-Holt
density dependence, `recruits = K * m / (K + m)`, where `m` is the
mechanistic mean recruitment (identical to
[regeneration](https://finnverse.github.io/FINN/reference/regeneration.md))
and `K` is a fitted carrying capacity (stems ha^-1 step^-1). As `m`
grows, recruitment saturates towards `K`. `regeneration_saturation()` is
the same function under its older name, kept as a deprecated alias for
one release.

## Usage

``` r
regeneration_capped(species, parReg, pred, light, debug = FALSE)

regeneration_saturation(species, parReg, pred, light, debug = FALSE)
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

- debug:

  logical If TRUE, return the intermediate components as a list.
  Defaults to FALSE.

## Value

torch.Tensor Mean recruitment, stems ha^-1 step^-1,
`[sites, patches, species]` – a per-HECTARE rate, not a per-patch count.
`forward()` multiplies it by `patch_size_ha` itself before the recruit
draw, and `reg` is scored on the per-hectare value.

## Details

The cap is the parameter `reg_logK` (`K = exp(reg_logK)`), which **FINN
creates itself** – there is nothing to declare.
[`createProcess(cap = )`](https://finnverse.github.io/FINN/reference/createProcess.md)
decides its resolution: `"all"` (the default) gives one cap shared by
every species, `"species"` one cap per species. Read the fitted value
back with `exp(as.numeric(m$reg_logK))`.

## The prior on the cap

A cap that only the data constrain runs to whatever the fitted window
supports, which on a national inventory is far above anything observed.
So the cap carries a prior, added to the objective as
`prior_lambda * (log K - log prior)^2` per cap – the same form as the L2
pull of a staged calibration. It is added **after** the
`loss_weights = "auto"` baseline scaling, so `prior_lambda` means the
same thing on every dataset.

[`createProcess(prior = "max_obs_recruits")`](https://finnverse.github.io/FINN/reference/createProcess.md),
the default, resolves the prior from the data at the start of
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md): the
largest observed recruitment rate per timestep,
`max(reg / period_length)`. The division matters on any inventory with a
remeasurement interval, because `reg` is the SUM over that interval
while the cap bounds the mean of one timestep – taking the column as
supplied would place the prior a factor `period_length` too high, and
the cap would never bind. A number sets the prior directly, on that same
per-timestep scale. `prior_lambda = 0` switches the prior off and leaves
the cap to the data alone. Scoring no recruitment at all is an error
rather than a silent default, because nothing else in the objective
constrains a cap.

The resolved value is also the cap's starting value, and it is resolved
only on the FIRST
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) of a model:
a second call continues from the trained cap instead of resetting it to
the prior.

With `cap = "species"` the prior is **hierarchical**: each species' log
cap is pulled towards the mean of all current log caps rather than
towards a fixed value, so the prior penalises how far the species spread
apart and leaves the overall level to the data. Every species still
starts at the observed maximum.

## See also

[regeneration](https://finnverse.github.io/FINN/reference/regeneration.md),
[`regeneration_adult()`](https://finnverse.github.io/FINN/reference/regeneration_adult.md),
[`regeneration_adult_effect()`](https://finnverse.github.io/FINN/reference/regeneration_adult_effect.md),
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md),
[recruitment](https://finnverse.github.io/FINN/reference/recruitment.md)
for the recruitment options and how to combine them

## Examples

``` r
library(data.table)
Nsp <- max(fia_obs_dt$species)
# the default: one cap for all species, prior from the data
m <- finn(N_species = Nsp, recruits_dbh = 12.9)

# one cap per species, a fixed prior (stems per ha and year), a stronger pull
m2 <- finn(N_species = Nsp, recruits_dbh = 12.9,
           regeneration_process = createProcess(~ temp + prec, regeneration_capped,
                                                cap = "species", prior = 500, prior_lambda = 1))

# the cap is created and resolved by fit(); read it back on the natural scale
sites <- 1:3
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
fit(m2, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 1L, patch_size = 0.06, plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.768)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%
exp(as.numeric(m2$reg_logK))   # one carrying capacity per species
#>  [1] 495.0253 495.0258 495.0256 495.0253 505.0248 495.0256 495.0253 495.0256
#>  [9] 495.0256 495.0253 495.0256
```
