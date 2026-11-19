# Plot the fitted parameters of a FINN model

Draws the fitted parameters themselves, per species, on the scale on
which they enter the model – no simulation is involved. By default the
process parameters; `pars = "env"` switches to the environmental niches.

## Usage

``` r
# S3 method for class 'finn_class'
plot(
  x,
  pars = c("process", "env"),
  species = NULL,
  env = NULL,
  species_names = NULL,
  parameter_names = NULL,
  ...
)
```

## Arguments

- x:

  A model built by
  [`finn()`](https://finnverse.github.io/FINN/reference/finn.md),
  usually fitted with
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md).

- pars:

  Which parameters: `"process"` (default) the species-specific process
  parameters, `"env"` the environmental niches. The argument and its two
  values are those of the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method of
  FINN 0.1.0, so a call such as `plot(m, pars = "env")` keeps working.

- species:

  Species numbers to draw; `NULL` (default) all.

- env:

  Table whose predictor ranges the niches span, in the format of the
  `env` of [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)
  (see
  [finn_data](https://finnverse.github.io/FINN/reference/finn_data.md));
  `NULL` (default) the environment the model was fitted on. Only used
  with `pars = "env"`.

- species_names:

  Species names for the labels: a character vector in species order, a
  vector named by species number, or a table with `species` and
  `species_name`. `NULL` (default): the names
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) kept from
  the `species_name` column of the data, else "species 1", "species 2",
  ...

- parameter_names:

  Names of the species-specific parameters of a process function of your
  own, as a list by process, e.g.
  `list(mortality = c("shade", "size", "growth"))`. `NULL` (default):
  the documented names for FINN's own process functions
  ([`growth()`](https://finnverse.github.io/FINN/reference/growth.md),
  [`mortality()`](https://finnverse.github.io/FINN/reference/mortality.md),
  [`regeneration()`](https://finnverse.github.io/FINN/reference/regeneration.md),
  [`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md),
  [`competition()`](https://finnverse.github.io/FINN/reference/competition.md)),
  else "parameter 1", "parameter 2", ...

- ...:

  Unused.

## Value

The `ggplot` object, invisibly (`NULL` when there is nothing to draw);
it is printed.

## Details

- **Process parameters** (`pars = "process"`, the default): one panel
  per parameter of the growth, mortality, regeneration and competition
  processes, and the recruit dispersion `theta`; one point per species,
  on the constrained scale the process function uses, the axis spanning
  the fitted values. A bound of the admissible range set by the
  arguments `lower` and `upper` of
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  (for `theta` the floor, its argument `theta_floor`) is drawn as a
  dashed line where a value comes near it, and a species whose parameter
  ended the fit at a bound is an open triangle (see
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  on what such a parameter means).

- **Environmental niches** (`pars = "env"`): for each process whose
  environmental response is linear (no `hidden`, no custom `NN`, not
  [`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md)),
  the fitted environmental term of every species across the range each
  predictor covers in `env`, with the other predictors held at their
  mean. The curve is the linear predictor the process function receives
  (`pred`): for FINN's default functions growth and regeneration use
  `exp(pred)` and mortality `plogis(pred + ...)`, so a curve higher up
  means faster growth, more recruits or more deaths. The x axis is on
  the original scale of the predictor; the stored standardisation of
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) is
  applied internally. A process with a network response is skipped with
  a warning: it has no coefficients to show, and
  [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) is the
  way to look at it. Competition does not depend on the environment.

## See also

[`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) for the
effect of a network response,
[`summary.finn_class()`](https://finnverse.github.io/FINN/reference/summary.finn_class.md)
for the skill, the convergence and the environmental importances of a
fit,
[`plot.finn_simulation()`](https://finnverse.github.io/FINN/reference/plot.finn_simulation.md)
for a simulation.

## Examples

``` r
library(data.table)
sites <- 1:3
Nsp   <- max(fia_obs_dt$species)
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
FINN.seed(1)
m <- finn(N_species = Nsp, recruits_dbh = 12.9,
          growth_process       = createProcess(~ temp + prec, growth),
          mortality_process    = createProcess(~ temp + prec, mortality),
          regeneration_process = createProcess(~ temp + prec, regeneration_capped))
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
    plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.768)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%
plot(m)

plot(m, pars = "env", species = 1:3)
```
