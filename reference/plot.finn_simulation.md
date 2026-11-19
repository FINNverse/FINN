# Plot a FINN simulation

Draws the trajectories returned by
[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
or
[`predict()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md):
one panel per response, one line per species, averaged over the sites
(dbh over the sites where the species has trees); with `sites` one site
(or the average over a few). To see every site, loop over them. The
`n_species` leading species are coloured and named, every other species
is a pale grey line, and the stand total of `ba` and `trees` is dashed.

## Usage

``` r
# S3 method for class 'finn_simulation'
plot(
  x,
  y = NULL,
  response = NULL,
  sites = NULL,
  ...,
  n_species = 6L,
  species_names = NULL
)

# S3 method for class 'finn_simulation'
print(x, ...)
```

## Arguments

- x:

  The object returned by
  [`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  or
  [`predict()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  (class `finn_simulation`).

- y:

  `NULL` (default), or observations to compare with (the `obs` of
  [`comparePlot()`](https://finnverse.github.io/FINN/reference/comparePlot.md)).

- response:

  Responses to draw, any of `"dbh"`, `"ba"`, `"trees"`, `"growth"`,
  `"mort"`, `"reg"`, `"r_mean_ha"`; `NULL` (default) the stand variables
  `ba`, `dbh` and `trees`. With `y`, as in
  [`comparePlot()`](https://finnverse.github.io/FINN/reference/comparePlot.md).

- sites:

  Site ids to include; `NULL` (default) all.

- ...:

  [`plot()`](https://rdrr.io/r/graphics/plot.default.html): passed to
  [`comparePlot()`](https://finnverse.github.io/FINN/reference/comparePlot.md)
  when `y` is given (`type`, `species`), otherwise unused.
  [`print()`](https://rdrr.io/r/base/print.html): passed to the print
  method of a list.

- n_species:

  Number of species drawn in colour: those reaching the highest value of
  the first response; every other species is a pale grey line. Default
  `6`. With `y`, the `n_species` of
  [`comparePlot()`](https://finnverse.github.io/FINN/reference/comparePlot.md)
  (its default when not given).

- species_names:

  Species names for the legend (see
  [`plot.finn_class()`](https://finnverse.github.io/FINN/reference/plot.finn_class.md));
  `NULL` (default): the names the simulation carries from the model.

## Value

[`plot()`](https://rdrr.io/r/graphics/plot.default.html): the `ggplot`
object, invisibly; it is printed.
[`print()`](https://rdrr.io/r/base/print.html): `x`, invisibly.

## Details

A bare-ground run (`simulateForest(bare_ground = TRUE)`) is drawn as a
succession: basal area by default, one panel per site, labelled with the
site's mean climate.

Given observations as the second argument, `plot(sim, obs)` compares the
simulation with them instead: it is
[`comparePlot()`](https://finnverse.github.io/FINN/reference/comparePlot.md),
with the same arguments.

[`print()`](https://rdrr.io/r/base/print.html) of a simulation shows the
list of tables as before 0.2.0, without the class attribute.

## See also

[`plot.finn_class()`](https://finnverse.github.io/FINN/reference/plot.finn_class.md)
for the fitted parameters,
[`comparePlot()`](https://finnverse.github.io/FINN/reference/comparePlot.md)
and
[`evaluate()`](https://finnverse.github.io/FINN/reference/evaluate.md)
for the comparison with observations.

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
sim <- simulateForest(m, env = fia_env_dt[siteID %in% sites], init_cohort = init)
plot(sim)

plot(sim, response = c("growth", "mort"), sites = 1)

# against the observations
plot(sim, fia_obs_dt[siteID %in% sites])
```
