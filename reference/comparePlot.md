# Compare a simulation with observations

Pairs a simulation with the observed data and draws, per response, how
well they agree. The simulation is aggregated the way
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) scores it:
the stand states `dbh`, `ba` and `trees` at the census year, `growth`
and `mort` compounded over the `period_length` years before the census
(the geometric mean of `1 + g` minus one; one minus the geometric-mean
survival), and the expected recruits (`r_mean_ha`) summed over them for
`reg` (see
[finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md)
and, for the observed columns,
[finn_data](https://finnverse.github.io/FINN/reference/finn_data.md)).
An observation is placed as in
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md): its `year`
is a year of the `env` the simulation was run with, and the simulated
step is the position of that year among the years of `env` (for a plain
site table, which does not carry those years, `year` is taken as the
step itself). Observations of a site or a species that is not in the
simulation, or of a year that is not simulated, are left out with a
message. Both tables may come in long (`variable`, `value`) or wide
form.

## Usage

``` r
comparePlot(
  sim,
  obs,
  type = c("scatter", "timeseries", "heatmap"),
  species = NULL,
  response = NULL,
  sites = NULL,
  tree_data = NULL,
  n_species = 12L,
  species_names = NULL,
  ...
)

# S3 method for class 'finn_heatmap'
print(x, ...)

# S3 method for class 'finn_comparison'
print(x, ...)
```

## Arguments

- sim:

  The result of
  [`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  or
  [`predict()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  for the sites of `obs`, started from the same initial state (or its
  site table).

- obs:

  Observations as
  [`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
  returns them (wide), or long with `variable` and `value`;
  `period_length`, `n_at_risk` and `n_died` are used where present.

- type:

  `"scatter"` (default), `"heatmap"` or `"timeseries"`.

- species:

  `NULL` (default) pools all species in one panel per response; species
  numbers draw one row of panels per species.

- response:

  Responses to draw, any of `"dbh"`, `"ba"`, `"trees"`, `"growth"`,
  `"mort"`, `"reg"` (and `"dbh_max"` with `tree_data`); `NULL` (default)
  all that are observed.

- sites:

  Site ids to include; `NULL` (default) all.

- tree_data:

  `NULL` (default), or the observed trees in the format of the
  `tree_data` of
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) (see
  [finn_data](https://finnverse.github.io/FINN/reference/finn_data.md):
  `siteID`, `year`, `species`, `dbh` at the start of the interval, and
  `growth`, `died`, `period_length`; the interval is taken from `obs`
  when the trees do not carry it): adds the largest diameter (`dbh_max`)
  to the comparison.

- n_species:

  Number of species in the recruitment panel. Default `12`.

- species_names:

  Species names for the labels: a character vector in species order, a
  vector named by species number, or a table with `species` and
  `species_name`. `NULL` (default): the names the simulation carries
  from the model, else the `species_name` column of `obs`.

- ...:

  Unused.

- x:

  A `finn_comparison` object (the scatter view of `comparePlot()`).

## Value

`type = "scatter"`: a `finn_comparison` object (a list of the `ggplot`
panels `scatter`, `mort` and `reg`, which
[`print()`](https://rdrr.io/r/base/print.html) arranges on one page),
invisibly; `type = "heatmap"`: a `finn_heatmap` object (the panels
`skill` and `ratio`, side by side when printed), invisibly;
`type = "timeseries"`: the `ggplot` object, invisibly. It is printed.
The paired table is attached as attribute `"matched"` and the skill
measures as `"skill"` (per response, and per species when `species` is
given; the measures are described in
[`summary.finn_class()`](https://finnverse.github.io/FINN/reference/summary.finn_class.md)).

## Details

- `type = "scatter"` (default): observed against simulated for `dbh`,
  `ba`, `trees` and `growth`, with the 1:1 line and R2 (`1 - SSE/SST`)
  in each panel; for mortality the ROC curve of the interval death
  probability with its AUC (each tree at risk is one case, `n_died`
  positives and `n_at_risk - n_died` negatives; without those two
  columns in `obs`, a scatter of the annual rates); for recruitment,
  whose site-level values are too noisy for a correlation, the ratio of
  simulated to observed recruitment summed over all sites, per species
  (the `n_species` species with most observed recruitment). With
  `tree_data`, the dbh panel shows the largest diameter instead of the
  mean: per site, census and species, or, all species pooled, per site
  and census (the maximum over the species), which is robust against a
  few large or many small trees. The simulation must have recorded its
  cohorts at the census steps
  (`simulateForest(..., return_cohorts = <steps>)`, a step being the
  position of the census year among the years of `env`;
  [`evaluate()`](https://finnverse.github.io/FINN/reference/evaluate.md)
  does this itself).

- `type = "heatmap"`: the skill per species at a glance, every species a
  row (by name): R2 of every response and the AUC of mortality on the
  left, the ratio of simulated to observed on the right (log scale,
  white at 1, so over- and underestimation show as colours).

- `type = "timeseries"`: the simulated trajectory over the years as a
  line and the observations as points at their census years, averaged
  over the selected sites. With all species pooled, `ba` and `trees` are
  the stand totals and the other responses the mean over the species.

`plot(sim, obs)` is the same as `comparePlot(sim, obs)` (see
[`plot.finn_simulation()`](https://finnverse.github.io/FINN/reference/plot.finn_simulation.md)).

## See also

[`evaluate()`](https://finnverse.github.io/FINN/reference/evaluate.md),
which returns the same comparison and keeps it on the model, and
[`summary.finn_class()`](https://finnverse.github.io/FINN/reference/summary.finn_class.md),
which prints it;
[`plot.finn_simulation()`](https://finnverse.github.io/FINN/reference/plot.finn_simulation.md)
for the simulation alone.

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
comparePlot(sim, fia_obs_dt[siteID %in% sites])

plot(sim, fia_obs_dt[siteID %in% sites], type = "timeseries", species = 1:2)
```
