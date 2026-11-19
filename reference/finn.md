# Build a FINN forest model

Builds a FINN (Forest Informed Neural Network) model: a forest simulator
that follows groups of similar trees (cohorts) on small patches of
forest through time, year by year, and that can learn its parameters
from forest inventory data. Each of the four demographic processes –
competition for light, growth, mortality and regeneration – is either a
classic, interpretable equation
([`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md))
or a neural network
([`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md)).
Whatever the mix, all processes are estimated together, in one go, by
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md).

`finn()` only assembles the model. Nothing is estimated yet: pass the
result to [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)
to estimate the parameters from your data, and to
[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
to simulate forest development.

## Usage

``` r
finn(
  N_species,
  mortality_process = createProcess(func = mortality),
  growth_process = createProcess(func = growth),
  regeneration_process = createProcess(func = regeneration_capped),
  competition_process = createProcess(func = competition),
  recruits_dbh = 1,
  growth_rate = c("annual", "period"),
  recruit_obs_weight = 1,
  ...
)
```

## Arguments

- N_species:

  (`integer(1)`)  
  Number of species in the data. Species are coded `1, ..., N_species`
  in every table
  ([`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
  does the coding).

- mortality_process:

  (`process`)  
  How trees die. Default `createProcess(func = mortality)`, the
  mortality equation of
  [`mortality()`](https://finnverse.github.io/FINN/reference/mortality.md)
  with a logistic response to every environment column. Build your own
  with
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  or
  [`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md).
  Its environmental response starts at about 2 % mortality per time step
  (an intercept of `qlogis(0.02)`, every other weight 0), so that an
  unfitted stand on yearly steps does not die out;
  `createProcess(initEnv = )` sets another start, see
  [`mortality()`](https://finnverse.github.io/FINN/reference/mortality.md).

- growth_process:

  (`process`)  
  How trees grow in diameter. Default `createProcess(func = growth)`,
  the growth equation of
  [`growth()`](https://finnverse.github.io/FINN/reference/growth.md)
  with a linear response to every environment column. Build your own
  with
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  or
  [`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md).

- regeneration_process:

  (`process`)  
  How new trees appear. Default
  `createProcess(func = regeneration_capped)`: the regeneration equation
  of
  [`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md),
  whose expected recruitment is capped (new in 0.2.0). Alternatives:
  [`regeneration()`](https://finnverse.github.io/FINN/reference/regeneration.md)
  (no cap, the equation of the paper),
  [`regeneration_adult()`](https://finnverse.github.io/FINN/reference/regeneration_adult.md)
  (recruitment limited by the adults of the same species on the site)
  and
  [`regeneration_adult_effect()`](https://finnverse.github.io/FINN/reference/regeneration_adult_effect.md)
  (adults raise recruitment, but a species without adults still
  recruits). See
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md)
  for the recruitment options and how they act together. Build your own
  with
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  or
  [`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md).

- competition_process:

  (`process`)  
  How trees shade each other. Default
  `createProcess(func = competition)`, the light model of
  [`competition()`](https://finnverse.github.io/FINN/reference/competition.md).
  Only
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  is accepted here, and competition cannot depend on the environment in
  this version.

- recruits_dbh:

  (`numeric(1)`)  
  Diameter (cm) at which new trees enter the simulation. Default `1`.
  Set it to the smallest diameter your inventory records (12.9 cm for
  the bundled FIA data), so that simulated recruits are the trees the
  inventory counts as recruits.

- growth_rate:

  (`character(1)`)  
  What the `growth` column of the data is, so that
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) compares
  it with the right simulated quantity. `"annual"` (default): a rate per
  time step – with annual steps the mean annual relative growth rate
  over the observation interval (each year's `dbh / dbh_before - 1`,
  compounded to one rate per year). FINN compounds its simulated
  per-step rates over the interval the same way, so the two are the same
  quantity (see
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)).
  `"period"`: the relative diameter increment over the whole observation
  interval (`dbh_end / dbh_start - 1`), which compounds the per-step
  rates; FINN then compounds its simulated rates the same way. The two
  only differ when an observation interval spans several time steps
  (`period_length` in the data). If unsure, keep the default: with
  annual time steps it is the annual growth rate.
  [`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
  returns the annual rate (as in the bundled FIA data); with
  [`makeObsData(rates = "interval")`](https://finnverse.github.io/FINN/reference/makeObsData.md)
  it returns the increment between two censuses, which is the per-step
  rate only when one step is one census interval. The terms follow
  Kohyama, Kohyama & Sheil (2018), who distinguish a rate per year from
  the change over a census interval.

- recruit_obs_weight:

  (`numeric(1)`)  
  Only for inventories that do not count recruits per hectare
  (angle-count sampling, small-tree subplots): the factor that turns a
  recruitment density into the counts such a design tallies. Default
  `1`, which leaves per-hectare data as they are. See
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md).

- ...:

  Arguments of earlier versions, still accepted with a warning:
  `reg_floor` (now
  [`createProcess(reg_floor = )`](https://finnverse.github.io/FINN/reference/createProcess.md)
  of the regeneration process) and `growth_period_scale` (`TRUE` is now
  `growth_rate = "period"`).

## Value

An object of class `finn_class` (a
[torch::nn_module](https://torch.mlverse.org/docs/reference/nn_module.html)):
the assembled, not yet fitted model, ready for
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) or
[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md).

## Details

The defaults work for a first model: `finn(N_species = 11)` builds a
model whose four processes are the package's own equations, with every
column of the environment table as a predictor. Pass a process built
with
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
to choose the predictors (its `formula`), to change the equation (its
`func`), or to fix parameters; pass one built with
[`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md)
to replace the whole process by a neural network.

![One time step of FINN: the cohorts of a patch and the site's
environment enter competition, growth, mortality and regeneration, which
produce the cohorts of the next step and the simulated stand
variables.](figures/FINN-overview.jpg)

*One simulation step.* The cohorts of a patch at time t – their size,
species, number of trees and cohort id – and the site's environment are
the inputs. Competition turns the cohorts into the light each of them
receives; growth, mortality and regeneration then act on every cohort,
each of them either a mechanistic equation (the gear) or, in a hybrid
model, a neural network, and each with its own response to the
environment. The distribution drawn beside a process is the one its
outcome is drawn from or scored with: Gaussian for growth, binomial for
the deaths, negative binomial for the recruits. Out come the cohorts of
the next step and the quantities FINN reports for every site and year –
mean diameter, basal area, number of trees, and the growth, mortality
and recruitment rates.

The figure is Figure 1 of Pichler & Käber (2026), reproduced from the
preprint version of that paper (arXiv:2508.01228), which is published
under CC BY 4.0. The caption above is ours, not the paper's.

**How FINN simulates a forest.** A site is a set of patches (small plots
of `patch_size` ha). On each patch the trees are kept as cohorts, each
with a species, a diameter at breast height (dbh, cm) and a number of
trees. Every simulated time step, usually one year, runs the four
processes in turn:

1.  *competition* works out how much light reaches each cohort from the
    basal area of the taller trees on its patch
    ([`competition()`](https://finnverse.github.io/FINN/reference/competition.md));

2.  *growth* increases each cohort's diameter, depending on its light,
    size, species and the site's environment
    ([`growth()`](https://finnverse.github.io/FINN/reference/growth.md));

3.  *mortality* removes a share of each cohort's trees, depending on its
    growth, light, size, species and environment
    ([`mortality()`](https://finnverse.github.io/FINN/reference/mortality.md));

4.  *regeneration* adds new cohorts of small trees (`recruits_dbh`),
    depending on light, species and environment
    ([`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md)).

Growth, mortality and regeneration each carry their own response to the
environment: for every species, the columns of the environment table
that the process' `formula` names are turned into one number that scales
the process. By default that response is a simple linear (growth,
regeneration) or logistic (mortality) function with one coefficient per
predictor and species;
[`createProcess(hidden = )`](https://finnverse.github.io/FINN/reference/createProcess.md)
turns it into a small neural network. The simulated stand variables
(mean dbh, basal area, number of trees) and the simulated rates (growth,
mortality, recruitment) of every site and year are what
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) compares
with the observations.

**The default processes** are, in 0.2.0,
[`createProcess(func = growth)`](https://finnverse.github.io/FINN/reference/createProcess.md),
`createProcess(func = mortality)`,
`createProcess(func = regeneration_capped)` and
`createProcess(func = competition)`: the equations of Pichler & Käber
(2026), with one change – recruitment is capped
([`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md)
is new in 0.2.0;
[`regeneration()`](https://finnverse.github.io/FINN/reference/regeneration.md)
is the uncapped process of the paper). A process written in the
signature as `createProcess(func = )` uses every column of the
environment table (`formula = ~.`) and its default settings. Passing
`NULL` for a process is still accepted and means the same default. The
competition process cannot depend on the environment in this version.

**Patches, and why there are several.** FINN is a gap model: how much
light a tree gets depends on the trees standing next to it, so the unit
of simulation is a patch small enough that one large tree shades it
(`patch_size`, in ha). A site is simulated as several such patches that
share the site's environment, and a site's stand variables are the mean
over its patches. The patches mean two different things depending on
where the simulation starts. With an `init_cohort` built from an
inventory they are the inventory's own subplots, each carrying its real
trees, and the patch count comes from that table rather than from
`patches`. From bare ground (`init_cohort = NULL`) they are stochastic
replicates of the same site: more of them means a smoother trajectory,
at a cost that grows with their number.

**What is estimated, and how.** Everything that is switched on, at once.
A process built with
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
estimates the weights of its environmental response by default
(`optimizeEnv = TRUE`), and its per-species parameters only when asked
to (`optimizeSpecies = TRUE`; by default they keep their starting
values). These, and where a process is a network that network's weights,
form one parameter vector, and
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) moves it by
gradient descent *through the simulation*: the model is run forward over
the observed years, the simulated stand variables and rates are compared
with the observations, and the gradient is carried back through the time
loop by automatic differentiation. No process is fitted on its own and
plugged in afterwards. That is what lets a process that has been
replaced by a network still be identified: it has to produce a stand
that the other three processes can also live with.

**Mechanistic, hybrid, or in between.** Each process is set
independently:

- *fully mechanistic* – the equation, with a linear or logistic response
  to the environment. Every parameter has a meaning and a unit.

- *level 1* – the equation, but its environmental response is a small
  network
  ([`createProcess(hidden = c(20, 20))`](https://finnverse.github.io/FINN/reference/createProcess.md)).
  The demography stays the package's equation, so the model still
  extrapolates the way that equation says it should, while the response
  to climate is free in shape.

- *level 2* – the equation itself replaced by a network
  ([`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md)).
  The data get the most freedom and the fitted process the least
  interpretation.

Competition is mechanistic at every level. A mixed model is the normal
case: the process you suspect is misspecified goes to level 1 or 2, the
others stay mechanistic, and the comparison between the two fits is the
result.

**Setting up a first model.** Three choices carry most of the outcome.
`recruits_dbh` must be the diameter at which your inventory starts
counting trees, or the simulation recruits trees the data never
recorded. `patch_size` must be the area one plot of that inventory
covers. And a process' `formula` should name the predictors you would
defend in a regression rather than every column of `env`, because each
one costs a coefficient per species.

**Checking a fitted model.**
[`fit(plot_progress = TRUE)`](https://finnverse.github.io/FINN/reference/fit.md)
and `$history` say whether the fit moved and which response drove it;
[`summary.finn_class()`](https://finnverse.github.io/FINN/reference/summary.finn_class.md)
reports the skill on the training data and a holdout (stored by
[`evaluate()`](https://finnverse.github.io/FINN/reference/evaluate.md),
drawn by
[`comparePlot()`](https://finnverse.github.io/FINN/reference/comparePlot.md)),
whether the fit has settled, and the environmental importances;
[`plot.finn_class()`](https://finnverse.github.io/FINN/reference/plot.finn_class.md)
draws the fitted parameters;
[`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) shows the
fitted environmental responses, including those of a network;
[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
from bare ground (`bare_ground = TRUE`) shows what the model does
outside the observation window. The last one is worth the time: a model
that predicts the holdout well and grows an implausible stand from bare
ground has usually learned the observed window rather than the dynamics
that produced it.

## References

Pichler, M., & Käber, Y. (2026). Inferring processes within dynamic
forest models using hybrid modelling. *Methods in Ecology and
Evolution*.
[doi:10.1111/2041-210x.70347](https://doi.org/10.1111/2041-210x.70347)

Kohyama, T. S., Kohyama, T. I., & Sheil, D. (2018). Definition and
estimation of vital rates from repeated censuses: choices, comparisons
and bias corrections focusing on trees. *Methods in Ecology and
Evolution*, 9, 809-821.
[doi:10.1111/2041-210X.12929](https://doi.org/10.1111/2041-210X.12929)

## See also

[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) to estimate
the parameters,
[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
to simulate,
[finn_data](https://finnverse.github.io/FINN/reference/finn_data.md) for
the input tables,
[recruitment](https://finnverse.github.io/FINN/reference/recruitment.md)
for the recruitment options

## Examples

``` r
library(data.table)

# --- the bundled data -------------------------------------------------
# 200 plots of the US Forest Inventory and Analysis programme in Oregon,
# 11 tree species, four subplots per plot, each plot measured twice about
# ten years apart. `year` counts years from the start of the panel.
#  * fia_obs_dt     one row per plot, year and species: the stand variables
#                   (dbh, ba, trees) and the rates (growth, mort, reg) that
#                   an inventory yields -- these are what fit() scores
#  * fia_env_dt     climate of every plot and year (temp, prec, ...)
#  * fia_init_trees the individual trees of the first census, the state the
#                   simulation starts from
#  * fia_species_dt the species codes 1..11 and their names
head(fia_obs_dt)
#>             species_name siteID  year        ba trees      reg      dbh
#>                   <char>  <int> <int>     <num> <num>    <num>    <num>
#> 1: Pseudotsuga menziesii      1    10 0.3416488   0.5  0.00000 92.83700
#> 2:       Pinus ponderosa      1    10 0.0000000   0.0  0.00000       NA
#> 3:        Pinus contorta      1    10 0.0000000   0.0  0.00000       NA
#> 4:        Abies concolor      1    10 0.0000000   0.0  0.00000       NA
#> 5:         Abies grandis      1    10 0.4447334   2.5 16.66667 45.98458
#> 6:    Tsuga heterophylla      1    10 0.0000000   0.0  0.00000       NA
#>         growth growth_n n_at_risk n_died  mort period_length species
#>          <num>    <int>     <int>  <int> <num>         <int>   <int>
#> 1: 0.003564897        2         2      0     0            10       1
#> 2:          NA        0         0      0    NA            10       2
#> 3:          NA        0         0      0    NA            10       3
#> 4:          NA        0         0      0    NA            10       4
#> 5: 0.028781159        6         6      0     0            10       5
#> 6:          NA        0         0      0    NA            10       6
fia_species_dt
#>     species            species_name
#>       <int>                  <char>
#>  1:       1   Pseudotsuga menziesii
#>  2:       2         Pinus ponderosa
#>  3:       3          Pinus contorta
#>  4:       4          Abies concolor
#>  5:       5           Abies grandis
#>  6:       6      Tsuga heterophylla
#>  7:       7       Tsuga mertensiana
#>  8:       8 Lithocarpus densiflorus
#>  9:       9             Alnus rubra
#> 10:      10          Abies amabilis
#> 11:      11                   other

op <- par(mfrow = c(1, 2))
# what the model has to reproduce: the stand of each plot, per year
stand <- fia_obs_dt[, .(ba = sum(ba, na.rm = TRUE)), by = .(siteID, year)]
plot(stand$year, stand$ba, pch = 16, col = "#2b6a4d40",
     xlab = "year of the panel", ylab = "basal area (m^2/ha)",
     main = "observed stands")
# and the size distribution it starts from
hist(fia_init_trees$dbh[which(fia_init_trees$living)], breaks = 40, col = "grey80",
     border = "white", xlab = "dbh at the first census (cm)",
     main = "initial trees")

par(op)

# --- the workflow -----------------------------------------------------
# The whole workflow on three Oregon FIA sites: build, estimate, simulate.
# A real calibration uses all sites and hundreds of epochs.
sites <- 1:3
Nsp   <- max(fia_obs_dt$species)

# 1. build the model: growth, mortality and regeneration respond to annual
#    temperature and precipitation; competition has no environment
m <- finn(N_species = Nsp, recruits_dbh = 12.9,
          growth_process       = createProcess(~ temp + prec, growth),
          mortality_process    = createProcess(~ temp + prec, mortality),
          regeneration_process = createProcess(~ temp + prec, regeneration_capped))
m
#> An `nn_module` containing 12 parameters.
#> 
#> ── Parameters ──────────────────────────────────────────────────────────────────
#> • par_theta_recruits_raw: Float [1:11]
#> • reg_logK: Float [1:1]
#> 
#> ── Buffers ─────────────────────────────────────────────────────────────────────
#> • par_mortality_upper: Float [1:3]
#> • par_mortality_lower: Float [1:3]
#> • par_mortality_unconstrained: Float [1:11, 1:3]
#> • par_growth_upper: Float [1:2]
#> • par_growth_lower: Float [1:2]
#> • par_growth_unconstrained: Float [1:11, 1:2]
#> • par_regeneration_upper: Float [1:1]
#> • par_regeneration_lower: Float [1:1]
#> • par_regeneration_unconstrained: Float [1:11, 1:1]
#> • par_competition_upper: Float [1:2]
#> • par_competition_lower: Float [1:2]
#> • par_competition_unconstrained: Float [1:11, 1:2]

# 2. estimate the parameters from the inventory: the observed stand
#    variables and rates (data), the climate of every year (env) and the
#    trees at the first census (init_cohort)
init <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
FINN.seed(1)
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
    plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.768)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%
tail(m$history, 1)   # the loss of every response in the last epoch
#> [[1]]
#> [1]  2.080208  0.541941  0.478983 24.514814  1.374413  1.427557  0.000000
#> 

# 3. simulate the three sites over the years of `env` with the fitted model
sim <- simulateForest(m, env = fia_env_dt[siteID %in% sites],
                      init_cohort = init, patches = 4, patch_size = 0.06)
head(sim$long$site)
#>    siteID  year species variable    value
#>     <int> <int>   <int>   <fctr>    <num>
#> 1:      1     1       1      dbh 81.00741
#> 2:      1     2       1      dbh 82.33847
#> 3:      1     3       1      dbh 83.58874
#> 4:      1     4       1      dbh 84.76591
#> 5:      1     5       1      dbh 85.87665
#> 6:      1     6       1      dbh 86.92680

# the same model with the defaults spelled out: every process uses every
# column of `env` as a predictor
m2 <- finn(N_species = Nsp, recruits_dbh = 12.9,
           growth_process       = createProcess(func = growth),
           mortality_process    = createProcess(func = mortality),
           regeneration_process = createProcess(func = regeneration_capped),
           competition_process  = createProcess(func = competition))
```
