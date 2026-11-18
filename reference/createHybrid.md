# Define a hybrid (deep-neural-network) demographic process for FINN

Configures a process (growth, mortality, or regeneration) in which the
entire process equation – not just its environmental-response function –
is replaced by a deep neural network (DNN), for use as `growth_process`,
`mortality_process`, or `regeneration_process` in
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md). This is
the second ("Level 2") level of hybridization described by Pichler &
Käber (2026); to replace only the environmental-response function while
keeping the rest of the process mechanistic ("Level 1"), use
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
instead. `competition_process` must always be created with
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md),
as the competition (light availability) process does not support full
replacement by a DNN.

## Usage

``` r
createHybrid(
  formula = NULL,
  optimize = TRUE,
  dispersion_parameter = 1,
  NN = NULL,
  dropout = 0.3,
  encoder_layers = 1L,
  hidden = c(50L, 50L),
  sample_regeneration = TRUE,
  transformer = TRUE,
  emb_dim = 20L,
  dim_feedforward = 256L,
  mortality_draw = c("binomial", "bernoulli", "poisson", "gamma"),
  regeneration_sampler = c("refined", "gaussian"),
  cap = c("all", "species"),
  prior = "max_obs_recruits",
  prior_lambda = 0.1,
  theta_floor = 0.5,
  optimizeTheta = TRUE,
  mortality_floor = 0.001,
  mortality_ceiling = 0.999
)
```

## Arguments

- formula:

  (`formula`)  
  Environmental predictors passed to the network, evaluated against the
  `env` data. Default `NULL` uses `~.` (all available environmental
  covariates).

- optimize:

  (`logical(1)`)  
  Should the network's weights be estimated during
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)? Default
  `TRUE`.

- dispersion_parameter:

  (`numeric(1)`)  
  Initial dispersion parameter of the negative binomial distribution
  used to sample recruits, or its fixed value with
  `optimizeTheta = FALSE`. Only used when this object is the
  regeneration process. See
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md).

- NN:

  (`nn_module`)  
  Not supported in 0.2.0. A custom module was accepted and then ignored,
  so a non-`NULL` value now raises an error; shape the default network
  with `hidden`, `dropout` and `transformer` instead.

- dropout:

  (`numeric(1)`)  
  Dropout rate of the network. Default `0.3`; the Barro Colorado Island
  case study of Pichler & Käber (2026) used `0.1`.

- encoder_layers:

  (`integer(1)`)  
  Number of transformer encoder layers. Only used when
  `transformer = TRUE`.

- hidden:

  ([`integer()`](https://rdrr.io/r/base/integer.html))  
  Hidden-layer sizes of the feed-forward network. Only used when
  `transformer = FALSE`.

- sample_regeneration:

  (`logical(1)`)  
  Should recruits actually be sampled from the negative binomial
  regeneration distribution (`TRUE`), or only its expectation used
  (`FALSE`)? Only used when this object is the regeneration process.
  With `FALSE` the expected number of recruits per patch is rounded, so
  wherever it is below 0.5 no recruit enters at all. See
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md).

- transformer:

  (`logical(1)`)  
  Use a transformer encoder architecture (`TRUE`, default) or a
  feed-forward network (`FALSE`).

- emb_dim:

  (`integer(1)`)  
  Embedding dimension for the species/cohort/ environment features.

- dim_feedforward:

  (`integer(1)`)  
  Dimension of the feed-forward sublayer within each transformer encoder
  layer. Only used when `transformer = TRUE`.

- mortality_draw:

  (`character(1)`)  
  How the number of dead trees is drawn, `"binomial"` (default since
  0.2.0), `"bernoulli"`, `"poisson"` or `"gamma"`; see
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md).
  Only used when this object is the mortality process.

- regeneration_sampler:

  (`character(1)`)  
  How recruits are drawn, `"refined"` (default since 0.2.0) or
  `"gaussian"`; see
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md).
  Only used when this object is the regeneration process.

- cap, prior, prior_lambda:

  Accepted for symmetry with
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md),
  but not used in 0.2.0: the hybrid regeneration network has no
  Beverton-Holt cap. See
  [`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md)
  for the cap of the mechanistic process.

- theta_floor:

  (`numeric(1)`)  
  Lower bound on the recruit dispersion `par_theta_recruits`, default
  `0.5` since 0.2.0, `0` or `NULL` to leave it unbounded as in 0.1.0;
  see
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md).
  Only used when this object is the regeneration process. The floor
  keeps the fit stable. In one converged pair of fits on the German
  national forest inventory it calibrated the recruitment level as well
  as a fixed `theta`, but only after the parameters had settled; read
  the level once `$param_history` shows them settled, or fix `theta`
  (`optimizeTheta = FALSE`) for short fits. See
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md).

- optimizeTheta:

  (`logical(1)`)  
  Should the recruit dispersion be estimated during
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)? Default
  `TRUE`; `FALSE` keeps it at `dispersion_parameter`, which is the
  choice when the recruitment level is to be calibrated; see
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  and
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md).
  Only used when this object is the regeneration process.

- mortality_floor, mortality_ceiling:

  (`numeric(1)`)  
  The interval the mortality probability is clamped to per simulation
  step before the death draw, defaults `0.001` and `0.999`; see
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md).
  Only used when this object is the mortality process.

## Value

A list of class `"hybrid"` containing the process definition and
associated parameters, to be passed as `growth_process`,
`mortality_process`, or `regeneration_process` to
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md).

## Details

The network receives the same cohort- and site-level information that
the corresponding mechanistic process equation would use – diameter at
breast height, number of trees, available light, species identity and
the site's environmental predictors (plus the growth rate, for the
mortality process) – and predicts the process output directly. Species
identity is passed through a learned embedding layer rather than treated
as a categorical covariate with per-species coefficients, so the network
can learn low-dimensional "contrasts" between species or plant
functional types (PFTs). An inverse-link function appropriate to each
process is applied to the network's output: a sigmoid for mortality (a
per-tree death probability, as for
[mortality](https://finnverse.github.io/FINN/reference/mortality.md)),
and an exponential for growth and regeneration (always-positive rates;
for regeneration this is the mean of the negative binomial distribution
from which recruits are drawn, as for
[regeneration](https://finnverse.github.io/FINN/reference/regeneration.md)).

Two network architectures are available, selected with `transformer`.
The default (`transformer = TRUE`) is a small transformer encoder
(`encoder_layers` layers, embedding dimension `emb_dim`, feed-forward
dimension `dim_feedforward`) that embeds each cohort, its species and
the site's environment and attends across the cohorts of a patch.
Setting `transformer = FALSE` instead uses a feed-forward network with
hidden layers `hidden` (default two layers of 50 units each) and
`dropout`, matching the architecture used by Pichler & Käber (2026) for
the Barro Colorado Island case study (where dropout was set to 10%).

As with
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md),
hybrid processes are calibrated jointly, end-to-end, with the remaining
mechanistic or hybrid processes via
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md), rather
than pre-trained in isolation and plugged in afterwards. `optimize`
controls whether the network's weights are estimated during fitting or
kept fixed at their (random) initial values, and `dispersion_parameter`/
`sample_regeneration` have the same meaning as in
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
and are only used when this object is the regeneration process.

## References

Pichler, M., & Käber, Y. (2026). Inferring processes within dynamic
forest models using hybrid modelling. *Methods in Ecology and
Evolution*.
[doi:10.1111/2041-210x.70347](https://doi.org/10.1111/2041-210x.70347)

## See also

[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md),
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md),
[recruitment](https://finnverse.github.io/FINN/reference/recruitment.md)

## Examples

``` r
# A hybrid process is a recipe, like createProcess(): no data and no torch
# are needed until the model is built. The default network is a small
# transformer encoder
growth_process <- createHybrid(formula = ~ temp + prec)
class(growth_process)
#> [1] "hybrid"

# a feed-forward network with two hidden layers instead, as in the Barro
# Colorado Island case study of the paper
mortality_process <- createHybrid(~ temp + prec, transformer = FALSE,
                                  hidden = c(50L, 50L), dropout = 0.1)

library(data.table)
# a model whose growth is a network while mortality, regeneration and
# competition stay mechanistic, fitted for two epochs on three sites of the
# bundled Oregon FIA data (?fia_obs_dt)
sites <- 1:3
Nsp   <- max(fia_obs_dt$species)
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
FINN.seed(1)
m <- finn(N_species = Nsp, recruits_dbh = 12.9,
          growth_process       = createHybrid(~ temp + prec, transformer = FALSE,
                                              hidden = c(16L, 16L), dropout = 0.1),
          mortality_process    = createProcess(~ temp + prec, mortality),
          regeneration_process = createProcess(~ temp + prec, regeneration_capped))
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
    plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.768)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%
# the network takes the place of the growth equation in the simulation
sim <- simulateForest(m, env = fia_env_dt[siteID %in% sites], init_cohort = init)
head(sim$long$site[variable == "growth" & year == 10 & value > 0])
#>    siteID  year species variable      value
#>     <int> <int>   <int>   <fctr>      <num>
#> 1:      1    10       1   growth 0.03613109
#> 2:      1    10       3   growth 0.04689414
#> 3:      1    10       4   growth 0.04193100
#> 4:      1    10       5   growth 0.03449138
#> 5:      1    10       6   growth 0.05649123
#> 6:      1    10      11   growth 0.03411257
```
