# Define a demographic process for FINN

Configures one demographic process (growth, mortality, regeneration or
competition) for use as `mortality_process`, `growth_process`,
`regeneration_process`, or `competition_process` in
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md), either
as a fully mechanistic process with an explicit, interpretable
functional form, or – if `hidden` is supplied – as the first ("Level 1")
level of hybridization described by Pichler & Käber (2026): only the
process' environmental-response function is replaced by a small
feed-forward neural network, while the rest of the process equation
remains mechanistic. To replace the entire process equation with a
neural network ("Level 2"), use
[`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md)
instead.

## Usage

``` r
createProcess(
  formula = NULL,
  func,
  initSpecies = NULL,
  initEnv = NULL,
  hidden = NULL,
  optimizeSpecies = FALSE,
  optimizeEnv = TRUE,
  inputNN = NULL,
  outputNN = NULL,
  dispersion_parameter = 1,
  NN = NULL,
  upper = NULL,
  lower = NULL,
  dropout = 0,
  sample_regeneration = TRUE,
  n_quantiles = 10L,
  continuous = FALSE,
  light_mask = TRUE,
  custom_parameters = NULL,
  mortality_draw = c("binomial", "bernoulli", "poisson", "gamma"),
  regeneration_sampler = c("refined", "gaussian"),
  growth_scale = NULL,
  cap = c("all", "species"),
  prior = "max_obs_recruits",
  prior_lambda = 0.1,
  theta_floor = 0.5,
  optimizeTheta = TRUE,
  mortality_floor = 0.001,
  mortality_ceiling = 0.999,
  reg_floor = 0.001
)
```

## Arguments

- formula:

  (`formula`)  
  Environmental predictors for this process' environmental-response
  function, evaluated against the `env` data. Default `NULL` uses `~.`
  (all available environmental covariates).

- func:

  (`function`)  
  Mechanistic process equation, e.g.
  [growth](https://finnverse.github.io/FINN/reference/growth.md),
  [mortality](https://finnverse.github.io/FINN/reference/mortality.md),
  [regeneration](https://finnverse.github.io/FINN/reference/regeneration.md),
  or
  [competition](https://finnverse.github.io/FINN/reference/competition.md),
  or a custom function with the same arguments. Required.

- initSpecies:

  (`matrix`)  
  Optional custom initial values for the species-specific process
  parameters. Default `NULL` draws random initial values.

- initEnv:

  (`list`)  
  Optional custom initial values for the parameters of the
  environmental-response network/regression. Default `NULL` draws random
  initial values. For the linear/logistic response (no `hidden`) this is
  a matrix with one row per species and one column per term of `formula`
  (the intercept first when the formula has one). Without it the
  mortality response starts at about 2 % per time step (intercept
  `qlogis(0.02)`, other weights 0; with `hidden` layers the output bias
  carries it) and the other processes at torch's random start; e.g.
  `initEnv = matrix(c(-3, 0), N_species, 2, byrow = TRUE)` for
  `~ 1 + temp` starts mortality at 5 %.

- hidden:

  ([`integer()`](https://rdrr.io/r/base/integer.html))  
  Hidden-layer sizes of the feed-forward neural network that replaces
  the environmental-response function (Level 1 hybrid). Default `NULL`
  keeps the environmental response mechanistic (linear/logistic).

- optimizeSpecies:

  (`logical(1)`)  
  Should the species-specific process parameters be estimated during
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)? Default
  `FALSE`.

- optimizeEnv:

  (`logical(1)`)  
  Should the parameters of the environmental-response function be
  estimated during
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)? Default
  `TRUE`.

- inputNN, outputNN:

  Not used in 0.2.0: the input dimension of the environmental-response
  network always follows from `formula` and `env`, and its output
  dimension is the number of species. Leave at `NULL`.

- dispersion_parameter:

  (`numeric(1)`)  
  Initial dispersion parameter of the negative binomial distribution
  used to sample recruits (the size `theta`), or its fixed value with
  `optimizeTheta = FALSE`. Only used when this process is the
  regeneration process. See
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md).

- NN:

  (`nn_module`)  
  Optional custom `torch` module overriding the default
  environmental-response network architecture.

- upper, lower:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html))  
  Boundaries (natural scale) for the species-specific process
  parameters, one value per parameter of the process. They are guard
  rails that keep the optimiser inside the range in which the equations
  stay meaningful, not ranges taken from the data, and
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) does not
  move them. The defaults, as `lower` / `upper`: growth `c(0.01, 0.01)`
  / `c(0.99, 4.0)` (light requirement of the shade response, size
  decay); mortality `c(-10, -10, -10)` / `c(10, 10, 10)` (the shade,
  size and growth coefficients of the logit); regeneration `0.01` /
  `0.99` (light requirement); competition `c(0.3, 0.0)` /
  `c(log(2.3), 2.0)` (height parameter, competition strength), where the
  height parameter is capped so that no species grows taller than 130 m
  (see
  [`height()`](https://finnverse.github.io/FINN/reference/height.md)). A
  parameter that ends a fit at its bound is not identified by the data;
  widening the bounds is a decision for 0.2.1.

- dropout:

  (`numeric(1)`)  
  Dropout rate of the environmental-response network. Ignored unless
  `hidden` is set.

- sample_regeneration:

  (`logical(1)`)  
  Should recruits actually be sampled from the negative binomial
  regeneration distribution (`TRUE`), or only its expectation used
  (`FALSE`)? Only used when this process is the regeneration process.
  With `FALSE` the expected number of recruits per patch is rounded, so
  wherever it is below 0.5 no recruit enters at all: a silent zero for
  sparse species on small patches. See
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md).

- n_quantiles:

  (`integer(1)`)  
  Number of height classes used to discretize cohorts when computing
  shading/light availability. Only used when this process is the
  competition process and `continuous = FALSE`.

- continuous:

  (`logical(1)`)  
  Compute shading continuously for every pair of cohorts (`TRUE`)
  instead of binning cohorts into `n_quantiles` height classes (`FALSE`,
  default). Only used when this process is the competition process.

- light_mask:

  (`logical(1)`)  
  Take the height quantiles of the binned competition over the LIVE
  cohorts only. The cohort dimension is padded to the widest site x
  patch in the batch and an empty slot has height 0.001, so with enough
  padding every live cohort is searchsorted into the top bin, where
  nothing shades anything: measured on a 251 m2/ha stand of ten live
  cohorts, 100 padded slots put all ten at light 1.00. Only the quantile
  INPUT is masked, so every cohort keeps its real height as the bin
  query. `TRUE` is the default since 0.2.0, because the padded width is
  set by the widest site in the BATCH: the same stand of ten cohorts
  gets light 1.00 throughout beside a neighbour holding a hundred, and
  its correct gradient down to 0.159 beside a neighbour holding ten.
  Batch composition is not part of the model, so the unmasked quantile
  made the fitted objective depend on `batchsize`. `FALSE` reproduces
  results published before 0.2.0. Ignored when `continuous = TRUE`,
  which compares pairs and is correct at any padding.

- custom_parameters:

  (`named list` or `NULL`)  
  Extra parameters that a custom `func` needs, declared here so they are
  registered on the model by
  [`finn()`](https://finnverse.github.io/FINN/reference/finn.md),
  optimised jointly by
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) (when
  trainable), and round-tripped by `torch_save()`/`torch_load()`. Each
  element is either a numeric vector (a trainable parameter with that
  starting value) or a list with `init` (numeric, length 1 or
  `N_species`) and optional `trainable` (`logical`, default `TRUE`;
  `FALSE` makes it a fixed buffer). The custom function accesses each as
  `self[[name]]`. Parameterise in the *unconstrained* space and apply
  any link inside the function (e.g. store `log K` and use
  `exp(self$logK)`), mirroring FINN's own `par_*_unconstrained`
  convention.

- mortality_draw:

  (`character(1)`)  
  How the number of dead trees is drawn from the cohort's mortality
  probability. Only used when this process is the mortality process.
  `"binomial"` (default since 0.2.0) is `torch`'s exact binomial
  sampler: integer-valued, exact in both moments, at constant memory. It
  has no gradient of its own, which costs nothing under
  `fit(rb_gradient = TRUE)` (the default), where the gradient is taken
  through the expected number of deaths. Under
  `fit(rb_gradient = FALSE)` there is no such path, and `"binomial"`
  falls back to the `"gamma"` draw. `"bernoulli"` gives each tree of the
  cohort its own reparameterised Bernoulli trial
  ([`binomial_from_bernoulli()`](https://finnverse.github.io/FINN/reference/binomial_from_bernoulli.md)),
  which is integer-valued, matches the binomial in both moments –
  measured 12.011 / 8.387 against an exact 12.000 / 8.400 at `n = 40`,
  `p = 0.3` – and keeps a gradient path to `p` through a
  straight-through sigmoid. It expands each cell to the largest cohort
  in the batch, so its memory grows with the largest `trees` value in
  the data. `"poisson"` draws a Bernoulli where the cohort holds one
  tree or fewer and a Poisson capped at the cohort size otherwise:
  integer-valued, exact in expectation for small `n*p`, and no pathwise
  gradient – which costs nothing now that `fit(rb_gradient = TRUE)`
  takes the gradient through the expectation anyway. The Bernoulli
  branch matters because capping a Poisson at a one-tree cohort delivers
  only 88 % of the calibrated probability at the `p = 0.26` of a
  suppressed stem. `"gamma"` (the behaviour of every published FINN
  result) matches a Gamma to the binomial's first two moments and
  `rsample()`s it, which gives a pathwise gradient – the hazard a
  discrete draw does not cover – but the forward pass rounds that draw
  straight-through, so the gradient flows through a quantity the
  simulation discards and the rounded count delivers a measured
  0.68-0.72 of the calibrated deaths for one-tree cohorts at annual
  probabilities of 0.01-0.10. See
  [`binomial_from_poisson()`](https://finnverse.github.io/FINN/reference/binomial_from_poisson.md)
  for the measured effect on the simulated stand.

- regeneration_sampler:

  (`character(1)`)  
  How recruits are drawn from the negative binomial. Only used when this
  process is the regeneration process and `sample_regeneration = TRUE`.
  `"refined"` (default since 0.2.0) samples exactly below `lambda = 2`
  and drops both guards below, neither of which applies to
  `torch_poisson`; see
  [`rnbinom_torch_refined()`](https://finnverse.github.io/FINN/reference/rnbinom_torch_refined.md).
  Above `lambda = 2` it is still the floored Gaussian draw, which
  delivers about half a recruit per draw too few (measured over 1e6
  draws: -0.43 at `lambda = 2.5`, -0.50 from `lambda = 20`, against
  about 0 with rounding). Rounding instead of flooring is planned for
  0.2.1; it changes results, so it is not in this release. `"gaussian"`
  (the behaviour of every published FINN result) offsets the Gamma draw
  by `+0.001`, clamps the probability at `1e-4` and pushes the result
  through `lambda + eps*sqrt(lambda)`. The offset and the clamp guard
  the Gaussian branch against a negative draw and against a NaN gradient
  through `log(0)`; together they are a recruitment FLOOR of 0.0011 per
  cell-year, which at a FINN patch of 0.0168 ha is 0.066 recruits ha^-1
  yr^-1 per species – over a pool of a hundred-odd species, as large as
  a temperate inventory's whole ingrowth.

- growth_scale:

  (`numeric(1)` or `NULL`)  
  Overrides the `growth_scale` default of the process `func`. Only used
  when this process is the mortality process and `func` declares that
  argument. `NULL` (the default) keeps the function's own default, which
  for
  [`mortality()`](https://finnverse.github.io/FINN/reference/mortality.md)
  is `0.01` – growth enters the logit in %/yr. `growth_scale = 1` feeds
  it the raw annual relative rate, the behaviour of 0.1.0. The value is
  set as the function's formal rather than passed down the forward pass,
  so a custom mortality `func` without the argument keeps working
  unchanged.

- cap, prior, prior_lambda:

  The Beverton-Holt cap of
  [`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md):
  `cap = "all"` (default, one cap shared by all species) or `"species"`;
  `prior = "max_obs_recruits"` (default, the largest observed
  recruitment rate in the fitted data, recruits per ha and year) or a
  positive number; `prior_lambda` (default `0.1`) the strength of the
  pull of `log K` towards `log prior`. See
  [`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md).
  Only used when this process is the regeneration process and its `func`
  uses the cap.

- theta_floor:

  (`numeric(1)`)  
  Lower bound on the negative-binomial recruit dispersion
  `par_theta_recruits`, applied wherever it is read: the recruit draw
  and the regeneration likelihood. Default `0.5` since 0.2.0; `0` or
  `NULL` leaves it unbounded, as in 0.1.0. Only used when this object is
  the regeneration process.

  theta is trained, and it collapses, because lowering it widens the
  regeneration likelihood and so cheapens that term without improving
  the predicted mean. The same parameter sets the variance of the
  *simulated* recruit draw, `mu + mu^2/theta`, which then turns bursty:
  on the German NFI it reached 0.035 by epoch 50 and the simulated stand
  overshot observed recruitment 133-fold, with basal-area and stem
  losses flat until epoch 25 and running away after. The Beverton-Holt
  cap does not prevent this – it bounds the mean, and the variance
  around a capped mean is still `K + K^2/theta`.

  The floor keeps the fit stable, and in one converged pair of fits on
  the German national forest inventory it calibrated the recruitment
  level as well as a fixed `theta`: out of sample on a spatially blocked
  holdout, 1.06 times the observed level after 2000 epochs (about 126
  000 optimizer steps) against 1.04 for `optimizeTheta = FALSE`. A fit
  stopped early is different: after 50 of those epochs (3 150 optimizer
  steps) the floored `theta` still put the level about 12 times too
  high, a fixed `theta` 2 to 3 times. Read the level only once the
  parameters have settled, or fix `theta` for short fits. See
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md).

- optimizeTheta:

  (`logical(1)`)  
  Should the recruit dispersion `par_theta_recruits` be estimated during
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)? Default
  `TRUE`. With `FALSE` it stays at `dispersion_parameter` for the whole
  fit (a buffer, not a parameter), in the recruit draw and in the
  regeneration likelihood alike; `theta_floor` still applies. Only used
  when this process is the regeneration process.

  Use `FALSE` for short fits, or when the recruitment level should not
  depend on how long the fit runs: a trained `theta` falls on sparse
  recruitment data, and with the floor (`theta_floor`) the level is
  right only once the fit has converged, while a fixed `theta` gets
  there sooner. `dispersion_parameter = 1` is FINN's start value; on the
  German national forest inventory it matched the floored `theta` at
  convergence. A higher value such as `50` brings the recruit draw close
  to a Poisson draw: the variance `mu + mu^2/theta` approaches `mu`
  while the expected count per patch `mu` stays small against `theta`
  (at `theta = 50` the excess is 2 % of the mean at `mu = 1` and 20 % at
  `mu = 10`). Per-patch recruit counts are usually small enough for
  that; larger counts need a larger value. See
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md).

- mortality_floor, mortality_ceiling:

  (`numeric(1)`)  
  The interval the per-cohort mortality probability is clamped to **per
  simulation step**, before the death draw. Only used when this process
  is the mortality process; the per-tree mortality loss is not bounded
  by either.

  The defaults are `0.001` and `0.999` (0.1 % and 99.9 % per step).
  Until 0.2.0 the interval was hard-coded to `[0.01, 0.99]`, which is
  why a simulated low-mortality stand – a mature unthinned stand, a
  slow-growing species, an equilibrium run – lost stems the observations
  did not: at a 1 % floor it loses 1 % of every cohort per step however
  low the fitted rate.

  The floor exists because
  [`binomial_from_gamma()`](https://finnverse.github.io/FINN/reference/binomial_from_gamma.md)'s
  shape `n*p/(1-p)` degenerates as `p -> 0`. Measured with that draw, it
  stays finite with an unbiased mean gradient down to `1e-4`, but
  single-draw gradients grow heavy-tailed as `n*p` shrinks – the largest
  per-tree gradient is about ten times the old default's at `1e-3` and
  about a hundred times at `1e-4` – so `1e-3` is the lowest tested value
  and is the new default.

- reg_floor:

  (`numeric`)  
  A small recruitment that every species gets in every patch and time
  step, whatever the light and the environment, in trees per ha and time
  step. Only used when this process is the regeneration process
  ([`regeneration()`](https://finnverse.github.io/FINN/reference/regeneration.md),
  [`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md),
  [`regeneration_adult()`](https://finnverse.github.io/FINN/reference/regeneration_adult.md)
  and
  [`regeneration_adult_effect()`](https://finnverse.github.io/FINN/reference/regeneration_adult_effect.md)
  add it to their expected recruitment:
  `mean = regP * exp(env) + reg_floor`). Default `1e-3`, which keeps the
  recruitment likelihood finite where the light response closes and adds
  next to nothing. One value for all species or one per species. Before
  0.2.0 the floor was hard-coded to `0.2`, which with many species adds
  a large recruitment that does not depend on anything (65 species x 0.2
  = 13 trees per ha and year); models saved before then keep that value
  when reloaded. `finn(reg_floor = )`, the argument of earlier versions,
  is still accepted with a warning and overrides this value.

## Value

A list of class `"process"` containing the process definition and
associated parameters, to be passed as `mortality_process`,
`growth_process`, `regeneration_process`, or `competition_process` to
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md).

## Details

Each demographic process in FINN is the product of (i) a process
equation `func` that operates on the cohort state (dbh, number of trees,
available light, species) and species-specific process parameters, and
(ii) a species- and process-specific environmental-response function
that maps site-level environmental predictors to a scalar effect on the
process (see
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md) for the
underlying equations). `createProcess()` configures both parts:

- `func` implements the mechanistic process equation itself. The
  package's default process functions
  ([growth](https://finnverse.github.io/FINN/reference/growth.md),
  [mortality](https://finnverse.github.io/FINN/reference/mortality.md),
  [regeneration](https://finnverse.github.io/FINN/reference/regeneration.md),
  [competition](https://finnverse.github.io/FINN/reference/competition.md))
  reproduce the equations described by Pichler & Käber (2026); a custom
  function with the same arguments can be passed instead to use a
  different functional form while keeping the process embedded in, and
  jointly calibrated with, the rest of the model.

- The environmental-response function is, by default (`hidden = NULL`),
  a linear/logistic niche function with one coefficient per
  environmental covariate (named in `formula`) and per species,
  comparable to a classic species distribution model. Setting `hidden`
  to a vector of hidden-layer sizes (e.g. `c(25L)`) instead replaces
  this function with a feed-forward neural network, while `func` and the
  species-specific process parameters remain mechanistic – the "Level 1"
  hybridization described in the paper.

`formula` selects which columns of the `env` data (passed to
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) or
[`predict.finn_class()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md))
enter the environmental-response function; `initSpecies`/`initEnv` allow
supplying custom starting values for the process
parameters/environmental-response model instead of the package's random
initialization; `optimizeSpecies`/`optimizeEnv` control whether these
parameters are estimated during
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) or held
fixed at their initial values; and `upper`/`lower` set box constraints
(on the natural process-parameter scale) within which the
species-specific process parameters are constrained during optimization.

## References

Pichler, M., & Käber, Y. (2026). Inferring processes within dynamic
forest models using hybrid modelling. *Methods in Ecology and
Evolution*.
[doi:10.1111/2041-210x.70347](https://doi.org/10.1111/2041-210x.70347)

## See also

[`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md),
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md),
[recruitment](https://finnverse.github.io/FINN/reference/recruitment.md)
for the recruitment options and how to combine them

## Examples

``` r
# A process is a recipe: createProcess() needs neither data nor torch, and the
# predictors named in `formula` are looked up in `env` when the model is
# fitted or simulated. Here: growth responding to temperature and precipitation
growth_process <- createProcess(formula = ~ temp + prec, func = growth)
class(growth_process)
#> [1] "process"
deparse(growth_process$formula)
#> [1] "~temp + prec"

# the same equation with its species parameters estimated as well (by default
# only the environmental response is), and a quadratic temperature response
growth_process <- createProcess(~ temp + I(temp^2) + prec, func = growth,
                                optimizeSpecies = TRUE)

# level-1 hybrid: the environmental response of mortality as a small neural
# network with two hidden layers; the mortality equation stays mechanistic
mortality_process <- createProcess(~ temp + prec, func = mortality,
                                   hidden = c(16L, 16L), dropout = 0.1)
#> Neural Network detected...

# the regeneration process with a smaller recruitment floor, a cap per
# species and a fixed recruit dispersion
regeneration_process <- createProcess(~ temp, func = regeneration_capped,
                                      reg_floor = 1e-4, cap = "species",
                                      optimizeTheta = FALSE)

# competition takes no environment; continuous = TRUE compares every pair of
# cohorts instead of binning them into height classes
competition_process <- createProcess(func = competition, continuous = TRUE)

# the processes become a model in finn(), here for the 11 species of the
# bundled Oregon FIA data (?fia_obs_dt), whose env table has `temp` and `prec`
m <- finn(N_species = 11, recruits_dbh = 12.9,
          growth_process       = growth_process,
          mortality_process    = mortality_process,
          regeneration_process = regeneration_process,
          competition_process  = competition_process)
# the species parameters of growth: light requirement and size decay
head(as.matrix(m$par_growth), 3)
#>           [,1]       [,2]
#> [1,] 0.4952571 0.05029631
#> [2,] 0.4675127 0.04575275
#> [3,] 0.5246698 0.04777756
```
