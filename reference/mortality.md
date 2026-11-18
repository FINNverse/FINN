# Mortality

FINN's mechanistic mortality process: the probability that a tree of a
cohort dies in one simulation step. It is a logistic function of four
terms,
`p = sigmoid(pred + light * parMort[, 1] + parMort[, 2] * dbh/100 + growth/growth_scale * parMort[, 3])`,
with one coefficient per species for each of shade, size and growth.
`pred` is the species' environmental response on the logit scale, so a
species with no environmental predictors dies at a constant rate. Growth
enters as the rate the growth process produced in the same step, which
is what makes mortality depend on a tree's recent performance rather
than on its size alone.

The probability is clamped to `[mortality_floor, mortality_ceiling]`
(default `0.001` / `0.999` per simulation step) and the number of deaths
is then drawn with `mortality_draw`; both are settings of
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md).

## Usage

``` r
mortality(
  dbh,
  species,
  trees,
  parMort,
  pred,
  light,
  base_steepness = 5,
  debug = FALSE,
  growth = NULL,
  growth_scale = 0.01
)
```

## Arguments

- dbh:

  torch.Tensor Diameter of each cohort, in cm.

- species:

  torch.Tensor Species index of each cohort.

- trees:

  torch.Tensor Number of trees in each cohort.

- parMort:

  torch.Tensor Mortality parameters per species, three columns: the
  shade, size and growth coefficients of the logit.

- pred:

  torch.Tensor Environmental response per species, on the logit scale.

- light:

  torch.Tensor Available light for each cohort.

- base_steepness:

  numeric Not used by this function; kept in the signature for custom
  mortality functions written against it.

- debug:

  logical If TRUE, return a list with the light (`light`), the mortality
  probability before any further step (`environment`) and the returned
  probability (`predM`, identical to `environment`). Defaults to FALSE.

- growth:

  torch.Tensor (Optional) Growth entering the mortality response;
  defaults to the model's current growth.

- growth_scale:

  numeric Growth enters the logit as `growth / growth_scale`. The
  default 0.01 expresses the annual relative growth in %/yr (about one
  standard deviation of observed growth), so that the growth coefficient
  `parMort[, 3]` acts per 1 %/yr and its bounds of +-10 allow a
  meaningful effect. With the raw rate (`growth_scale = 1`, the
  behaviour before 0.2.0) the term spans only a few tenths on the logit
  between slow and fast trees even at the bound. Set it for a fitted
  model through
  `createProcess(func = FINN::mortality, growth_scale = 1)`, not through
  a wrapper function: the kernel reads `self$record_raws` from its own
  environment and a wrapper leaves that unbound.

## Value

A `torch` tensor of per-cohort mortality probabilities per simulation
step, `[sites, patches, cohorts]`.

## The starting value

The environmental response is `sigmoid(pred)`. Before any fitting it
starts at about 2 % mortality per time step: the intercept of the
response is `qlogis(0.02)` and its other weights are 0 (with
`createProcess(hidden = )` the bias of the output layer carries the
intercept), and the species terms for light, size and growth move it by
at most about one logit unit. Until 0.2.0 the weights started at torch's
random values, i.e. at `sigmoid(0) = 0.5` per step, which on yearly
steps emptied an unfitted stand within a decade.
`createProcess(initEnv = )` sets another start; a formula without an
intercept (`~ 0 + temp`) leaves torch's start in place.

## See also

[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
(`mortality_draw`, `mortality_floor`, `growth_scale`),
[`growth()`](https://finnverse.github.io/FINN/reference/growth.md),
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md),
[`predictTrees()`](https://finnverse.github.io/FINN/reference/predictTrees.md)
to evaluate the fitted function on observed trees

## Examples

``` r
# mortality as the mortality process of a model, responding to one predictor;
# growth enters in %/yr (growth_scale = 0.01, the default)
m <- finn(N_species = 2,
          mortality_process = createProcess(~ 1 + temp, mortality, growth_scale = 0.01))

# a 30-year simulation from bare ground; an unfitted model uses `env` as
# given (no scaling constants yet), so the predictor is passed standardised
env <- data.frame(siteID = 1L, year = 1:30, temp = 0)
sim <- simulateForest(m, env = env, patches = 5, patch_size = 0.1)
tail(sim$wide$site)
#>    siteID species  year      dbh        ba trees    growth       mort   reg
#>     <int>   <int> <int>    <num>     <num> <num>     <num>      <num> <num>
#> 1:      1       2    25 31.39856 0.1588265   1.2 0.3833489 0.02072886   0.0
#> 2:      1       2    26 34.06740 0.1730944   1.2 0.3382641 0.02045126   0.0
#> 3:      1       2    27 25.05285 0.1888763   1.8 0.2854405 0.02013500   0.6
#> 4:      1       2    28 27.46775 0.2065060   1.8 0.4132742 0.02094093   0.0
#> 5:      1       2    29 27.18676 0.2260164   2.0 0.3639072 0.02065206   0.2
#> 6:      1       2    30 27.18188 0.2473062   2.2 0.3631251 0.02065450   0.2
#>    r_mean_ha
#>        <num>
#> 1: 0.6073149
#> 2: 0.5978763
#> 3: 0.5862898
#> 4: 0.5724127
#> 5: 0.5558798
#> 6: 0.5365475

# the start of the environmental response: qlogis(0.02) and 0 per species
as.matrix(m$nn_mortality$parameters[[1]])
#>          [,1] [,2]
#> [1,] -3.89182    0
#> [2,] -3.89182    0
```
