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
  continuous = FALSE
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
  initial values.

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

- inputNN:

  (`integer(1)`)  
  Input dimension of the environmental-response network. Default `NULL`
  infers it from `formula`/`env`.

- outputNN:

  (`integer(1)`)  
  Output dimension of the environmental-response network. Default `NULL`
  uses the number of species.

- dispersion_parameter:

  (`numeric(1)`)  
  Initial dispersion parameter of the negative binomial distribution
  used to sample recruits. Only used when this process is the
  regeneration process.

- NN:

  (`nn_module`)  
  Optional custom `torch` module overriding the default
  environmental-response network architecture.

- upper:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html))  
  Upper boundaries (natural scale) for the species-specific process
  parameters.

- lower:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html))  
  Lower boundaries (natural scale) for the species-specific process
  parameters.

- dropout:

  (`numeric(1)`)  
  Dropout rate of the environmental-response network. Ignored unless
  `hidden` is set.

- sample_regeneration:

  (`logical(1)`)  
  Should recruits actually be sampled from the negative binomial
  regeneration distribution (`TRUE`), or only its expectation used
  (`FALSE`)? Only used when this process is the regeneration process.

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
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md)

## Examples

``` r
growth_process <- createProcess(formula = ~temperature + precipitation, func = growth)
```
