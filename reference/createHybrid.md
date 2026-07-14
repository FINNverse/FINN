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
  dim_feedforward = 256L
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
  used to sample recruits. Only used when this object is the
  regeneration process.

- NN:

  (`nn_module`)  
  Optional custom `torch` module overriding the default network
  architecture entirely.

- dropout:

  (`numeric(1)`)  
  Dropout rate of the network.

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
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md)

## Examples

``` r
if (FALSE) { # \dontrun{
growth_process <- createHybrid(formula = ~temperature + precipitation)
} # }
```
