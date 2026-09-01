# FINN 0.2.0

This release adds the model features developed for the **FIA** (US Forest
Inventory and Analysis) and **BWI** (German Bundeswaldinventur) calibration
projects. Every new argument defaults to the previous behaviour, so code written
for 0.1.0 runs unchanged.

## New features

* **Per-predictor environmental scaling (`env_autoscale`).** Previously every
  environmental predictor was standardised the same way. `env_autoscale` now takes
  a per-predictor specification — `"auto"` (z-score, the default), `"identity"`
  (leave as supplied), `"0to1"` (min-max), or a user-supplied `function(x)`
  returning `list(center=, scale=)` — given either as a single value applied to
  all predictors or as a named vector/list per predictor. Every mode is affine, so
  ALE interpretation stays invertible. The fitted constants remain accessible via
  `compute_env_scaling()` / `apply_env_scaling()`. Dispatch is fully backward
  compatible: the previous `TRUE` / `FALSE` / `NULL` forms behave exactly as before.

* **Regeneration saturation (`finn(regeneration_saturation = )`).** An optional
  Beverton-Holt cap on recruitment that keeps regeneration bounded at high
  propagule pressure. The fitted carrying capacity `K` (stems ha⁻¹ step⁻¹) is read
  back with the new exported `reg_saturation_K()`. Off by default.

* **Recruitment observation weight (`finn(recruit_obs_weight = )`).** Scales the
  regeneration term in the loss so recruitment observations can be up- or
  down-weighted relative to the other responses during calibration.

* **Period-scale growth scoring (`finn(growth_period_scale = )`).** Scores growth
  over a multi-year remeasurement interval rather than a single step, for
  inventories whose growth is only observed between remeasurements.

* **Custom process parameters (`createProcess(custom_parameters = )`).** Lets a
  user-supplied process function declare and calibrate its own named parameters,
  beyond FINN's built-in per-species and per-environment coefficients.

* **Regeneration observation operator.** An internal hook mapping latent
  recruitment to what an inventory actually records — a building block for the
  inventory simulator.

## Bug fixes

* `createProcess(NN = ...)` is no longer silently dropped inside `create_nn()`.


# FINN 0.1.0

First public release.

FINN is a differentiable forest gap model: a cohort-based dynamic vegetation
model whose demographic processes (competition, growth, mortality, regeneration)
can each be a mechanistic function, a neural network, or a mixture of the two,
all calibrated end-to-end by gradient descent.

Highlights of the current interface:

* `finn()` assembles a model from one process per demographic component, each
  built with `createProcess()` (mechanistic) or `createHybrid()` (neural
  network).
* `simulateForest()` runs a model forward; `fit()` calibrates one to data and
  `predict()` scores it, returning patch- and site-level results.
* Per-response likelihoods: Gaussian (`"mse"`/`"gaussian"`), Poisson, negative
  binomial, and a `"binomial"` likelihood for mortality that takes a
  closed-cohort count pair (`n_at_risk`, `n_died`).
* `weights = "auto"` (the default in `fit()`) scales each loss by its
  intercept-only baseline, so the six responses are commensurable and every term
  reads as a fraction of its own null deviance.
* Data helpers `makeObsData()`, `resolveSiteIDs()` and `makeInitCohorts()` turn a
  raw tree list into FINN's input tables.
* Model interpretation with `ALE()`, `summary()`, `feature_importance()` and
  `conditionalEffects()`.

See the vignettes for a guided tour: *Introduction to FINN*, *Plausible
succession from a handful of species*, *Preparing your data for FINN*, *Fitting
FINN to forest inventory data*, and *Mortality*.
