# FINN 0.2.0

Integrates the calibration work from the **FIA** (US Forest Inventory and
Analysis) and **BWI** (German Bundeswaldinventur) projects. The FIA
data-preparation pipeline shipped in 0.1.0; this release adds the model features
those projects needed. New arguments default to the previous behaviour, so
existing 0.1.0 code runs unchanged.

## New features

* `finn()` gains `recruit_obs_weight`, `growth_period_scale` and
  `regeneration_saturation` arguments. `regeneration_saturation` adds an optional
  Beverton-Holt cap on regeneration; the new exported `reg_saturation_K()` reads
  the fitted carrying capacity `K` (stems/ha/step) back out of a fitted model.
* New `growth_period_scale` option for period-scale growth scoring.
* `env_autoscale` now supports per-predictor scaling — `"auto"`, `"identity"`,
  `"0to1"`, or user-supplied functions.
* New regeneration observation operator (a hook for the inventory simulator).
* `createProcess()` gains a `custom_parameters` interface.

## Bug fixes

* `createProcess(NN = ...)` is no longer dropped in `create_nn()`.


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
