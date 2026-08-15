# FINN 0.2.0.9000 (development)

Adds a **prescriptive management layer** and a worked **PROFOUND** example.

## New features

* New **parametric management operator**: every strategy — thinning, target-diameter
  harvest, species conversion, and each WET2024 forest development type — is one point
  in a single continuous parameter set (`management_params()`, the "management
  genome"), never a discrete case. `management_removal()` evaluates the per-cohort
  removal-fraction field (Zielstärkennutzung + Durchforstung); `apply_management()`
  applies it and reports removed stems, basal area, and the mean diameter of the
  removed trees (`Dg_aus`) for validation against yield tables. `thinning_size_bias`
  (Durchforstungsart) may be size-dependent via `graduated_bias()`, reproducing the
  NW-FVA graduated thinning (from below when young, from above near the target).
  Presets `thin_from_below()`, `thin_from_above()`, `target_diameter_harvest()`,
  `species_removal()` and `clearfell()` return parameter points over this operator.
* `simulate_managed()` (in progress) threads a management schedule through a fitted
  model by simulating in segments and re-initialising from the thinned state — no
  FINN core change and no refit.
* New vignette *Fitting FINN to PROFOUND data and simulating management*, calibrating
  FINN to one PROFOUND site and running the regimes above on the fitted model.
* New growth process `growth_env()`: an environment-controlled size-decline. Where
  `growth()` slows diameter growth with a fixed per-species constant
  (`k = parGrowth[,2]`), `growth_env()` uses `k_eff = k_base * exp(-k_env * pred)`
  — a species baseline combined with the environment via the growth predictor, so
  productive sites both grow faster and sustain growth longer. Supplied opt-in via
  `createProcess(custom_parameters = list(k_env = ...))`; `k_env = 0` reduces
  exactly to `growth()`. Motivated by fitting the NW-FVA yield tables, where a
  single species `k` over-/under-shoots the best/poorest site classes.
* NW-FVA growth-and-yield reference tables now ship with FINN as an open CSV,
  loaded via `nwfva_yield_tables()` (5 species, yield classes -1..4; source Albert
  et al. 2021, CC-BY-4.0; see `inst/extdata/nwfva_gy_tables_SOURCE.md`).


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
