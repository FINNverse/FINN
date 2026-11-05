# FINN 0.2.0

This release adds the model features developed for the **FIA** (US Forest
Inventory and Analysis) and **BWI** (German Bundeswaldinventur) calibration
projects. Every new argument defaults to the previous behaviour, so code written
for 0.1.0 runs unchanged.

## New features

* **Bundled Oregon FIA data as package datasets.** The calibration sample that
  the `D-Fit_to_FIA` and `E-Mortality` vignettes fit on now ships as lazy-loaded
  datasets in FINN's input format: `fia_obs_dt`, `fia_env_dt`, `fia_init_trees`
  (200 training sites), `fia_obs_test`, `fia_env_test`, `fia_init_test` (200
  disjoint holdout sites) and `fia_species_dt` (the shared species coding). They
  are ready for `fit()` / `makeInitCohorts()` without any file reading; see
  `?fia_obs_dt`. The `fia_*.csv` files previously under `inst/extdata` are gone
  (the xz-compressed datasets are ~120 KB versus ~2.2 MB of CSV); the
  `example_*.csv` raw inputs for the data-preparation vignette stay.

* **Scoring individual trees, two ways.** `fit()` can score `growth` and
  `mortality` per observed tree instead of per site, species and census. Both
  options take the same tree table in `tree_data` (one row per tree and census:
  `siteID`, `year`, `patchID`, `species`, `growth`, `died`, ...) and are chosen
  explicitly: `loss_aggregation = c(growth = "tree", mortality = "tree")`
  predicts each OBSERVED tree (its measured `dbh`, light from the observed
  stand; the responses can be chosen separately and the initial cohorts may be
  binned), while `loss_resolution = "tree"` predicts each SIMULATED cohort that
  started as an observed tree along the model's own trajectory (needs a
  `cohortID` column and one initial cohort per tree; combines with `reanchor`
  and `teacher_forcing`). `loss_resolution = "tree"` no longer takes a list of
  per-timestep target arrays — fit() builds those from the table and refuses
  any row it cannot place on an initial cohort. See `?fit`, section *Scoring
  individual trees*, for the comparison. The two cannot be combined on one
  response. `predictTrees()` evaluates the first option on a fitted model.

* **Training regimes: `reanchor` and `teacher_forcing`.** Multiple shooting
  (`reanchor = TRUE` with `anchor_cohorts`, optionally `anchor_mode = "merge"`
  below `anchor_min_dbh`) rebuilds the simulated stand from the observed tree
  list at every census; `teacher_forcing = TRUE` resets it every year. The
  development branch's `fit(shooting = "multiple")`, which re-anchored to the
  tree table itself, was dropped before release in favour of these.

* **`fit(nonfinite = )`.** When one of the six loss terms is not finite at a
  timestep, `"drop_term"` (the default) leaves only that term out, still
  backpropagates the other responses there and counts the occurrences in
  `model$nonfinite_counts`; `"skip_step"` discards the whole timestep, which is
  what earlier versions did.

* **Per-predictor environmental scaling (`env_autoscale`).** Previously every
  environmental predictor was standardised the same way. `env_autoscale` now takes
  a per-predictor specification — `"auto"` (z-score, the default), `"identity"`
  (leave as supplied), `"0to1"` (min-max), or a user-supplied `function(x)`
  returning `list(center=, scale=)` — given either as a single value applied to
  all predictors or as a named vector/list per predictor. Every mode is affine, so
  ALE interpretation stays invertible. The fitted constants remain accessible via
  `compute_env_scaling()` / `apply_env_scaling()`. Dispatch is fully backward
  compatible: the previous `TRUE` / `FALSE` / `NULL` forms behave exactly as before.

* **Saturating regeneration process (`FINN::regeneration_saturation`).** A new,
  exported regeneration process function that adds a Beverton-Holt density cap
  `recruits = K*m/(K+m)` on top of the standard mechanistic recruitment, keeping
  regeneration bounded at high propagule pressure. It is a drop-in alternative to
  `FINN::regeneration`, used as the `regeneration_process` in `finn()`. The fitted
  carrying capacity is a normal process parameter, declared through
  `createProcess(custom_parameters = list(reg_logK = ...))`: a length-1 init gives
  one shared `K`, a length-`N_species` init a per-species `K` (`K = exp(reg_logK)`,
  read back with `exp(as.numeric(m$reg_logK))`). This replaces the earlier
  (unreleased) `finn(regeneration_saturation = )` model-level option, so the cap
  now composes like any other process rather than as a special case.

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

## Behaviour changes

* **Growth enters `mortality()` in %/yr.** The default mortality process now
  uses `growth / growth_scale` with `growth_scale = 0.01`, so the growth
  coefficient `parMort[, 3]` acts per 1 %/yr of annual relative growth. On the
  raw rate the term spanned only a few tenths on the logit between slow and fast
  trees even at the coefficient bound of -10 (in the BWI calibration 13 of 15
  species sat at that bound). This is the one default that does not reproduce
  0.1.0: pass `growth_scale = 1` through a custom process function to recover
  the old behaviour. Saved models keep the process function they were fitted
  with.

## Performance

* **`fit()` set-up no longer scales with sites x rows.** `climateDF2array()`
  (behind every response and environment array `fit()` builds, ~140 calls per
  fit with 15 species) looked each row's site and year up by a linear scan, so
  a call cost O(rows x sites): ~90 s at 32,686 BWI sites, ~20 min of set-up per
  fit. The lookup is now a hash map, O(rows); the arrays are bit-identical.

* **Observed-tree tensors have the trees' own patches.** For
  `loss_aggregation = "tree"` and `predictTrees()`, the tree tensors were sized
  by `patches` (e.g. 25) although the trees occupy fewer patches (all BWI trees
  have `patchID` 1), and competition, growth and mortality ran on the empty
  ones in every batch. They are now sized by the number of distinct `patchID`s
  in `tree_data`. Light is computed per patch, so no prediction changes;
  `predictTrees(patches = )` is now ignored.

## Bug fixes

* `fit(epochs = 0)` ran two epochs (`1:0`); it now runs the set-up only
  (data tensors, `loss_baseline` / `loss_weights`) and trains nothing.

* `fit()` resets the observed-tree, size-class and quantile tensors of an
  earlier call (or a `torch_load()`ed model) before building those of the
  current one; stale tensors were otherwise still evaluated in `forward()`.

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
