# FINN 0.2.0

This release adds the model features developed for the **FIA** (US Forest
Inventory and Analysis) and **BWI** (German Bundeswaldinventur) calibration
projects. Every new argument defaults to the previous behaviour, with three
exceptions that do not reproduce 0.1.0 and are marked in the table below:
`mortality(growth_scale = 0.01)`, `createProcess(light_mask = TRUE)` and
`finn(reg_floor = 1e-3)`. A fourth change reproduces 0.1.0 only on a machine
whose `LC_COLLATE` was already `C`: `resolveSiteIDs()` now numbers sites,
patches and species by radix order instead of the collation locale's.

The two tables are the overview; the sections below them carry the detail.

### What is new since 0.1.0

| what | where | default |
|---|---|---|
| Oregon FIA calibration sample as lazy-loaded datasets | `fia_obs_dt`, `fia_env_dt`, `fia_init_trees`, `fia_obs_test`, `fia_env_test`, `fia_init_test`, `fia_species_dt` | — |
| Per-response likelihood, renamed | `fit(loss_family = )` | as `loss` had it; `loss` still accepted, deprecated |
| How a response is reduced before it is scored: `"mean"`, `"sum"`, `"size_classes"`, `"quantiles"`, `"tree"` | `fit(loss_aggregation = )` | `NULL` (per-cell means) |
| Tree table feeding both per-tree routes | `fit(tree_data = )` | `NULL` |
| Score growth and mortality on each simulated cohort that started as an observed tree | `fit(loss_resolution = )` | `"site"` |
| Evaluate the observed-tree route on a fitted model | `predictTrees()` | new export |
| Per census and epoch: `frac_killed`, `mean_kill_year`, `frac_ghost`; and the per-tree predictions | `$tree_deaths`, `$tree_pred` | under `loss_resolution = "tree"` |
| Multiple shooting: rebuild the stand from the observed tree list at every anchor | `fit(reanchor = , anchor_cohorts = , anchor_mode = , anchor_min_dbh = )` | `FALSE`, `NULL`, `"replace"`, `NULL` |
| Reset the stand to the observation every year | `fit(teacher_forcing = )` | `FALSE` |
| Per-site aggregation windows for mixed remeasurement intervals | automatic from `period_length`; override with `model$period_mode` | `"auto"` |
| A non-finite loss term drops that term, or the whole timestep; occurrences counted | `fit(nonfinite = )`, `$nonfinite_counts` | `"drop_term"` |
| Set-up only, no training | `fit(epochs = 0)` | — |
| Per-predictor scaling: `"auto"`, `"identity"`, `"0to1"` or a function | `fit(env_autoscale = )` | `TRUE`, which z-scores as before; the new strings are opt-in |
| Beverton-Holt cap on recruitment | `regeneration_saturation()` | new export |
| Recruitment limited by conspecific adults | `regeneration_adult()` | new export |
| A process function declares and calibrates its own named parameters | `createProcess(custom_parameters = )` | `NULL` |
| **Additive recruitment floor, was hard-coded at 0.2.** A default that does **not** reproduce 0.1.0 — see below | `finn(reg_floor = )` | **`1e-3`** (was 0.2) |
| Weight of the regeneration term in the loss | `finn(recruit_obs_weight = )` | `1.0` |
| Score growth over the whole remeasurement interval | `finn(growth_period_scale = )` | `FALSE` |
| **Binned competition takes its height quantiles over the live cohorts only.** A default that does **not** reproduce 0.1.0 — see below | `createProcess(light_mask = )` | **`TRUE`** (was unmasked) |
| Bernoulli/Poisson death draw instead of the rounded Gamma | `createProcess(mortality_draw = )`, `createHybrid()`, `binomial_from_poisson()` | `"gamma"` |
| Recruit draw without the `+0.001` offset and the `1e-4` clamp | `createProcess(regeneration_sampler = )`, `createHybrid()`, `rnbinom_torch_refined()` | `"gaussian"` |
| A cohort id identifies the initial cohort it came from, so a simulated cohort traces back to its tree | `simulateForest(return_cohorts = TRUE)`, `$pred` | — |
| **Growth enters `mortality()` in %/yr.** The one default that does **not** reproduce 0.1.0 — existing code changes its answer | `mortality(growth_scale = )`, overridable per model with `createProcess(growth_scale = )` | **`0.01`** (was the raw rate); `growth_scale = 1` restores 0.1.0 |

### What was broken

Most of these were silent. The last column says whether a result changes for
someone who was affected.

| defect | results change? |
|---|---|
| `fit(epochs = 0)` ran two epochs (`1:0`) | yes — it trained when asked not to |
| `createProcess(upper = , lower = )` were overwritten with `NULL`, so every fit ran on the process defaults | yes — a fit asking for tighter bounds ran unconstrained |
| `createProcess(NN = )` was dropped inside `create_nn()` | yes — the supplied network was not used |
| `createProcess(growth_scale = )` was unvalidated: a bare `NA` divides by 1 and silently restores the 0.1.0 objective, `0` saturates mortality | yes, and with nothing to signal it; now an error |
| A `torch_load()`ed model's optimizer was reused although it no longer owned the model's tensors | yes — a crash (`optim_ignite_*`) or no training at all (`optim_adam`) |
| Observed-tree, size-class and quantile tensors of an earlier `fit()` were still evaluated in `forward()` | yes — a stale term entered the loss |
| `fit()` set-up looked site and year up by linear scan, O(rows x sites) | no — arrays bit-identical; ~20 min of set-up per fit at 32,686 sites |
| Observed-tree tensors were sized by `patches`, not by the patches the trees occupy | no — light is per patch; `predictTrees(patches = )` is now ignored |
| Per-tree growth and mortality were evaluated at the interval's last year only, against interval observations | yes, for the per-tree routes (new in 0.2.0) |
| `predictTrees()` returned rows in its own order, and failed on a model that had never been run | yes — misaligned predictions |
| A tracked cohort the simulation killed left the scored set mid-window, so the scored sample was a function of the parameters | yes, for `loss_resolution = "tree"` (new in 0.2.0) |
| A retained dead cohort, and any empty padding slot, moved the height quantiles of the binned competition and so the light of the living | yes — masking is now the default; `light_mask = FALSE` keeps the old behaviour |
| `torch_nanquantile` over a patch with no live cohort returned `NaN` into the likelihood | yes — the fit aborted; the quantiles are now 0 |
| A per-tree term was scored under a weight taken from the per-tree null while the cell term ran | yes, on data with a constant census period |
| With sites observed in different years, each response's loss summed one mean per observation year, so a year with ten trees weighed as much as one with thousands | yes, for `period_mode = "per_site"` — the loss is now one pooled mean per batch |
| `period_mode = "constant"` with a per-tree term was silently overridden, with a false reason | yes — now an error naming both halves |
| Cohort ids were enumerated flat over the whole array, losing the tree each cohort came from | no — the link was recoverable, just not without a positional formula |
| `resolveSiteIDs()` numbered sites, patches and species through `LC_COLLATE`, so the same tables gave different integer ids on a German and on a `C` machine | yes, off a `C` locale — the species number indexes `par_growth` / `par_mortality`, so a parameter row belonged to another species |
| `$tree_deaths` counted an earlier window's ghosts among this window's kills, and recorded liveness before the death draw | diagnostic only — `frac_killed` was identically 0 on a one-year window |
| `$tree_pred` was not reset by `fit()`, and `$tree_deaths` reported `NaN` for a census with no fate target | diagnostic only |
| A U+039B and a U+00D7 in `?fit`, a stray `\cr`, two dangling `\link` targets, a `@param` on the wrong function | no — the LaTeX manual and two `R CMD check` WARNINGs |

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
  `loss_resolution = "tree"` reports `$tree_deaths` per census and epoch: the
  share of inventory-alive trees the simulation killed inside the window
  (`frac_killed`), the mean year it killed them (`mean_kill_year`), and the share
  that entered the window already dead from an earlier one (`frac_ghost`).
  Liveness is recorded after the death draw as well as before it: on the
  entering count alone a cohort killed in the window's LAST year is
  indistinguishable from a survivor, so both numbers were undercounts and
  `frac_killed` was identically 0 on a one-year window. On the three-site,
  24-tree fixture of `tests/testthat/test-tree-loss.R` the correction moves
  `frac_killed` from 0 to 0.167 and 0.208 on two one-year windows, from 0.250 and
  0.042 to 0.292 and 0.125 on two three-year windows, and from 0.458 to 0.625 on
  a single four-year window, with `mean_kill_year` rising 1.82 to 2.40 there.
  `frac_ghost` is unchanged.

* **Training regimes: `reanchor` and `teacher_forcing`.** Multiple shooting
  (`reanchor = TRUE` with `anchor_cohorts`, optionally `anchor_mode = "merge"`
  below `anchor_min_dbh`) rebuilds the simulated stand from the observed tree
  list at every census; `teacher_forcing = TRUE` resets it every year. The
  development branch's `fit(shooting = "multiple")`, which re-anchored to the
  tree table itself, was dropped before release in favour of these.

* **Two opt-in numerical repairs to the stochastic draws
  (`createProcess(mortality_draw = , regeneration_sampler = )`).** Both default to
  the previous behaviour, which is bit-identical (`tests/testthat/test-identity.R`,
  config `sample_reg`). `mortality_draw = "poisson"` replaces the Gamma matched to
  the binomial's first two moments by a Bernoulli draw where the cohort holds one
  tree or fewer and a Poisson capped at the cohort size otherwise: the default
  draw is rounded straight-through by the forward pass and then delivers 0.68-0.72
  of the calibrated deaths for one-tree cohorts at annual probabilities of
  0.01-0.10, and its pathwise gradient runs through a quantity the simulation
  discards. `regeneration_sampler = "refined"` drops the `+0.001` Gamma offset and
  the `1e-4` probability clamp (both guard the Gaussian branch against a negative
  draw and a NaN gradient through `log(0)`; neither hazard applies to
  `torch_poisson`) and samples exactly below `lambda = 2`. The offset and the clamp
  are a recruitment floor of 0.0011 per cell-year, which at a patch of 0.0168 ha
  is 0.066 recruits ha^-1 yr^-1 per species -- over a hundred-odd species, as
  large as a temperate inventory's whole ingrowth. Both repairs move the
  simulated stand at equilibrium by a factor of several in basal area and stem
  number while barely moving the fitted loss: the holdout cannot separate them,
  the simulated stand can.

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

* **The recruitment floor is `1e-3`, not 0.2.** `finn(reg_floor = )` replaces a
  hard-coded additive floor of 0.2 on the expected recruitment of every species
  in every cell-year. At 0.2 the floor alone delivered about 13 stems per hectare
  per year over a 65-species FIA pool, against roughly 0.7 observed -- the floor
  was the recruitment. It exists to keep the negative-binomial likelihood away
  from a zero mean, which `1e-3` does. `finn(reg_floor = 0.2)` restores 0.1.0.

* **The binned competition masks padding out of its height quantiles.**
  `createProcess(light_mask = )` now defaults to `TRUE`. The cohort dimension is
  padded to the widest site x patch **of the batch** and an empty slot has height
  ~0, so unmasked the bin edges are set by how many empty slots happen to sit
  beside a stand. Measured on ten live cohorts of 55 to 8 cm: beside a batch
  neighbour holding ten cohorts their light runs 1.00, 0.76, 0.58, 0.44, 0.34,
  0.27, 0.22, 0.19, 0.17, 0.16; beside a neighbour holding a hundred it is 1.00
  for all ten -- competition disappears, with the stand itself unchanged. Batch
  composition is not part of the model, so the unmasked quantile made the fitted
  objective a function of `batchsize`. This is a defect rather than a modelling
  choice, which is why it changes a default instead of arriving as an option.
  `createProcess(light_mask = FALSE)` reproduces results published before 0.2.0.
  `continuous = TRUE` compares pairs and was never affected.

* **Site, patch and species ids are radix-ordered, not locale-ordered.**
  `resolveSiteIDs()` assigned them through `as.factor()` and `sort()`, both of
  which order by `LC_COLLATE`: `c("Acer", "abies_x")` is that order under `C`
  and the reverse under `de_DE.UTF-8`. The species number indexes the fitted
  parameter matrices, so the same input file gave `par_growth[1, ]` to one
  species on one machine and to another species on the next, with nothing in the
  output to show it. Ordering is now `sort(..., method = "radix")`, which ignores
  the locale. Results are unchanged for anyone who ran under `C`; everyone else
  gets a different, and now reproducible, numbering, so ids in artefacts built
  before this release do not necessarily line up with ids built after it.

* **Growth enters `mortality()` in %/yr.** The default mortality process now
  uses `growth / growth_scale` with `growth_scale = 0.01`, so the growth
  coefficient `parMort[, 3]` acts per 1 %/yr of annual relative growth. On the
  raw rate the term spanned only a few tenths on the logit between slow and fast
  trees even at the coefficient bound of -10 (in the BWI calibration 13 of 15
  species sat at that bound). This is the one default that does not reproduce
  0.1.0: `createProcess(func = FINN::mortality, growth_scale = 1)` recovers the
  old behaviour. Saved models keep the process function they were fitted with.

* **Per-site path: one mean per batch, not one per observation year.** With
  sites observed in different years (`period_mode = "per_site"`), each
  response's loss was the MEAN over the units of each observation year, and the
  years' means were SUMMED. A year with ten trees then weighed as much as one
  with thousands (in the BWI split-census set-up the pre-drought half of 10 % of
  the trees carried 75 % of the mortality objective, the drought half 3 %). The
  loss is now one mean over all units scored in the batch (each year's mean
  weighted by its share of trees / `n_at_risk` / `growth_n` / cells), which is
  also what the `weights = "auto"` baselines measure; the reported history is
  that pooled mean. A batch that scores one observation year is bit-identical
  to before; the constant-period path is unchanged.

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

* `createProcess(growth_scale = )` is validated. The value is used as a bare
  divisor inside `mortality()`, so every wrong value was silent: a bare `NA` is
  logical and coerces to `TRUE`, so `growth / NA` divides by 1 and reproduces the
  0.1.0 objective exactly, with no warning and no `NaN`; `NA_integer_` divides by
  `-2^31`; `growth_scale = 0` saturates mortality to 0 or 1 across the whole
  stand. Anything but a single positive finite number, or `NULL`, is now an error
  at `createProcess()` time. Saved models are unaffected.

* `resolveSiteIDs()` orders site, patch and species names with
  `sort(..., method = "radix")` instead of the collation locale; see *Behaviour
  changes*. Covered by `tests/testthat/test-locale-stable-ids.R`, which runs the
  function under `C` and `de_DE.UTF-8` on names the two locales order
  differently.

* `fit()` resets `$tree_pred` along with `$tree_obs` and `$tree_deaths`. It is
  written only by the per-tree branch of `forward()`, so a plain fit after a
  per-tree one left the earlier call's predictions in place.

* `$tree_deaths` skips a census that carries no fate target at all (`died` is
  `NA` for every row of it) instead of reporting `NaN` for `frac_killed`,
  `window` and `frac_ghost`, which are `0/0` there.


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
