# FINN 0.2.1

A bug-fix release. No argument, default or default result changes; each fix
below changes numbers only on the path it names.

## Bug fixes

* `fit()` collects garbage once per epoch. Torch tensors are freed only when R's
  garbage collector runs, which R triggers on its own heap, not on torch's
  allocations, so the autograd graphs of several epochs piled up. Long
  full-batch fits with hybrid processes (`createHybrid()`) on large data sets
  ran out of memory (FIA CONUS, hybrid growth: 37 to 85 GB in five epochs; with
  the fix flat at 38-39 GB). Results are unchanged (#79).
* `makeObsData()` with `Nspecies` or `NspeciesQuantile` lumped the species
  outside the selection into `"other"` for every column except the mortality
  counts, which are grouped by the species of the previous census. Those deaths
  and trees at risk fell out of `"other"` (on the example tree list in
  `inst/extdata` with `Nspecies = 2`: 6 trees at risk and 0 deaths instead of
  162 and 16). Without a species cap nothing changes (#73).
* `fit()` with `loss_aggregation` `"quantiles"` (dbh) or `"size_classes"`
  (trees): a missing value in one of the observed quantile or class columns
  dropped its row, so every later value of that species was scored against
  another site or year. Missing values now stay in their own cell and are not
  scored. Data without missing values in these columns are unaffected (#75).
* `fit()` on the per-site aggregation path (mixed `period_length`) skipped the
  whole step of a census at which no site closes an observation window, so the
  dbh, ba and tree terms of that census were not scored. Data in which every
  census row carries a `period_length` are unaffected (#74).
* A second `fit()` call on the same model, after a call with a learning-rate
  schedule, reused the optimizer at the rate the schedule had left instead of
  the `lr` it was given. Every call now starts at its own `lr`. A model loaded
  with `torch_load()` was not affected (#76).
* `?makeInitCohorts` says how the patch dimension is set: by the `patchID`s
  present in `init_trees`, not by `fit(patches = )`; a patch that starts empty
  needs a row with `trees = 0`. Help pages that named 0.2.1 for the recruit
  draw rounding and the parameter bounds now refer to a later feature release.

# FINN 0.2.0

## Changes that alter results

These change the numbers a fit or a simulation produces. Everything else in this
release is new behaviour behind a new argument, a rename that keeps the old name,
or documentation. Each row names the setting that restores the 0.1.0 result.

| change | new default | restore 0.1.0 |
|---|---|---|
| Offset of the `trees` likelihood at a simulated zero count | `fit(trees_offset = 0.01)` | `trees_offset = NULL` |
| Count gradients taken through the expected deaths and recruits, not through the draw | `fit(rb_gradient = TRUE)` | `rb_gradient = FALSE` |
| Exact binomial death draw instead of the moment-matched Gamma | `createProcess(mortality_draw = "binomial")` | `mortality_draw = "gamma"` |
| Exact recruit draw below `lambda = 2` instead of the Gaussian approximation | `createProcess(regeneration_sampler = "refined")` | `regeneration_sampler = "gaussian"` |
| Limits of the mortality probability, per simulation step | `createProcess(mortality_floor = 0.001, mortality_ceiling = 0.999)` | `0.01` / `0.99` |
| Recruitment capped by a Beverton-Holt term, with a prior on the cap | `finn()`'s default regeneration process is `regeneration_capped(cap = "all", prior = "max_obs_recruits", prior_lambda = 0.1)` | `FINN::regeneration` |
| Lower bound on the recruit dispersion, in the draw and the likelihood alike | `createProcess(theta_floor = 0.5)` | `theta_floor = 0` or `NULL` |
| Growth enters `mortality()` in %/yr instead of as the raw rate | `createProcess(growth_scale = 0.01)` | `growth_scale = 1` |
| Binned competition takes its height quantiles over the live cohorts only | `createProcess(light_mask = TRUE)` | `light_mask = FALSE` |
| Additive recruitment floor | `createProcess(reg_floor = 1e-3)` | `reg_floor = 0.2` |
| A numeric `loss_weights` vector multiplies the baseline-scaled objective instead of replacing it | `fit(weight_scale = "auto")` | `weight_scale = FALSE` |
| `predict()` and `simulateForest()` leave `requires_grad` as they found it instead of switching every parameter on (#55) | always | - |
| `resolveSiteIDs()` numbers sites, patches and species by radix order instead of the collation locale's | always | reproduces 0.1.0 where `LC_COLLATE` was already `C` |
| A hybrid growth process reads the growth predictors and a hybrid mortality process the mortality predictors (#72, B2) | always | - (0.1.0 and 0.2.0-rc fed each the other's; identical where the two formulas are the same) |
| `predict()` / `simulateForest()` from bare ground without `patches` use the patch count of the fit instead of 100 | `patches` not given: `model$patches` | `patches = 100` |

The sections *The stem-count gradient*, *The death draw*, *The recruit draw*, *The
mortality limits* and *The recruitment cap* under **The release in detail** carry
the measurements behind the first six rows.

## The interface of `finn()` and `fit()`

Agreed for the release so that the two central functions read for an
ecologist who does not know modelling terms. Every renamed or moved argument
is still accepted under its old name, with a warning that names the new form,
and gives the same result.

* `finn()` writes its default processes in the signature instead of `NULL`:
  `createProcess(func = mortality)`, `createProcess(func = growth)`,
  `createProcess(func = regeneration_capped)` and
  `createProcess(func = competition)`. They are exactly what `NULL` used to
  build, and `NULL` is still accepted. `?finn` carries the model overview
  figure, says which equation each default is, that competition cannot depend
  on the environment, and has a complete build-estimate-simulate example.

* `finn(growth_rate = c("annual", "period"))` replaces
  `finn(growth_period_scale = )`: it says what the `growth` column of the data
  is -- the mean annual relative growth rate (default; `TRUE` was `"period"`,
  the relative increment over the whole interval).

* The recruitment floor moves to the regeneration process:
  `createProcess(reg_floor = 1e-3)`. `finn(reg_floor = )` is accepted with a
  warning and overrides the process' value.

* `finn(recruit_obs_weight = )` is the last argument; it matters only for
  sampling designs that do not count recruits per hectare, and its
  explanation lives in `?recruitment`.

* `fit()` is titled *Estimate the parameters of a FINN model* and says that
  all processes are calibrated jointly. Its loss arguments stand together,
  each with one line and a link to the new help page `?finn_loss`
  (likelihoods, aggregation, resolution, weights, the per-tree routes,
  `loss_recruit`); the input tables have the new help page `?finn_data`
  (mandatory columns, missing values, incomplete series, what FINN fills in,
  `env`, `disturbance`, `init_cohort`, `tree_data`, the helpers).

* `fit(loss_weights = )` replaces `fit(weights = )`, so the name says these
  are the weights of the loss terms. `fit(weight_scale = )` also takes `TRUE`
  (= `"auto"`, the default) and `FALSE` (= `"raw"`).

* `fit(loss_recruit = )` replaces `fit(recruits = )`; `predict()` and
  `simulateForest()` keep `recruits =`, where it is not a loss setting.

* `fit(control = finnControl())` collects the optimizer settings:
  `optimizer`, `lr_scheduler`, `lr_scheduler_params`, `clip_norm` and extra
  optimizer arguments such as `weight_decay`. `lr`, `epochs` and `batchsize`
  stay arguments of `fit()`. The four are still accepted as arguments of
  `fit()` with a warning.

* `fit(loss = )` is deprecated in favour of `loss_family`; it is kept as the
  last formal argument so that it does not partially match the `loss_*`
  arguments. `loss_family` also takes a partial named vector
  (`c(trees = "nbinom")` keeps the other five defaults).

* `fit()` checks its input before anything runs and says what is wrong: a
  missing key or response column of `data` (a response with
  `loss_family = "none"` may be absent), species codes outside
  `1..N_species`, duplicated rows, a plain `data.frame` is accepted; `env`
  without a predictor column (#64), with a missing value, with sites that
  do not share the same years, or without a year that `data` observes;
  `disturbance` without `intensity`, with values outside `[0, 1]`, or over
  other sites and years than `env`. `fit()` and `predict()` warn when
  `patches` is given with `init_cohort` and differs (it is ignored, #66),
  and when `patch_size` differs from the one the model was last fitted or
  simulated with; `predict()`/`simulateForest()` default `patch_size` to the
  fitted one. `fit(start_time = )` warns that it has no effect in 0.2.0.

* The progress plot of `fit(plot_progress = TRUE)` shows every loss term
  over the whole range: it is redrawn each epoch with the y-axis fitted to
  the data, so a term that rises above its starting value is no longer cut
  off, and a term that is switched off is left out.

* `predict()`/`simulateForest()` share one help page with a description;
  their examples fit two epochs first, so they run on the bundled data
  (an unfitted model has no environment scaling yet).

* Removed from the exported interface, as nothing in the package, the tests
  or the vignettes used them: `aggregate_results_old()`, `np_runif()`,
  `sample_poisson_relaxed()`, and the internal helpers `groupby_mean()`
  and `index_species()` are no longer exported.

## Documentation and data preparation

* The process functions `growth()`, `mortality()`, `regeneration()`,
  `competition()` and `regeneration_capped()`, `summary()` of a model,
  `plot()` of an `ALE()` object, `compute_env_scaling()` and `dbh2ba()` have
  runnable examples in the context of the fitting workflow; the `print()`
  method of `finnControl()` is documented with it. `?fit` says in a few lines
  what makes a fit faster instead of reporting the benchmark behind it, the
  `env_autoscale` item is short and points to `?finn_data`, and the help
  pages no longer cite issue numbers.

* Examples for the fitting workflow (`finn()`, `fit()`, `predict()`,
  `simulateForest()`, `predictTrees()`, `regeneration_capped()`,
  `regeneration_adult()`, gated on an installed torch backend) and for the data
  helpers (`makeObsData()`, `resolveSiteIDs()`, `makeInitCohorts()`,
  `climateDF2array()`); the fitting vignettes gain sections on loss balancing,
  convergence checks, persistence and blocked holdouts, staged calibration, the
  recruitment cap and per-tree mortality scoring.

* New datasets `fia_tree_dt` and `fia_tree_test`: the individual trees behind
  the mortality counts of `fia_obs_dt` and `fia_obs_test` (dbh, growth, fate
  and interval length per tree, all macroplot trees). The FIA vignette uses
  them to score mortality per tree (`loss_aggregation = c(mortality = "tree")`)
  and compares the size-mortality curve of that fit with the site-level fit.

* `resolveSiteIDs()` returns `env_dt` with one row per site and year; it used to
  repeat each row once per patch, which `climateDF2array()` now rejects, so
  `fit()` failed on the output of the documented data-preparation pipeline.

* `fit()` stops when the per-site path (`loss_aggregation = "tree"` or
  `model$period_mode = "per_site"`) meets a `data` without a finite
  `period_length`. Every site's window was then empty and all six loss terms
  were zero, so the fit ran to the end without scoring anything (#46). The
  `data` documentation no longer claims that `period_length` must be the same
  for every site.

* `fit()` gains `threads`, torch's intra-op thread count during the fit. The
  default `NULL` leaves torch's setting alone, so `OMP_NUM_THREADS` or
  `torch_set_num_threads()` keep working. FINN's forward pass is a chain of
  small tensor operations, and extra threads rarely pay (#48):
  - on macOS one thread was up to 30 % faster at the size of a BWI calibration
    batch (512 sites x 25 patches);
  - on a 64-core Linux node the thread count changed little, and only a single
    2000-site batch gained from 4 threads.

  `threads = 1` is the considerate choice when several fits share a machine.
  Where torch uses OpenMP the previous count is restored afterwards. torch's
  native backend (macOS) allows one change per session, so there the session
  keeps the count, and `fit()` says so when an earlier setting prevents the
  change.

## The release in detail

This release adds the model features developed for the **FIA** (US Forest
Inventory and Analysis) and **BWI** (German Bundeswaldinventur) calibration
projects. Every new argument defaults to the previous behaviour, with three
exceptions that do not reproduce 0.1.0 and are marked in the table below:
`mortality(growth_scale = 0.01)`, `createProcess(light_mask = TRUE)` and
`createProcess(reg_floor = 1e-3)`. Two further exceptions repair the stem-count gradient
and are on by default: `fit(trees_offset = 0.01)` and `fit(rb_gradient = TRUE)`,
both described under *The stem-count gradient* below. One further change
reproduces 0.1.0 only on a machine whose `LC_COLLATE` was already `C`: `resolveSiteIDs()` now numbers sites,
patches and species by radix order instead of the collation locale's.

The two tables are the overview; the sections below them carry the detail.

### The stem-count gradient

Two defects in how the `trees` term reached the parameters, both repaired here
and both on by default (FINNverse/FINN#53).

The count likelihood added a constant to the predicted stem count so that
`log_prob` stays finite where a species' simulated count is exactly zero and the
observed one is not. That constant was `1e-5` for `"poisson"` and `0.001` for
`"nbinom"`, which leaves the value finite and the gradient at `y/offset` -- four
orders of magnitude above every other term in the objective. Where the gradients
are clipped per group, one such cell then rescales the whole mechanistic block.
`fit(trees_offset = )` exposes the constant and defaults to `0.01`, on the scale
of the predicted count. `NULL` restores the two historical values.

The simulator draws deaths and recruits and rounds them, and the gradient ran
through both the reparameterised draw and the rounding's straight-through
surrogate, so the sampler's noise reached the parameters. On the German
inventory 97 % of the mortality gradient was draw noise (signal-to-noise 0.52
against 43 without it), which made a fit diverge or not depending on the host's
float kernels. `fit(rb_gradient = TRUE)` keeps the drawn, rounded count in the
forward pass -- a simulation is unchanged -- and takes the derivative through
`n_at_risk * m` and the expected recruitment instead. Against a 200-draw
reference gradient the repaired direction agrees to a cosine of 1.000 in every
parameter block at `update_step = 1`, where the unrepaired objective's cosine is
negative. `FALSE` restores the previous path.

Only fits that score `trees` are affected; a model scoring growth and mortality
alone has no such term.

### The death draw

`mortality_draw` gains `"binomial"`, which is the new default, and `"bernoulli"`.

`"binomial"` is the exact binomial, and it resolves by whether the draw has to
carry a gradient. With `rb_gradient = TRUE` (also the new default) the drawn
count is detached and the derivative comes from the expectation, so the draw is
torch's own binomial: exact in both moments at constant memory, and
non-differentiable, which no longer matters. With `rb_gradient = FALSE` the
gradient reaches the mortality parameters *through the state* -- deaths enter
`trees`, `trees` enters `ba` -- so the draw must be differentiable, and it
resolves to the moment-matched Gamma instead. Measured on that state path, the
Gamma gives -39.36 against the expectation's -39.32, while a straight-through
Bernoulli gives -22.59 and the capped Poisson gives exactly 0.

Consequence worth knowing: under `"binomial"`, switching `rb_gradient` changes
the sampler as well as the gradient, so the two no longer simulate the same
stand. Hold `mortality_draw` fixed to compare them.

`"poisson"` (a Bernoulli for one-tree cohorts, otherwise a capped Poisson) stays
available. It is integer-valued and exact in expectation for small `n*p`, but
its variance is `n*p` rather than `n*p*(1-p)` -- measured 12.044 against 8.400 at
`n = 40`, `p = 0.3` -- and it has no pathwise gradient. Combining it with
`rb_gradient = FALSE` leaves the mortality parameters with exactly zero gradient,
so that combination is now refused rather than silently fitted.

`"bernoulli"` is `binomial_from_bernoulli()`, repaired in this release
(FINNverse/FINN#62): it returned exactly zero for every input, through an
inverted trial mask and a reparameterised draw that was passed to a second
`torch_bernoulli()`. It is exact in both moments and keeps a gradient path, but
expands each cell to the largest cohort in the batch, so its memory grows with
the largest `trees` value in the data.

### The recruit draw

`createProcess(regeneration_sampler = "refined")` becomes the default. It samples exactly below
`lambda = 2` and drops the `+0.001` offset and the `1e-4` clamp, which together
imposed a recruitment floor of 0.0011 per cell-year -- at a 0.0168 ha patch,
0.066 recruits per hectare and year **per species**.
`regeneration_sampler = "gaussian"` restores the previous draw.

### The mortality limits

The simulation clamps the per-cohort mortality probability before the death
draw. That interval was hard-coded to `[0.01, 0.99]`; it is now
`createProcess(mortality_floor =, mortality_ceiling =)` and defaults to
`[0.001, 0.999]`.

The bound is per **simulation step**, not per year, and the old floor is why a
simulated low-mortality stand -- a mature unthinned stand, a slow-growing
species, an equilibrium run -- lost stems the observations did not: at a 1 %
floor every cohort is at risk at 1 % per step however low the fitted rate says
it should be. The per-tree mortality loss was never bounded by it, so a fit
could be correct while the trajectory it simulated was not.

The floor cannot simply be dropped: `binomial_from_gamma()`'s shape
`n*p/(1-p)` degenerates as `p -> 0`. Measured with that draw it stays finite
with an unbiased mean gradient down to `1e-4`, but single-draw gradients grow
heavy-tailed as `n*p` shrinks -- the largest per-tree gradient is about ten
times the old default's at `1e-3` and about a hundred times at `1e-4` -- so
`1e-3` is the lowest tested value and is the new default.

### The recruitment cap

`regeneration_capped()` replaces `regeneration_saturation()` (kept as a
deprecated alias for one release) and becomes the default regeneration process.
`FINN::regeneration`, unbounded as before, stays available and is now the opt-in.

Two things change beyond the name. The cap parameter `reg_logK` is **created by
FINN**, so nothing has to be declared through `custom_parameters` any more;
`createProcess(cap = )` sets its resolution, `"all"` for one cap shared by every
species and `"species"` for one per species.

And the cap carries a prior, because nothing else in the objective constrains
it: a cap fitted freely runs to whatever the calibration window supports, which
on a national inventory is far above anything observed. The penalty is
`prior_lambda * (log K - log prior)^2` per cap, the same form as the L2 pull of
a staged calibration. `prior = "max_obs_recruits"` (the default) resolves it at
the start of the first `fit()` from the largest observed recruitment rate per
timestep, `max(reg / period_length)` -- `reg` is the sum over an observation's
interval while the cap bounds one timestep's mean -- which is also the cap's
starting value; a later `fit()` continues from the trained cap rather than
resetting it; a number sets the prior directly, on the same per-timestep scale; `prior_lambda = 0` switches
it off. With `cap = "species"` the prior is hierarchical -- each species' log cap
is pulled towards the differentiable mean of all of them, so it penalises how far
the species spread apart and leaves the overall level to the data.

The penalty is added **after** the `loss_weights = "auto"` baseline scaling, so
`prior_lambda` means the same thing on every dataset.

### What is new since 0.1.0

| what | where | default |
|---|---|---|
| Offset of the `trees` count likelihood where a simulated count is exactly zero | `fit(trees_offset = )` | `0.01`; **not** 0.1.0 behaviour, `NULL` restores it |
| Exact binomial death draw instead of the moment-matched Gamma | `createProcess(mortality_draw = )`, `createHybrid()` | `"binomial"`; **not** 0.1.0 behaviour, `"gamma"` restores it |
| Exact recruit draw below `lambda = 2` instead of the Gaussian approximation | `createProcess(regeneration_sampler = )`, `createHybrid()` | `"refined"`; **not** 0.1.0 behaviour, `"gaussian"` restores it |
| Lower and upper limit of the mortality probability, per simulation step | `createProcess(mortality_floor =, mortality_ceiling =)`, `createHybrid()` | `0.001` / `0.999`; **not** 0.1.0 behaviour, `0.01` / `0.99` restores it |
| Lower bound on the negative-binomial recruit dispersion, applied to the draw and the likelihood alike | `createProcess(theta_floor = )`, `createHybrid()` | `0.5`; **not** 0.1.0 behaviour, `0` or `NULL` restores it |
| Recruit dispersion kept fixed at `dispersion_parameter` instead of estimated | `createProcess(optimizeTheta = )`, `createHybrid()` | `TRUE` (estimated, as in 0.1.0) |
| Recruitment capped by a Beverton-Holt term, with a prior on the cap | `finn()`'s default regeneration process, `createProcess(cap =, prior =, prior_lambda =)` | `regeneration_capped(cap = "all", prior = "max_obs_recruits", prior_lambda = 0.1)`; **not** 0.1.0 behaviour, `FINN::regeneration` restores it |
| Count gradients taken through the expected deaths and recruits instead of through the draw | `fit(rb_gradient = )` | `TRUE`; **not** 0.1.0 behaviour, `FALSE` restores it |
| Oregon FIA calibration sample as lazy-loaded datasets | `fia_obs_dt`, `fia_env_dt`, `fia_init_trees`, `fia_obs_test`, `fia_env_test`, `fia_init_test`, `fia_species_dt` | — |
| Per-response likelihood, renamed | `fit(loss_family = )` | as `loss` had it; `loss` still accepted, deprecated |
| How a response is reduced before it is scored: `"mean"`, `"sum"`, `"size_classes"`, `"quantiles"`, `"tree"` | `fit(loss_aggregation = )` | `NULL` (per-cell means) |
| Tree table feeding both per-tree routes | `fit(tree_data = )` | `NULL` |
| Score growth and mortality on each simulated cohort that started as an observed tree | `fit(loss_resolution = )` | `"site"` |
| Evaluate the observed-tree route on a fitted model | `predictTrees()` | new export |
| Per census and epoch: `frac_killed`, `mean_kill_year`, `frac_ghost`; and the per-tree predictions | `$tree_deaths`, `$tree_pred` | under `loss_resolution = "tree"` |
| Multiple shooting: rebuild the stand from the observed tree list at every anchor | `fit(reanchor = , anchor_cohorts = , anchor_mode = , anchor_min_dbh = )` | `FALSE`, `NULL`, `"replace"`, `NULL` |
| Reset the stand to the observation every year | `fit(teacher_forcing = )` | `FALSE` |
| Drive the stand with the observed recruitment (`reg / period_length` per ha and year in every year of its interval, drawn by the recruit sampler) instead of the regeneration process, which is frozen and unscored for that call: stage 2 of a staged calibration | `fit(loss_recruit = )`; `simulateForest(recruits =, obs_recruits =)`, `predict()` | `"simulated"`, as in 0.1.0 |
| Per-site aggregation windows for mixed remeasurement intervals | automatic from `period_length`; override with `model$period_mode` | `"auto"` |
| A non-finite loss term drops that term, or the whole timestep; occurrences counted | `fit(nonfinite = )`, `$nonfinite_counts` | `"drop_term"` |
| Set-up only, no training | `fit(epochs = 0)` | — |
| Per-predictor scaling: `"auto"`, `"identity"`, `"0to1"` or a function | `fit(env_autoscale = )` | `TRUE`, which z-scores as before; the new strings are opt-in |
| Recruitment limited by conspecific adults | `regeneration_adult()` | new export |
| Conspecific adult density as an effect on recruitment, with a share `f0` recruited without adults | `regeneration_adult_effect()` | new export |
| A process function declares and calibrates its own named parameters | `createProcess(custom_parameters = )` | `NULL` |
| **Additive recruitment floor, was hard-coded at 0.2.** A default that does **not** reproduce 0.1.0 — see below | `createProcess(reg_floor = )` on the regeneration process | **`1e-3`** (was 0.2) |
| Observation operator mapping predicted and observed recruitment per ha into the tallied-count space of the sampling design (applied under the `"nbinom"` regeneration likelihood only, #65) | `finn(recruit_obs_weight = )` | `1.0` |
| Score growth over the whole remeasurement interval | `finn(growth_rate = )` | `"annual"` |
| **Binned competition takes its height quantiles over the live cohorts only.** A default that does **not** reproduce 0.1.0 — see below | `createProcess(light_mask = )` | **`TRUE`** (was unmasked) |
| A cohort id identifies the initial cohort it came from, so a simulated cohort traces back to its tree | `simulateForest(return_cohorts = TRUE)`, `$pred` | — |
| **A numeric `weights` vector multiplies the baseline-scaled objective instead of replacing it.** A default that does **not** reproduce 0.1.0 | `fit(loss_weights = )`, `fit(weight_scale = )` | **`"auto"`**; `weight_scale = FALSE` restores the old meaning |
| **Growth enters `mortality()` in %/yr.** A default that does **not** reproduce 0.1.0 — existing code changes its answer | `mortality(growth_scale = )`, overridable per model with `createProcess(growth_scale = )` | **`0.01`** (was the raw rate); `growth_scale = 1` restores 0.1.0 |

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

* **`plot()` of a fitted model draws its parameters.** `plot(model)` now runs
  from a bare call and shows the fitted parameters themselves, per species:
  the species-specific parameters of every process and the recruit dispersion
  `theta` on the scale they act on, each against its admissible range from
  `createProcess(lower, upper)`, so a parameter that ended the fit at a bound
  stands out. `pars = "env"` switches to the environmental niches: for every
  process whose environmental response is linear, the fitted environmental
  term over the range each predictor covers in the training data, on the
  predictor's original scale; a process with a neural network response is
  skipped with a warning that points to `ALE()`. `species =` selects species.
  The argument `pars` and its values `"process"` / `"env"` are those of 0.1.0,
  where the method was undocumented, failed on a bare call (a length-2
  condition in `if()`, an error since R 4.2.0) and drew the first-layer
  weights of the environmental network as bars.

* **`plot()` of a simulation.** `simulateForest()` and `predict()` return a
  list of class `finn_simulation` (the tables inside are unchanged), and
  `plot(sim)` draws basal area, mean dbh and stem number per species over the
  simulated years, averaged over the sites, or for one site with
  `sites =`; `response =` takes any of the seven simulated variables. The
  leading species are coloured and named, every other species is a pale grey
  line, and the stand total is dashed. A bare-ground run is drawn as a
  succession, one panel per site labelled with its mean climate.

* **`ALE()` and the conditional effects for any process function.** FINN now
  records the inputs of every process step itself while it computes the
  conditional effects, so a process function of one's own (which does not
  record them) is covered as well; a model saved with an earlier version keeps
  the methods it was saved with, and a process without recorded inputs is left
  out with a message instead of an error. `plot(ale)` names the species.

* **`plot(model, parameter_names = )`** names the parameters of a process
  function of one's own; FINN's own process functions are named from their
  documentation.

* **Species by their names.** `fit()` keeps the names of the `species_name`
  column of the data (`model$species_names`); the plots label species by
  name, and every plot takes `species_names =` for a model without them.

* **`comparePlot()`: simulation against observations.** `comparePlot(sim, obs)`,
  or `plot(sim, obs)`, pairs a simulation with the observed data, aggregated
  the way `fit()` scores it (states at the census, growth and mortality
  compounded over the interval, recruits summed). It draws observed against
  simulated for `dbh`, `ba`, `trees` and `growth` with R2 in each panel, the
  ROC curve of the interval death probability with its AUC for mortality, and
  for recruitment the ratio of simulated to observed recruitment over all
  sites per species. With the observed trees (`tree_data =`, as for `fit()`)
  the dbh panel shows the largest diameter per site and species (all species
  pooled: per site), which is robust against a few large or many small trees;
  the simulation records its cohorts at the census years for it.
  `type = "heatmap"` shows the skill of every species at a glance: R2 (and the
  AUC for mortality) per response and species, next to the ratio of
  simulated to observed. `type = "timeseries"` draws the simulated trajectory
  with the observations as points at their census years. `species =` pools all
  species (default) or draws chosen species; long and wide tables are both
  accepted, and observations without a simulated counterpart are left out with
  a message. An observation is placed as in `fit()`, by the position of its
  `year` among the years of `env` (the simulation keeps them as the attribute
  `"env_years"`), so `env` and the observations may carry calendar years.

* **`evaluate()`: an optional step after the fit.** It simulates the training
  data, a holdout (`holdout = list(data =, env =, init_cohort =)`) and,
  optionally, a bare-ground run (`bare_ground = TRUE`, `years = 500`), and
  returns the comparison as an object of class `finn_evaluation`, which the
  model also keeps (`model$evaluation`) for `summary()`. `print()` shows the
  skill tables; `plot()` draws the comparison (`type = "scatter"` or
  `"heatmap"`, from the stored tables, no new simulation) or the bare-ground
  runs (`type = "bare_ground"`). `fit(..., holdout = )` runs it right after the
  fit. A model fitted with an earlier version takes its training data through
  `train =`.

* **`?FINN` shows how to start.** The package help page lists the functions in
  the order of an analysis (prepare the data, build, fit, evaluate, look
  inside, simulate) with a runnable example, and the README walks through a
  fit of the bundled FIA data with its evaluation and a bare-ground run.

* **`summary()` of a model judges the fit.** It now shows, before the
  environmental importance and the average conditional effects, the skill per
  response on the training data and the holdout (R2 and the correlation of
  observed and simulated first, then RMSE, bias,
  the ratio of the means, and the AUC for mortality) and a convergence check:
  by how many percent each loss term changed over the last `epochs_back`
  epochs, and each parameter: every species-specific parameter by its name on
  the scale it acts on, every network block as a whole, the largest change of
  each named in a sentence (`print(s, epochs_back = 10)`). For this `fit()`
  keeps the parameters of the last 51 epochs; for an older model the parameter
  snapshots of the fit are used. `summary()` returns an
  object of class `summary.finn_class` that `print()` shows, instead of
  printing as it computes.

* **Bare-ground runs in one call.** `simulateForest(model, bare_ground = TRUE)`
  runs the model from an empty stand for `years = 500` years under the mean
  climate of the fitted environment held constant, or under each site's mean
  climate with `sites =`.

* **The patch count of the fit is the default.** `fit()` stores the patch
  count it used (`model$patches`), together with the training observations;
  `predict()`, `simulateForest()` and `evaluate()` take it when `patches` is
  not given.

* **`regeneration_adult_effect()`: conspecific adults as an effect, not a
  gate.** `regeneration_adult()` multiplies the whole potential recruitment by
  `A / (A + A_half)`, so a species without adults on a site recruits only at
  `reg_floor` and bare ground cannot be colonised. The new process uses
  `f(A) = f0 + (1 - f0) A / (A + A_half)` with
  `f0 = plogis(reg_adult_logitf0)`, a new custom parameter (shared or per
  species) for the share of the potential recruitment without conspecific
  adults (seed rain from outside, advance regeneration). `f0 -> 1` is
  `regeneration_capped()`, `f0 -> 0` is `regeneration_adult()`, which is
  unchanged. See `?regeneration_adult_effect` and `?recruitment`.

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

* **Two repairs to the stochastic draws
  (`createProcess(mortality_draw = , regeneration_sampler = )`).** Both are on by
  default in 0.2.0: `mortality_draw = "binomial"` and
  `regeneration_sampler = "refined"`. `mortality_draw = "gamma"` and
  `regeneration_sampler = "gaussian"` restore the draws of 0.1.0, which are
  bit-identical to it (`tests/testthat/test-identity.R`, config `sample_reg`).
  The integer-valued alternative `mortality_draw = "poisson"` replaces the Gamma matched to
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

* **Recruitment observation weight (`finn(recruit_obs_weight = )`).** An
  observation operator, not a weight of the loss term: it multiplies predicted
  and observed recruitment per ha alike before the regeneration likelihood, so
  that the likelihood scores the counts a sampling design actually tallies. In
  0.2.0 it is applied under `loss_family["regeneration"] = "nbinom"` only (#65).

* **Period-scale growth scoring (`finn(growth_rate = "period")`).** Scores growth
  over a multi-year remeasurement interval rather than a single step, for
  inventories whose growth is only observed between remeasurements.

* **Custom process parameters (`createProcess(custom_parameters = )`).** Lets a
  user-supplied process function declare and calibrate its own named parameters,
  beyond FINN's built-in per-species and per-environment coefficients.

## Behaviour changes

* **Calibration runs on annual simulation steps: data preparation, bundled
  data and vignettes.** One simulation step per ten-year remeasurement
  interval distorts the dynamics: recruitment, competition and growth act
  once per decade, and the per-step mortality floor and ceiling become limits
  per decade. On FIA Oregon a model calibrated that way grew to a bare-ground
  equilibrium of 100-170 m2/ha that was still rising; calibrated on annual
  steps it levels off at 45-60 m2/ha. So everything the package ships now
  uses one step per year:
  - `makeObsData(rates = c("annual", "interval"))`: `growth` and `mort` are
    yearly rates by default, `(dbh / dbh_before)^(1 / p) - 1` per tree and
    `1 - (1 - n_died / n_at_risk)^(1 / p)` per cohort over an interval of `p`
    years, and the table carries `period_length` (years since the previous
    inventory). `n_at_risk`, `n_died`, `growth_n` and `reg` stay quantities of
    the interval, as `fit()` uses them: it averages the simulated yearly
    growth and mortality over the `period_length` years before a census,
    scores `mort` with the binomial weighted by `n_at_risk`, and sums the
    simulated recruits. The returned tree list gains `growth`, the tree's rate
    in the same unit (`rel_growth` stays the interval increment).
    `rates = "interval"` returns the rates over the whole interval, as before.
  - `resolveSiteIDs(time_axis = c("annual", "census"))`: `year` is the number
    of years since the first inventory (0, 10, 20 rather than 0, 1, 2), the
    first inventory is the initial state (initial cohorts, no observation
    row), and `env_dt` comes back with one row per site and year `1..Tmax`,
    so it needs one row per calendar year (a static climate is repeated).
    `env_dt` now holds only `siteID`, `year` and the predictors, which is
    what `fit(env = )` accepts (it carried the character `siteName`, which
    `fit()` refuses). `time_axis = "census"` keeps the inventory index of
    0.1.0, but builds the initial cohorts from the first inventory, as the
    annual axis does; 0.1.0 took the second inventory (`year == 1`) for them,
    so 0.1.0 results are not reproduced exactly on the census axis either.
  - The bundled FIA datasets (`fia_obs_dt`, `fia_env_dt`, `fia_obs_test`,
    `fia_env_test`) are rebuilt on the annual axis from the same plots and
    sample: observations at years 10 and 20 (21 and 22 for two sites whose
    third inventory came late) with `period_length`, yearly growth and
    mortality, and `env` over years 1..22. Growth is now also tree-weighted
    across patches as `makeObsData()` defines it; the old table carried the
    unweighted mean of the patch means. `ba`, `dbh`, `trees`, `reg` and the
    counts are unchanged; `example_env_dt.csv` has one row per calendar year.
  - The vignettes, examples and tests that use these data simulate the 20
    years between the inventories. The vignettes score a fitted model the
    way `fit()` does: the stand at the census year, growth and mortality
    averaged over the `period_length` years before it, and the expected
    recruits per hectare (`r_mean_ha`) summed over them.
  Results computed on the bundled data or with `makeObsData()` /
  `resolveSiteIDs()` defaults change.

* **The simulated rates of an interval are compounded as the observed rates
  are.** `makeObsData(rates = "annual")` derives a yearly growth rate as
  `(dbh / dbh_before)^(1 / p) - 1` and a yearly mortality as
  `1 - (1 - d / n)^(1 / p)` over an interval of `p` years, but `fit()` compared
  them with the arithmetic mean of the simulated yearly rates over the
  `period_length` window. With one disturbance year in a window the two differ
  a lot: a mortality of 0.5 in one year and 0.01 in nine is observed as
  0.0754, and the mean is 0.059. `fit()` now compounds the window the same
  way, growth as `exp(mean(log1p(g))) - 1` and mortality as `1 -
  exp(mean(log1p(-m)))`, with the same masks as before: a year outside a
  site's window or without trees counts as growth 0 and survival 1, and the
  denominator is the full `period_length`, since the observed rate is derived
  from the counts over the whole interval. Mortality is clamped below `1 -
  1e-4` before the log (sigmoid saturates to exactly 1 in float32 above a
  logit of 17.3), growth is not clamped (the package's growth is non-negative).
  The per-tree route follows for growth (`loss_aggregation = "tree"`,
  `loss_resolution = "tree"`, `predictTrees()`); its mortality was already the
  compounded interval probability against the tree's 0/1 fate. The constant
  increment of `finn(growth_rate = "period")` is unchanged, as is the census
  axis, where one step is one interval. Vignettes D and E score the fitted
  models with the same functional. Fits on data with `period_length > 1`
  change; the identity snapshots of the test suite were regenerated.

* **Mortality starts at about 2 % per year.** The environmental response of
  `mortality()` started at torch's random weights, i.e. at `sigmoid(0) = 0.5`
  per step at the average environment, so on yearly steps an unfitted stand
  lost most of its basal area within ten years and the first epochs of a fit
  began from a dying stand. Without `createProcess(initEnv = )` the response
  now starts with an intercept of `qlogis(0.02)` and every other weight 0
  (with `hidden` layers the output bias carries the intercept); the species
  terms for light, size and growth move the start by at most about one logit
  unit. A formula without an intercept (`~ 0 + temp`) keeps torch's start.
  `?finn`, `?createProcess` and `?mortality` say so. Fits and simulations
  that did not set `initEnv` for mortality start from another point and so
  change.

* **Bare ground is bare: `init_cohort = NULL` starts with no stems.**
  `simulateForest()` / `predict()` and `fit()` without an `init_cohort` used to
  start every patch with one 1-cm stem of every species, although the help said
  "empty patches". That planted species the regeneration process would not
  recruit, so a bare-ground (equilibrium) run did not test the regeneration
  process: on FIA conterminous US a species of moist eastern sites grew to
  62-73 % of 112-168 m2/ha at dry western sites where it barely recruits; from
  an empty start the same model ends at 30-61 m2/ha. Now every tree of such a
  run comes from the regeneration process. Results of runs without an
  `init_cohort` change; runs with one are unaffected.

* **Trees can grow up to 130 m, not 101 m.** The upper bound of the fitted height
  parameter (first column of the competition parameters) is `log(2.3)` instead of
  0.7. `height()` approaches `(exp(parHeight) - 1) * 100` m as the diameter grows,
  so 0.7 capped every species at 101 m, below the 122-130 m physical limit of tree
  height (Koch et al. 2004), and a species given exactly 0.7 (as in the
  introduction vignette) was clamped with a warning. Fits that train the
  competition parameters change slightly, because the bound sets the scale on
  which the parameter is optimised. `createProcess(upper = c(0.7, 2.0))` on the
  competition process restores 0.1.0.

* **The recruitment floor is `1e-3`, not 0.2.** `createProcess(reg_floor = )` replaces a
  hard-coded additive floor of 0.2 on the expected recruitment of every species
  in every cell-year. At 0.2 the floor alone delivered about 13 stems per hectare
  per year over a 65-species FIA pool, against roughly 0.7 observed -- the floor
  was the recruitment. It exists to keep the negative-binomial likelihood away
  from a zero mean, which `1e-3` does. `createProcess(reg_floor = 0.2)` on the regeneration process restores 0.1.0.

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

* **A numeric `weights` vector multiplies the scaled objective, it no longer
  replaces it.** `fit(weights = )` used to bypass the `"auto"` baseline scaling
  entirely, so the six numbers acted on losses whose raw magnitudes differ by
  three orders of magnitude. Measured on BCI with the weight vector of the
  published configuration, `c(0.1, 10, 1, 1, 1, 1)`: the raw baseline MSEs stand
  at dbh : ba : growth = 1246 : 26 : 1, so growth held 0.26 % of the objective
  at the start and 4.5 % at convergence, and `ba`'s nominal factor of ten was a
  factor of ten on a scale already 26x growth's. Weights now multiply the
  baseline-scaled losses, where a weight reads directly as that response's
  share. `fit(weight_scale = FALSE)` restores the old behaviour exactly and is
  what reproduces an analysis fitted before 0.2.0. Every fit now prints each
  response's share of the objective at the intercept-only baseline, in both
  modes -- the silence was the actual defect.

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
  also what the `loss_weights = "auto"` baselines measure; the reported history is
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

* **Fewer torch calls per time step in `forward()`.** On the tensors of one
  batch (~30k elements) a torch op costs ~25 us of fixed overhead and the time
  loop made ~2000 of them per step, so the op count, not the arithmetic, sets
  the speed. The per-species aggregation of the rates, the recruits and the
  stand (growth, mortality, dbh, basal area, trees, regeneration) now scatters
  all its samples in one call (`aggregate_species()`, internal) and writes the
  masked tree-weighted means with `torch_where()` instead of two masked index
  assignments; `trees$sum() > 0` is read with `$item()`; whether the gradients
  are yearly (no `period_length` anywhere) is decided once per call instead of
  by a reduction over `y` at every step. The results are bit-identical
  (checked on two BWI stage-2 calls: parameters, loss history, baselines,
  `param_history`, `predictTrees()`).

* **The per-tree scoring evaluates the window in one pass.**
  `loss_aggregation = "tree"` and `predictTrees()` applied growth and mortality
  to the observed trees once per year of the remeasurement window (10 calls of
  each process per batch on the BWI data); the years are now a dimension of
  the tensors and the processes run once (mechanistic processes; hybrid ones
  keep the year loop). Same formula, but the sums over the years and over the
  trees are taken in another order, so the result is equal only to floating
  point precision: on two BWI stage-2 calls the fitted parameters after 3
  epochs differ by at most 2e-6 (absolute, float32), the loss history by 6e-8
  and `predictTrees()` by 3e-6 (FINN-bwi fits/2026-11-06-finn-speed-review).

## Bug fixes

* Initial cohorts survive `saveRDS()` / `readRDS()` (#77). A `CohortMat` (e.g.
  `initCohorts` of the list `resolveSiteIDs()` returns) keeps a plain R copy of
  its arrays and rebuilds its torch tensors, which `saveRDS()` cannot store,
  when `fit()`, `predict()` or `simulateForest()` use it, or with
  `x$restore()`. Before, every tensor failed with "external pointer is not
  valid" after reading.

* An environment with no covariate column gives the process an intercept
  instead of an all-`NA` design (#64): `model.matrix()` on a frame with no
  columns returned zero rows, and the recruitment of an intercept-only
  process became `NaN`.

* `createProcess(initSpecies = , initEnv = )` accept a data.frame and use the
  matrix it holds, instead of failing with "Can't convert argument" (#11).

* `fit()` no longer creates or overwrites an object `counter` in the global
  environment. The batch counter was advanced with `<<-`, which also left the
  local counter at 1: with several batches per epoch (`batchsize` smaller than
  the number of sites) `$history`, the progress plot and the plateau scheduler
  saw the last batch's loss instead of the mean over the epoch's batches. With
  one batch per epoch, the default, nothing changes.

* A hybrid growth process is given the growth design matrix and a hybrid
  mortality process the mortality one. The environment list is ordered
  mortality, growth, regeneration, and the two were read by position, so
  `createHybrid()` growth was fitted and simulated on the mortality formula and
  vice versa. Regeneration and `predictTrees()` always read theirs by name, so
  the two routes disagreed for the same model. Nothing changes where the growth
  and mortality formulas are the same, which is every example in the package;
  formulas of different width stopped with "mat1 and mat2 shapes cannot be
  multiplied".

* `createHybrid(hidden = )` accepts a single hidden layer. `hidden = 16L` used
  to build the layer stack with `2:length(hidden)`, which counts backwards from
  2 to 1, and the first `fit()` stopped with "Expecting a single value:
  [extent=0]". The default `hidden = c(50L, 50L)` is unaffected.

* `predictTrees()` checks that `tree_data` and `env` hold the same sites, that
  every observed year is in `env` and that the species ids are those of the
  model. The tree tensors follow the sites of `tree_data` and the environment
  array those of `env`, and torch's gather accepts the shorter index, so a tree
  could be scored against another site's climate without any message.

* `predictTrees()` on a model that has never been fitted or simulated asks for
  `patch_size` instead of ending the R session. `BA_stand()` divided by a `NULL`
  patch size, which segfaults `torch_div()`.

* `conditionalEffects()`, `ALE()` and `summary()` pair each stand state with its
  own site's and year's predictors when `env` is not sorted by `siteID` and
  `year`. The recorded states follow the simulation's sorted order and the
  predictors were taken in order of appearance.

* `conditionalEffects()` and `ALE()` recompute when `env` or `init_cohort` is
  passed explicitly. The cache key is a fingerprint (row and column count,
  column names, column sums), so two different inputs could share it and the
  first call's result was returned for the second.

* `createProcess()` and `createHybrid()` find a formula that is held in a
  variable of the calling function (`f <- ~ temp; createProcess(f, growth)`
  inside a function used to stop with "object 'f' not found").

* `predictTrees()` puts the model in evaluation mode and on the scoring device
  before it predicts, as `simulateForest()` does. Straight after `fit()` the
  model was still in training mode, so a process with `dropout` returned
  different predictions on every call.

* `makeInitCohorts(dbh_binsize = )` adds up `trees` with `na.rm = TRUE`. A
  record without a stem count (`trees = NA`, e.g. a dead tree) used to turn its
  whole bin into `NA`, which lost the living trees of that bin (349 of the 5350
  stems of `fia_init_trees` with `dbh_binsize = 5`).

* `mortality(debug = TRUE)` returns the list of its components (`light`,
  `environment`, `predM`) instead of stopping with "object 'shade' not found".

* `predict()` and `simulateForest()` left `requires_grad = TRUE` on every
  parameter of the model they had just simulated, so a model frozen for a staged
  calibration silently unfroze as soon as it was simulated (#55). Both now
  restore the state they found.

* `createHybrid(NN = )` was documented as replacing the network architecture but
  was never stored, so a custom module was accepted and then ignored (#59). A
  non-`NULL` `NN` now raises an error; shape the default network with `hidden`,
  `dropout` and `transformer`.

* **A stand that dies out before a census is now scored.** When every cohort
  of a batch had died, the padding pass left no cohort slot and the census
  step was not scored at all (the scoring required at least one slot), so a
  die-off before a census went unpenalised. An empty stand now predicts zero
  basal area, zero trees and the recruits of the window against the
  observations; `dbh`, `growth` and `mortality` have no prediction without a
  tree and stay unscored at such a step. The zero predictions carry no
  gradient (nothing is left to move), but the loss is reported, so a fit that
  kills its stands shows it in `$history` instead of a flat zero.

* A model that never went through `fit()` simulated with the 0.1.0 Gamma death
  draw instead of the `mortality_draw = "binomial"` default, because the flag
  that selects the exact binomial (`rb_gradient`) was set only inside `fit()`.
  `finn()` now sets it to `fit()`'s default, so `simulateForest()` on a freshly
  built model, and a model rebuilt from saved parameters, draw deaths as a fitted
  model does (#70).

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

* `loss_weights = "auto"` no longer explodes on a response that is constant over the
  scored cells. The intercept-only baseline is then a few float ulps rather than
  zero, so `1/baseline` passed the `<= 0` guard and produced weights of order
  1e17 on a package fixture. Such a response is now left unscaled, with a
  warning naming it.

* `fit()` resets `$nonfinite_counts`, so the end-of-fit warning reports this
  fit's dropped terms rather than the total since the model was built.

* Arguments that gate on `isTRUE()` are validated as logicals: `fit(reanchor = )`,
  `fit(teacher_forcing = )` and `createProcess(light_mask = )` accepted `1` and
  then silently behaved as `FALSE`. `createProcess(reg_floor = )` is validated the way
  `growth_scale` is. `fit(loss_aggregation = )` rejects a response named twice
  instead of keeping the first entry, and `fit(env_autoscale = )` rejects a named
  entry that matches no predictor instead of silently scaling that predictor
  with the default.

* **`mortality()`'s growth coefficient is initialised on the scale growth enters
  on.** `parMort[, 3]` was drawn from `U(-0.2, 0.2)` regardless of
  `growth_scale`, so at the 0.01 default the same draw acted a hundred times more
  strongly. Measured on the four-site FIA fixture: an anchored stand of 1667
  stems/ha at 25 cm drops to 1.0-7.5 stems/ha in a single year (mortality 0.93 to
  1.00), against 28-75 stems/ha at `growth_scale = 1`. The range now scales with
  `growth_scale`; a custom mortality kernel without that formal keeps the old
  range. Random initialisations therefore differ from earlier versions'.

* The identity fixtures gained a configuration with unequal scored units per
  observation year (`per_site_unbal`). Every existing per-site configuration
  carried the same number of cells in every year -- 6 against 6 -- and with equal
  counts the per-unit pooling of the per-site loss and the former sum of per-year
  means agree exactly, so no snapshot could see a change to that weighting. The
  new configuration is one species-census short in the first scored year, 5
  against 6.

* `fit(folder = , checkpoints = )` writes each checkpoint to a temporary name and
  renames it. A kill or a full disk during the save used to leave a truncated
  `.pt` that `torch_load()` could not read.

* `fit()` warns when the regeneration process is being calibrated while
  `loss_family["regeneration"]` is `"none"`. No other response can train it:
  recruits reach `dbh` through a sigmoid of scale `1e-3`, whose derivative is
  zero almost everywhere, and `trees` through a rounding whose forward value is 0
  below 0.5 recruits per patch and species.

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

* `pred2DF()` selects its index columns with `with = FALSE` instead of the
  `..dim_names` prefix, which `R CMD check` reported as an undefined global.

* `fit()` and `predict()` check that `data`, `env` and `init_cohort` hold the
  same sites and stop with the offending ids and the table they are missing
  from (#33). Before, a mismatch ended in a torch index error that named
  nothing, or silently paired a site with another site's environment. The ids
  need not be `1:N`: every array is indexed by a site's position among the
  sorted ids. `siteID` must be integer-valued in `data` and `env` (the error
  points to `resolveSiteIDs()`). Covered by `tests/testthat/test-site-sets.R`.

* Site ids that are not `1:N` keep their identity end to end (#34).
  `obsDF2arrays()` numbers sites in radix order (independent of the locale,
  as `resolveSiteIDs()`), says so in a message instead of the former warning
  when the ids are not `1:N`, and records the user ids
  as attribute `"site_ids"`, which `CohortMat` keeps as `$site_ids`; a fitted
  model keeps the lookup from array position to `siteID` as `$site_lookup`;
  and `predict()` / `simulateForest()` return the caller's own `siteID`s
  instead of positions. For ids that are `1:N` nothing changes. Patch and
  cohort ids are renumbered in radix order as well. Covered by
  `tests/testthat/test-site-sets.R`.

* `climateDF2array()` stops on any (`siteID`, `year`[, `month`, `day`])
  combination that appears more than once, identical or not, naming it (#35);
  the last row used to win silently. A missing combination is still `NA`.
  Covered by `tests/testthat/test-climateDF2array.R`.

* The bundled `fia_env_dt` and `fia_env_test` hold one row per site x year
  (600 rows each). They repeated every row four times, once per patch of the
  source plot, all identical, so every environment array built from them is
  unchanged (`dev/make_extdata.R` now deduplicates; the two `.rda` files were
  deduplicated in place). The `env_autoscale` constants do move: the standard
  deviation (n - 1 denominator) over 600 instead of 2400 rows is 0.063 % larger,
  and a seeded fit on the z-scaled climate amplifies that, so the re-knitted
  vignettes D and E report slightly different fits.

* `feature_importance()` permutes over distinct `env` rows. On an `env` with
  repeated site-year rows the copies received different permuted values, of
  which the environment array kept one at random, so importance values computed
  on such data change.

* `createProcess(custom_parameters = )` checks names against a fixed list of
  fields FINN sets itself (#42). A name that a later process sets, such as
  `light_mask`, used to pass, became a trainable nobody read, and turned the
  later assignment into a silent no-op. Covered by
  `tests/testthat/test-processes.R`.

## Documentation

* New help topic `?recruitment`: the recruitment processes, the recruit draw,
  the dispersion `theta` and its floor, the regeneration likelihood and
  `recruit_obs_weight`, feeding observed recruitment, and a guide to which
  options to combine; the recruitment arguments of `createProcess()`,
  `createHybrid()`, `finn()`, `fit()` and `simulateForest()` point to it.


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
* `loss_weights = "auto"` (the default in `fit()`) scales each loss by its
  intercept-only baseline, so the six responses are commensurable and every term
  reads as a fraction of its own null deviance.
* Data helpers `makeObsData()`, `resolveSiteIDs()` and `makeInitCohorts()` turn a
  raw tree list into FINN's input tables.
* Model interpretation with `ALE()`, `summary()`, `feature_importance()` and
  `conditionalEffects()`.

See the vignettes for a guided tour: *Introduction to FINN*, *Plausible
succession from a handful of species*, *Preparing your data for FINN*, *Fitting
FINN to forest inventory data*, and *Mortality*.
