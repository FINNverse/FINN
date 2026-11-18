# Estimate the parameters of a FINN model

Estimates the parameters of a model built with
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md) from
forest inventory data. The model simulates every site from its initial
trees through the years of `env`, the simulation is compared with the
observed stand variables and demographic rates in `data`, and the
parameters of all processes – the species parameters of the mechanistic
equations (for a process built with
`createProcess(optimizeSpecies = TRUE)`), their responses to the
environment, and the weights of any neural-network process – are
adjusted together, in one joint optimization, until the simulation
matches the observations as well as it can. The processes are never
calibrated one at a time: this joint calibration is what makes the
estimated processes consistent with each other.

The defaults are robust in most cases. A first call needs `data`, `env`,
`init_cohort`, `patch_size` and `epochs`; everything else can stay as it
is. [finn_data](https://finnverse.github.io/FINN/reference/finn_data.md)
explains the input tables,
[finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md) how
the simulation is compared with the observations.

## Usage

``` r
fit(
  model,
  data = NULL,
  env,
  disturbance = NULL,
  patches = 100L,
  patch_size = 0.1,
  init_cohort = NULL,
  epochs = 20L,
  lr = 0.01,
  batchsize = NULL,
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality =
    "binomial", regeneration = "nbinom"),
  loss_aggregation = NULL,
  loss_resolution = c("site", "tree"),
  loss_weights = "auto",
  weight_scale = TRUE,
  loss_recruit = c("simulated", "observed"),
  tree_data = NULL,
  control = finnControl(),
  device = c("cpu", "gpu"),
  threads = NULL,
  update_step = 1L,
  start_time = 1L,
  plot_progress = TRUE,
  folder = NULL,
  checkpoints = 100L,
  shuffle = TRUE,
  record_gradients = FALSE,
  env_autoscale = TRUE,
  nonfinite = c("drop_term", "skip_step"),
  reanchor = FALSE,
  anchor_cohorts = NULL,
  anchor_mode = "replace",
  anchor_min_dbh = NULL,
  teacher_forcing = FALSE,
  rb_gradient = TRUE,
  trees_offset = 0.01,
  ...,
  loss = NULL
)
```

## Arguments

- model:

  (`finn_class`)  
  The model built by
  [`finn()`](https://finnverse.github.io/FINN/reference/finn.md).

- data:

  (`data.frame`)  
  The observations: one row per site, census year and species, with the
  columns `siteID`, `year`, `species`, `dbh`, `ba`, `trees`, `growth`,
  `mort` and `reg`; optionally `period_length`, `n_at_risk` and
  `growth_n`. `NA` marks a response that was not observed and is left
  out of the comparison.
  [`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
  and
  [`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
  build this table from a tree list. The columns, their units, missing
  values and incomplete series are explained in
  [finn_data](https://finnverse.github.io/FINN/reference/finn_data.md),
  with the bundled
  [fia_obs_dt](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
  as the example.

- env:

  (`data.frame`)  
  The environment of every site and simulated year: one row per `siteID`
  and `year`, plus one column per predictor (for the bundled FIA data
  `temp` and `prec`, annual mean temperature and precipitation; see
  [fia_env_dt](https://finnverse.github.io/FINN/reference/fia_env_dt.md)).
  The rows are the simulated time steps, so `env` must cover every year
  from the first step to the last observation; the predictors are in
  natural units and standardised internally (`env_autoscale`). See
  [finn_data](https://finnverse.github.io/FINN/reference/finn_data.md).

- disturbance:

  (`data.frame` or `NULL`)  
  Optional disturbances: one row per `siteID` and `year` with the column
  `intensity`, the probability (0 to 1) that a patch of that site loses
  all its trees in that year. Must have the same sites and years as
  `env`. `NULL` (default): no disturbance. See
  [finn_data](https://finnverse.github.io/FINN/reference/finn_data.md).

- patches:

  (`integer(1)`)  
  Number of patches per site when the fit starts from bare ground
  (`init_cohort = NULL`). Ignored, with a warning when it differs, as
  soon as `init_cohort` is given: the patch count then comes from the
  initial cohorts.

- patch_size:

  (`numeric(1)`)  
  Area of one patch in ha. Default `0.1`. Use the plot size of the
  inventory (the bundled FIA data: `0.06`). A model that was fitted or
  simulated with another patch size before gives a warning.

- init_cohort:

  (`CohortMat`)  
  The trees of every site at the first census, from
  [`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md).
  `NULL` (default) starts every site from bare ground with `patches`
  empty patches (no stems). See
  [finn_data](https://finnverse.github.io/FINN/reference/finn_data.md).

- epochs:

  (`integer(1)`)  
  Number of passes over the data. Default `20`; a real calibration takes
  hundreds, until the parameters stop moving. `epochs = 0` runs the
  set-up only (data tensors, `loss_baseline`, `loss_weights`) and trains
  nothing.

- lr:

  (`numeric(1)`)  
  Learning rate, the size of a parameter step. Default `0.01`. Lower it
  when the losses oscillate; see also
  [`finnControl()`](https://finnverse.github.io/FINN/reference/finnControl.md)
  for a learning-rate schedule.

- batchsize:

  (`integer(1)`)  
  Number of sites simulated between two parameter updates. `NULL`
  (default) takes all sites at once. Smaller batches use less memory and
  update more often per epoch.

- loss_family:

  (named `character`)  
  The likelihood each response is compared with: `dbh` and `ba` by
  squared error (`"mse"`), `trees` by a Poisson, `growth` by squared
  error, `mortality` by a binomial, recruitment (`regeneration`) by a
  negative binomial. The defaults are robust in most cases; name only
  what departs from them, e.g. `c(trees = "nbinom")`, or `"none"` to
  leave a response out. See
  [finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md).

- loss_aggregation:

  (named `character` or `list`)  
  How the simulated cohorts are summarised before they meet an
  observation: the default is the mean over cohorts for `dbh`, `growth`
  and `mortality` and the sum for `ba`, `trees` and `regeneration`,
  which matches the tables
  [`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
  builds. The alternatives (`"size_classes"`, `"quantiles"`, `"tree"`)
  are described in
  [finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md).

- loss_resolution:

  (`character(1)`)  
  Whether `growth` and `mortality` are scored per site, year and species
  (`"site"`, default) or per individual simulated tree (`"tree"`, which
  needs `tree_data` with a `cohortID` column). See
  [finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md).

- loss_weights:

  (`"auto"` or `numeric(6)`)  
  The weight of each response in the sum of the six losses, in the order
  `dbh`, `ba`, `trees`, `growth`, `mortality`, `regeneration`. `"auto"`
  (default) makes the six comparable by dividing each by its
  intercept-only baseline, so no response drowns the others; a numeric
  vector then multiplies these balanced terms (`c(1, 1, 1, 2, 1, 1)`
  gives `growth` twice the weight). Before 0.2.0 the argument was called
  `weights`, which is still accepted with a warning. See
  [finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md).

- weight_scale:

  (`logical(1)` or `character(1)`)  
  Whether a numeric `loss_weights` vector applies to the
  baseline-balanced losses (`TRUE` or `"auto"`, default) or to the raw,
  unbalanced losses (`FALSE` or `"raw"`, the meaning of `weights` before
  0.2.0). See
  [finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md).

- loss_recruit:

  (`character(1)`)  
  Where the recruits of the simulated stand come from during the fit.
  `"simulated"` (default): from the regeneration process, which is
  estimated with the rest. `"observed"`: the observed recruitment (`reg`
  in `data`) is fed into the simulation, the regeneration process is
  frozen and left out of the loss, and growth and mortality are
  estimated on a stand that recruits what was observed – the first stage
  of a staged calibration, see
  [recruitment](https://finnverse.github.io/FINN/reference/recruitment.md).
  Before 0.2.0 the argument was called `recruits`, which is still
  accepted with a warning. See
  [finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md).

- tree_data:

  (`data.frame` or `NULL`)  
  Individual observed trees, one row per tree and census, only needed to
  score growth or mortality per tree (`loss_aggregation = "tree"` or
  `loss_resolution = "tree"`). Unlike `data`, which holds one row per
  site, year and species, it holds the trees themselves. Columns and the
  two per-tree routes:
  [finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md)
  and
  [finn_data](https://finnverse.github.io/FINN/reference/finn_data.md).

- control:

  (`finnControl`)  
  Settings of the optimizer: the optimizer itself, a learning-rate
  schedule and its parameters, gradient clipping. The defaults of
  [`finnControl()`](https://finnverse.github.io/FINN/reference/finnControl.md)
  are robust in most cases.

- device:

  (`character(1)`)  
  `"cpu"` (default) or `"gpu"`; `"cuda"` is accepted as a synonym for
  `"gpu"`. Only NVIDIA cards are supported, and a card does not always
  pay: see the section on speed below.

- threads:

  (`integer(1)` or `NULL`)  
  torch's intra-op thread count during the fit. `NULL` (the default)
  leaves torch's setting alone, including one made through
  `OMP_NUM_THREADS` or
  [`torch::torch_set_num_threads()`](https://torch.mlverse.org/docs/reference/threads.html).
  Extra threads rarely pay; `threads = 1` is the considerate choice when
  several fits share a machine. The previous count is restored
  afterwards where torch allows it (OpenMP builds, e.g. Linux). torch's
  native backend (macOS) accepts one change per session, so there the
  session keeps the count `fit()` set, and a count fixed earlier in the
  session cannot be changed (`fit()` then says so).

- update_step:

  (`integer(1)`)  
  Every how many simulated years the gradient is cut (truncated
  backpropagation through time). Default `1`: a loss sees the processes
  of its own year only. The parameters are updated once per batch, with
  the gradients of all cuts added up. A loss is only scored at a year
  that is a multiple of `update_step`, so an observation year off that
  grid is not scored, and `fit()` stops when none is on it. The
  simulated state runs on uncut. Larger values let the loss see a longer
  stretch of the trajectory at a higher cost and with fading gradients.

- start_time:

  (`integer(1)`)  
  Reserved for a burn-in: the simulated year from which observations are
  scored. In 0.2.0 it has no effect – every observation year in `data`
  is scored – and a value other than `1` gives a warning. Start `env` at
  the first year you want simulated and leave out of `data` what should
  not be scored; see
  [finn_data](https://finnverse.github.io/FINN/reference/finn_data.md).

- plot_progress:

  (`logical(1)`)  
  Draw the losses of every response per epoch while fitting, each as a
  fraction of its value in the first epoch. Default `TRUE`.

- folder:

  (`character(1)` or `NULL`)  
  Folder in which the model is saved every `checkpoints` epochs
  (`epoch_<n>model.pt`, readable with
  [`torch::torch_load()`](https://torch.mlverse.org/docs/reference/torch_load.html)).
  `NULL` (default): no checkpoints.

- checkpoints:

  (`integer(1)`)  
  Every how many epochs the parameters are recorded in
  `model$param_history` and, with `folder`, the model is saved. Default
  `100`.

- shuffle:

  (`logical(1)`)  
  Draw the batches in a new random order every epoch. Default `TRUE`.

- record_gradients:

  (`logical(1)`)  
  Also record the gradients of every parameter at each checkpoint, in
  `model$gradients`. Default `FALSE`; the record can get large.

- env_autoscale:

  (`logical(1)` or `character`)  
  How the predictors in `env` are standardised before they enter the
  processes. `TRUE` (default): every predictor to mean 0 and standard
  deviation 1. `FALSE`: as supplied. Or one entry per predictor (named,
  or in column order): `"auto"`, `"identity"`, `"0to1"` or a function
  `f(x)` returning `list(center =, scale =)`, e.g.
  `c(temp = "auto", management = "identity")` to keep an ordinal
  predictor on its own scale. The constants are stored on the model
  (`model$env_scaling`) and re-applied by
  [`predict()`](https://rdrr.io/r/stats/predict.html)/[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md),
  so `env` is always passed in natural units. See
  [finn_data](https://finnverse.github.io/FINN/reference/finn_data.md).

- nonfinite:

  (`character(1)`)  
  What happens when one of the six loss terms is not finite at a scored
  timestep (an overflowing likelihood, a target outside the family's
  support). `"drop_term"` (default) leaves only that term out of the
  timestep's gradient and reported loss; the other responses are still
  backpropagated there, and the occurrences are counted per response in
  `model$nonfinite_counts` and summarised in a warning when the fit
  ends. `"skip_step"` discards the WHOLE timestep – none of the six
  terms is backpropagated or reported there – which is what every
  version before `drop_term` existed did; use it to reproduce such a
  fit. The difference only shows once a term overflows: with finite
  losses the two are identical.

- reanchor:

  (`logical(1)`)  
  Fit by **multiple shooting** instead of single shooting. With the
  default `FALSE` the cohort state is initialised once from
  `init_cohort` and then runs free over the whole horizon, so simulation
  error compounds across every observation. With `TRUE` the simulated
  state is discarded at every timestep for which `anchor_cohorts`
  supplies an observed tree list and rebuilt from it, after that
  timestep's loss has been backpropagated – drift is then bounded by one
  inter-census segment. This is not the same as the gradient truncation
  controlled by `update_step`, which cuts the gradient but keeps the
  state; set `update_step` to the segment length to let the gradient
  span the whole segment.

- anchor_cohorts:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  Named list of
  [CohortMat](https://finnverse.github.io/FINN/reference/CohortMat.md)
  objects, the observed state to restart from. Names are the simulation
  timesteps to re-anchor at (integers as characters, and entries of the
  model's internal `year_sequence`, i.e. timesteps that carry an
  observation). Every anchor must have the same number of species, sites
  and patches as `init_cohort`; the number of cohorts may differ between
  anchors and from `init_cohort`. The last observation year needs no
  anchor – nothing is simulated after it – and an anchor may not fall
  inside another observation's `period_length` window, because the
  aggregated rates would then mix pre- and post-anchor state. `fit()`
  copies the anchors onto its own device and never modifies them, so one
  list can be reused across fits. Ignored unless `reanchor = TRUE`.

  **Build the arrays directly** –
  `CohortMat(dbh = , trees = , species = , dims = , sp = )` – rather
  than via `CohortMat(obs_df = )`. The `obs_df` path runs
  [`obsDF2arrays()`](https://finnverse.github.io/FINN/reference/obsDF2arrays.md),
  which rounds `trees` to whole stems and silently renumbers any
  `siteID`/`patchID` that is not a complete `1:n` sequence; on a
  simulated or partly empty stand that corrupts the anchor, and a
  renumbered site no longer aligns with `init_cohort`.

- anchor_mode:

  (`character(1)`)  
  What an anchor does to the simulated state. `"replace"` (default)
  overwrites it wholesale. `"merge"` keeps the simulated cohorts with
  `dbh < anchor_min_dbh` and concatenates them with the full anchor. Use
  `"merge"` whenever the inventory has a minimum-dbh threshold (FIA:
  12.7 cm): an anchor built from it contains no small trees, so
  `"replace"` deletes FINN's entire regeneration pool at every census –
  on FIA that inflates 200-year equilibrium basal area by 41%. Nothing
  is double counted, because the observation is complete at and above
  the threshold. Ignored unless `reanchor = TRUE`.

- anchor_min_dbh:

  (`numeric(1)`)  
  The inventory's minimum-dbh threshold, in the units of the `dbh`
  column. Required by `anchor_mode = "merge"`, ignored by `"replace"`.

- teacher_forcing:

  (`logical(1)`)  
  Fit by **direct estimation** instead of by simulation, and never run a
  trajectory at all. The cohort state is reset at **every** timestep to
  the observed stand at the census that STARTS that timestep's
  observation interval – `init_cohort` for the first interval, then
  `anchor_cohorts` – so each year's demographic rates are evaluated at
  an observed state and the simulated state is never carried forward.
  The interval rates are then aggregated over the window and scored
  exactly as `fit()` already scores them, so growth, mortality and
  regeneration keep their usual meaning. The three STATE responses
  (`dbh`, `ba`, `trees`) **leave the loss**: predicting them requires a
  trajectory, and there is none here.

  Mutually exclusive with `reanchor`, which is a different estimator:
  re-anchoring still free-runs the whole segment between two censuses
  and scores where it lands. `update_step` is inert under teacher
  forcing – there is no multi-year trajectory to truncate – and the
  anchor validation refuses any `update_step` that would change which
  observations are scored.

  `anchor_cohorts` must supply an anchor for every observed timestep
  except the last (nothing is conditioned on the final census, because
  nothing follows it). Only `fit()` teacher-forces;
  [`predict()`](https://rdrr.io/r/stats/predict.html)/[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  always simulate, so a teacher-forced model is still validated by
  single shooting.

- rb_gradient:

  (`logical(1)`)  
  Take the gradient of the drawn death and recruitment counts through
  their expectation instead of through the draw. The stand still dies
  and recruits discretely, but the reparameterised draw and the
  rounding's surrogate no longer carry the sampler's noise into the
  gradient (on a national inventory most of the mortality gradient was
  draw noise without it). Default `TRUE` since 0.2.0; set `FALSE` to
  reproduce a fit made before that. The setting also decides which death
  draw the default `createProcess(mortality_draw = "binomial")` uses:
  the exact binomial with `TRUE`, the Gamma draw of earlier versions
  with `FALSE` (see
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)).
  It is stored on the model, so later simulations draw the same way as
  the last fit.

- trees_offset:

  (`numeric(1)`)  
  Offset added to the predicted stem count before the `trees` count
  likelihood, which is `-Inf` (poisson) or arbitrarily steep (nbinom)
  where a species' simulated count is exactly zero and the observed one
  is not. On the scale of the predicted count, so per patch. The default
  `0.01` keeps that cell's gradient comparable with every other term;
  the historical `1e-5` (poisson) and `0.001` (nbinom) left it at
  `y/offset`, four orders above the rest. `NULL` restores those two
  historical constants and so reproduces a fit made before 0.2.0.

- ...:

  Arguments of earlier versions, still accepted with a warning:
  `weights` (now `loss_weights`), `recruits` (now `loss_recruit`), and
  `optimizer`, `lr_scheduler`, `lr_scheduler_params`, `clip_norm` (now
  members of `control = finnControl()`). Any other named argument is
  passed on to the optimizer, as before; `finnControl(...)` is the place
  for it now.

- loss:

  Deprecated name of `loss_family`, still accepted with a warning.

## Value

The fitted `model`, invisibly. `fit()` trains the model in place, so the
returned object is the one passed in, now carrying the training results
(e.g. `$history`, `$loss_weights`, `$loss_baseline`, `$param_history`).

## Details

**What happens in a fit.** Every epoch, `fit()` takes the sites in
batches of `batchsize`, simulates each batch over the years of `env`
starting from `init_cohort`, and at every year that has observations
compares the simulated mean dbh, basal area, number of trees, growth,
mortality and recruitment with the observed ones, response by response,
through the likelihoods in `loss_family`. The six comparisons are added
up into one number, the loss, and every parameter is moved a small step
(`lr`) in the direction that lowers it. FINN is written in `torch`, so
these directions (gradients) come from automatic differentiation; no
knowledge of `torch` is needed. After `epochs` passes over the data the
model holds the estimated parameters; `fit()` changes the model in place
and also returns it. The per-response losses of every epoch are in
`model$history`.

**Rates over an observation interval.** An inventory interval of
`period_length` years is that many simulated years. The simulated yearly
rates of those years are compounded into one rate per year before they
meet the observation, with the same functional
[`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
uses for the observed rate: growth as the geometric mean of `1 + g`
minus one, mortality as one minus the geometric mean of the yearly
survival `1 - m`; the recruits are summed and the stand variables are
taken at the census year. A year in which a species has no trees in the
simulation counts as growth 0 and survival 1, and the rate is spread
over the full `period_length`, because the observed rate is derived from
the counts over the whole interval. With `finn(growth_rate = "period")`
growth is the compounded increment over the interval instead.

**Reading the progress.** With the default `loss_weights = "auto"` every
response's loss is reported as a fraction of the loss a constant
prediction (the mean of the observations) would get: `1` means no better
than the mean, below `1` better. `model$history` holds these fractions,
one vector per epoch in the order `dbh`, `ba`, `trees`, `growth`,
`mortality`, `regeneration` (a seventh, unused entry is 0); the plot
drawn with `plot_progress = TRUE` shows each of them relative to its
value in the first epoch. A fit has converged when the parameters stop
moving (`model$param_history`), not merely when the losses look flat;
see the vignette *Fitting FINN to FIA data*.

**Long simulated series.** Taking the gradient through a long simulation
is slow and the signal fades. `fit()` therefore cuts the gradient every
`update_step` years; the gradients of the cuts are added up and the
parameters are updated once per batch of sites. The simulated state
itself runs on uncut. The default of one year is the fastest and the
most short-sighted; larger values let a loss see further back at a
higher cost.

**Advanced estimators.** `reanchor` restarts the simulated stand from
the observed trees at every census (multiple shooting),
`teacher_forcing` evaluates the rates on the observed stands and never
simulates a trajectory, and
`loss_resolution`/`loss_aggregation = "tree"` score single trees. They
are for specific questions; the default single-shooting fit is what the
vignettes use.

## What makes a fit faster, and what does not

A fit is bound by the number of small tensor operations a simulated year
issues, not by the arithmetic in them. So a GPU does not always pay
(measure one configuration before committing a long fit to it), more
threads rarely help (`threads = 1` is usually fastest and the
considerate setting on a shared machine), and several fits run side by
side use a machine better than one fit given all of it. What does
shorten a fit is fewer epochs: watch `model$param_history` and stop when
the parameters settle instead of fixing a large `epochs` in advance.

## Bounds on simulated mortality

Before the death draw, the simulation clamps the per-cohort mortality
probability of every step to the interval
`createProcess(mortality_floor = 0.001, mortality_ceiling = 0.999)` of
the mortality process, so at least 0.1 % and at most 99.9 % of a cohort
is at risk per simulated step however extreme the fitted rate. The bound
is per *step*, not per year: on five-year steps the floor is 0.02 %/yr,
on decadal ones 0.01 %/yr. Until 0.2.0 the interval was hard-coded to
`[0.01, 0.99]`, which made simulated low-mortality stands lose stems the
observations did not. The per-tree mortality loss is not affected by
either bound; it clamps at `1e-6`. See
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md).

## The arguments that shape the loss, and how they combine

Ten arguments change what `fit()` optimises. They are not alternatives
to one another: each sits on its own axis, and a fit picks one value on
every axis.
[finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md)
carries the detail; this is the map.

|  |  |  |  |
|----|----|----|----|
| axis | argument | what it decides | default |
| which responses are scored | `loss_family = c(<response> = "none")` | leaves that response out of the loss entirely | all six scored |
| which likelihood | `loss_family` | the distribution each response is compared with | `mse`, `mse`, `poisson`, `mse`, `binomial`, `nbinom` |
| how the simulation is summarised | `loss_aggregation` | what the cohorts are reduced to before they meet an observation: `"mean"`, `"sum"`, `"size_classes"`, `"quantiles"`, `"tree"` | mean for `dbh`, `growth`, `mortality`; sum for `ba`, `trees`, `regeneration` |
| at what resolution growth and mortality are scored | `loss_resolution` | `"site"`: per site, year and species. `"tree"`: per simulated cohort that started as an observed tree | `"site"` |
| how the six terms are combined | `loss_weights`, `weight_scale` | `"auto"` divides each term by its intercept-only baseline; a numeric vector then multiplies those balanced terms, or the raw ones under `weight_scale = "raw"` | `"auto"`, `"auto"` |
| where the recruits come from | `loss_recruit` | `"observed"` feeds the observed recruitment into the simulation, freezes the regeneration process and drops its loss term | `"simulated"` |
| how the state is propagated | `reanchor` + `anchor_cohorts`, `anchor_mode`, `anchor_min_dbh`; `teacher_forcing` | the estimator – see the table below | single shooting |
| how far the gradient reaches | `update_step` | truncated backpropagation through time: the gradient is cut every `update_step` years, the state is not | `1` |

**The three estimators** differ in how much of the trajectory each
gradient sees, and therefore in what the fit is asked to get right:

|  |  |  |  |
|----|----|----|----|
|  | state over the horizon | scored responses | needs |
| single shooting (default) | initialised once, free-running to the end; error compounds across every census | all six | `init_cohort` |
| multiple shooting (`reanchor = TRUE`) | free-running within a census interval, rebuilt from the observed trees at each anchor after that step's loss | all six | `anchor_cohorts` per observed timestep |
| direct estimation (`teacher_forcing = TRUE`) | never carried forward: reset to the observed stand at every timestep | `growth`, `mortality`, `regeneration` only – `dbh`, `ba` and `trees` leave the loss, because predicting them needs a trajectory | `anchor_cohorts` for every observed timestep but the last |

`update_step` is a different axis from all three: it cuts the gradient
but keeps the state, so under single shooting the stand still compounds
while the parameters only see `update_step` years at a time. It is inert
under teacher forcing, where there is no multi-year trajectory to
truncate.

**What cannot be combined:**

- `reanchor` and `teacher_forcing` – different estimators, not settings
  of one.

- `loss_aggregation = "tree"` and `loss_resolution = "tree"` – the two
  per-tree routes score different objects (the observed tree at its
  measured diameter against the simulated cohort it became) and cannot
  run in one fit;
  [finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md)
  compares them side by side.

**What a value silently implies:** `loss_recruit = "observed"` removes
the regeneration term, so a fit under it cannot be compared with one
under `"simulated"` by its total loss; `teacher_forcing = TRUE` removes
three more; and `loss_weights = "auto"` rescales whatever is left, so
the total loss is only comparable between fits that score the same set
of responses.

## See also

[finn_data](https://finnverse.github.io/FINN/reference/finn_data.md) for
the input tables,
[finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md) for
how the simulation is compared with the observations,
[`finnControl()`](https://finnverse.github.io/FINN/reference/finnControl.md)
for the optimizer settings,
[recruitment](https://finnverse.github.io/FINN/reference/recruitment.md)
for the recruitment options (`loss_recruit`, the regeneration
likelihood),
[`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
to simulate with the fitted model

## Examples

``` r
library(data.table)
# three Oregon FIA sites over two epochs; a real calibration runs hundreds
sites <- 1:3
Nsp   <- max(fia_obs_dt$species)
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
FINN.seed(1)
m <- finn(N_species = Nsp, recruits_dbh = 12.9,
          growth_process       = createProcess(~ temp + prec, growth),
          mortality_process    = createProcess(~ temp + prec, mortality),
          regeneration_process = createProcess(~ temp + prec, regeneration_capped))
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
    plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.768)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%

# the per-response losses over the epochs, and the baselines they were
# scaled by to make `loss_weights = "auto"` comparable across responses
tail(m$history, 2)
#> [[1]]
#> [1]  2.1722367  0.7110548  0.5395446 28.7657299  1.3195336  1.3620560  0.0000000
#> 
#> [[2]]
#> [1]  0.9194432  0.3529012  0.5377366 24.7921906  1.3246045  1.4133602  0.0000000
#> 
m$loss_baseline
#>          dbh           ba        trees       growth    mortality regeneration 
#> 4.117755e+02 8.303134e-02 1.033310e+00 5.429424e-05 2.596614e-02 1.302194e+00 

# give growth twice the weight of every other response, and lower the
# learning rate every 50 epochs
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
    loss_weights = c(1, 1, 1, 2, 1, 1),
    control = finnControl(lr_scheduler = "step",
                          lr_scheduler_params = list(step_size = 50, gamma = 0.5)),
    plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.77)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=14%, ba=14%, trees=14%, growth=29%, mortality=14%, regeneration=14%
```
