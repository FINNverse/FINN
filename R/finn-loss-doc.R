#' @title How fit() compares the simulation with the observations: the loss
#'
#' @description
#' [FINN::fit()] estimates the parameters by making the simulation match the
#' observations. "Match" is measured by the loss: for each of the six
#' responses -- `dbh`, `ba`, `trees`, `growth`, `mortality`, `regeneration` --
#' the simulated value of every site, census year and species is compared
#' with the observed one through a likelihood, the six comparisons are added
#' up, and `fit()` lowers that sum. This page explains the arguments of
#' [FINN::fit()] that set the comparison: `loss_family`, `loss_aggregation`,
#' `loss_resolution`, `loss_weights` and `weight_scale`, `loss_recruit`, and
#' `tree_data`. Their defaults are robust in most cases: a first fit needs
#' none of them.
#'
#' @section The six responses and their likelihoods (`loss_family`):
#' Each response has its own likelihood, chosen for the kind of number it is.
#' The defaults, following Pichler & Käber (2026):
#'
#' | response | what is compared | default | alternatives |
#' |---|---|---|---|
#' | `dbh` | mean diameter, cm | `"mse"` (squared error) | `"gaussian"`, `"none"` |
#' | `ba` | basal area per patch | `"mse"` | `"gaussian"`, `"none"` |
#' | `trees` | number of trees per patch | `"poisson"` (a count) | `"nbinom"`, `"mse"`, `"none"` |
#' | `growth` | mean annual relative growth rate | `"mse"` | `"gaussian"`, `"none"` |
#' | `mortality` | share of the trees that died, 0 to 1 | `"binomial"` | `"mse"`, `"gaussian"`, `"poisson"`, `"none"` |
#' | `regeneration` | recruits per hectare over the interval | `"nbinom"` (an overdispersed count) | `"poisson"`, `"mse"`, `"gaussian"`, `"none"` |
#'
#' `loss_family` is a named character vector; name only what departs from the
#' defaults, e.g. `loss_family = c(trees = "nbinom")`. `"mse"` is the squared
#' error, equivalent to a Gaussian likelihood with a fixed variance;
#' `"gaussian"` additionally estimates the variance. `"binomial"` is the
#' Bernoulli/binomial likelihood for a proportion in \[0, 1\] (the model
#' predicts the probability of dying through a sigmoid), weighted by the
#' number of trees at risk when `data` carries `n_at_risk`; it respects the
#' bounds and the mean-variance link of a proportion, which `"mse"` ignores.
#' `"nbinom"` is implemented for `trees` and `regeneration` only. `"none"`
#' switches a response off: it contributes nothing and needs no column in
#' `data` -- use it when a response was not observed, or when its information
#' enters through another term (e.g. `dbh = "none"` and `ba = "none"` once the
#' size distribution is scored through size classes, below).
#'
#' The rates are compared over the observation interval: with
#' `period_length` in `data`, simulated growth and mortality are averaged
#' over the interval's steps and simulated recruitment is summed, before they
#' meet the observation (see [finn_data]). Growth is compared as a *relative*
#' rate (`dbh / dbh_before - 1`), the model's own quantity; training against
#' the absolute increment was tested and was worse and more seed-variable,
#' because it couples the growth parameters to whichever trees happen to be
#' present. `finn(growth_rate = )` says whether the observed growth is a rate
#' per time step or the increment over the whole interval.
#'
#' For `regeneration`, `"nbinom"` (the default) uses the recruit dispersion
#' `theta` of the regeneration process and applies `finn(recruit_obs_weight
#' = )`; `"poisson"` has no dispersion parameter, so it cannot collapse, but in
#' 0.2.0 it ignores `recruit_obs_weight` and scores per-hectare values
#' (FINNverse/FINN#65), with `trees_offset` added to the prediction. See
#' [recruitment].
#'
#' A missing observation (`NA`) is left out of its term; the other responses
#' of that cell are still scored. A term that is not finite at a step (an
#' overflowing likelihood) is dropped from that step with the default
#' `fit(nonfinite = "drop_term")` and counted in `model$nonfinite_counts`.
#'
#' @section Balancing the six terms (`loss_weights`, `weight_scale`):
#' The six terms are summed, and their raw sizes differ by orders of
#' magnitude: on the bundled FIA data `dbh` is a squared error in cm2 (about
#' 430 at the start) while `growth` is a squared error on a ratio (about
#' 0.03). Summed as they are, `dbh` would take 87 % of the objective and
#' `growth` 1 %, and growth would never improve.
#'
#' `loss_weights = "auto"` (the default) fixes this without knowing the units.
#' Each term is divided by its **intercept-only baseline**: the loss that
#' predicting the single best constant (the mean of the observations) would
#' get. Every term then measures the same thing on the same scale, the
#' fraction of its own null deviance. For a squared error the baseline is the
#' variance of the observations, so this is scaling by `1 / var(y)`; for the
#' Poisson, negative-binomial and binomial terms it is the corresponding null
#' deviance, where a standard deviation would be the wrong scale. The
#' baselines are stored in `model$loss_baseline`, the resulting weights in
#' `model$loss_weights`, and every fit reports each response's share of the
#' objective at the baseline. It is a weighting for a balanced objective, not
#' the joint likelihood of the six responses.
#'
#' A useful side effect: the per-response values in `model$history` read
#' directly as "how much better than the mean" -- `1` is no better than the
#' intercept, below `1` is better, above `1` is worse.
#'
#' A numeric `loss_weights` vector (six numbers, in the order `dbh`, `ba`,
#' `trees`, `growth`, `mortality`, `regeneration`) multiplies the balanced
#' terms, so each number reads as that response's relative share:
#' `c(1, 1, 1, 2, 1, 1)` gives `growth` twice the weight of every other
#' response. On the FIA example, balancing the objective is worth about +0.24
#' Spearman correlation on held-out `growth`, more than tuning `lr` or
#' `epochs` (vignette *Fitting FINN to FIA data*). `weight_scale = FALSE` (or
#' `"raw"`) applies the vector to the raw, unbalanced terms instead -- the
#' meaning the argument `weights` had before 0.2.0, kept for reproducing
#' earlier analyses; the shares are reported in both modes.
#'
#' @section How cohorts become the compared quantity (`loss_aggregation`):
#' The simulation works on cohorts; the observation is one number per site,
#' year and species. `loss_aggregation` says how the cohorts are summarised
#' before the comparison, per response, as a named character vector or -- for
#' a type that needs settings -- a named list. Anything not named keeps its
#' default, so the argument only ever states what departs from it:
#'
#' * `"mean"` -- the abundance-weighted mean over the cohorts of a site and
#'   species. The default for `dbh`, `growth` and `mortality`.
#' * `"sum"` -- the sum over cohorts. The default for `ba`, `trees` and
#'   `regeneration`.
#' * `"size_classes"` -- stem counts per diameter class, for `trees` only.
#'   Needs `breaks` (the class boundaries in cm) and optionally `tau` (the
#'   softness of the class edges in cm, default 1). `K` breaks give `K + 1`
#'   classes, read from the columns `n_class_1`, ..., `n_class_K+1` of `data`.
#'   With a Poisson family this scores total abundance and size composition at
#'   once (independent Poisson class counts factorise into a Poisson total
#'   times a multinomial composition). Class membership is a difference of
#'   sigmoids rather than a hard cut, so the term stays differentiable in the
#'   diameter; small `tau` approaches hard binning.
#' * `"quantiles"` -- the abundance-weighted diameter quantiles instead of
#'   the mean, for `dbh` only. Needs `probs` (e.g. `c(0.1, 0.5, 0.9)`) and
#'   optionally `tau` (the softness on the cumulative-abundance scale,
#'   default 0.05). The observation is read from the columns `dbh_q10`,
#'   `dbh_q50`, `dbh_q90` of `data` (`dbh_q<100*p>`). Where size classes ask
#'   how many stems fall in an interval, this asks how thick the stem at a
#'   given position of the distribution is: the upper quantile speaks about
#'   the big trees that carry the basal area, which a count likelihood weighs
#'   like saplings.
#' * `"tree"` -- no aggregation: the growth and mortality equations are
#'   evaluated on the **observed** trees at their observed diameters, with the
#'   light of the observed stand, and compared tree by tree. For `growth` and
#'   `mortality` only; needs `tree_data`. See *Scoring individual trees*.
#'
#' Example -- score the size distribution through stem counts and the two
#' rates per tree:
#'
#' ```r
#' fit(m, data = obs_dt, env = env, tree_data = tree_obs_dt,
#'     loss_family      = c(dbh = "none", ba = "none"),
#'     loss_aggregation = list(trees  = list(type = "size_classes",
#'                                           breaks = seq(20, 60, by = 10)),
#'                             growth = "tree", mortality = "tree"))
#' ```
#'
#' @section Scoring individual trees (`tree_data`, `loss_resolution`):
#' By default `growth` and `mortality` are scored per site, species and
#' census, as means over cohorts. Two options score them per tree instead;
#' both take the tree table in `tree_data` ([finn_data] lists its columns),
#' but they ask different questions.
#'
#' `loss_aggregation = c(growth = "tree", mortality = "tree")` predicts each
#' OBSERVED tree: the equations run on its measured dbh at the start of the
#' interval with the light of the observed stand, and only the yearly
#' environment changes across the interval. It measures the error of the
#' growth and mortality functions alone. `loss_resolution = "tree"` predicts
#' each SIMULATED cohort that started as an observed tree, along the model's
#' own trajectory (simulated growth, competitors, light). It measures the
#' error of the functions plus the drift of the simulated stand.
#'
#' | | `loss_aggregation = "tree"` | `loss_resolution = "tree"` |
#' |---|---|---|
#' | the tree is | observed: measured dbh, observed competitors, yearly environment | simulated: the cohort it started as, along the model's trajectory |
#' | measures | the growth and mortality functions alone | the functions plus the drift of the simulated stand |
#' | responses | chosen separately, e.g. only `mortality = "tree"` | growth and mortality together |
#' | `init_cohort` | any, incl. binned (`makeInitCohorts(dbh_binsize = )`) | one cohort per tree (`makeInitCohorts()` without `dbh_binsize`) |
#' | `tree_data` needs | `dbh` at the interval start | `cohortID` of the initial cohort |
#' | combines with | the simulated `dbh`/`ba`/`trees`/`regeneration` terms | `reanchor` and `teacher_forcing` |
#' | diagnostics | [predictTrees()] | `$tree_deaths`, `$tree_pred` |
#'
#' Use the first to estimate the process functions from the trees as
#' observed, the second to test whether simulated trajectories follow the
#' real trees. The two cannot be combined in one fit.
#'
#' Under `loss_resolution = "tree"` the two terms **replace** the per-cell
#' terms rather than adding to them (scoring a response twice under two
#' weights is a redundancy the optimiser can trade against itself); the loss
#' vector, `loss_weights` and `$history` keep their six positions. The
#' correspondence is carried by the cohort ids: a cohort keeps the id of the
#' slot it was created in, so cohort *j* of a (site, patch) is tree *j* of the
#' list that generation started from; recruits get ids above that block. A
#' generation starts at `init_cohort` and restarts at every anchor under
#' `reanchor`/`teacher_forcing`. Growth is scored over the interval as the
#' aggregate term is, mortality as the compounded probability
#' `1 - prod(1 - m)` against the tree's 0/1 fate -- never the sampled death,
#' which at one tree is a coin flip. A cohort the simulation kills inside a
#' window is kept with `trees = 0` until its generation ends, so both terms
#' score every tree and the scored sample does not depend on the parameters;
#' `$tree_deaths` reports per census and epoch how often the simulation
#' killed a tree the inventory found alive (`frac_killed`), in which year
#' (`mean_kill_year`), and the share that entered the window already dead
#' (`frac_ghost`). With the binned competition a killed cohort still occupies
#' a slot of the height quantiles, so `createProcess(light_mask = TRUE)` (the
#' default) or `continuous = TRUE` is needed for the competition process;
#' `fit()` warns otherwise.
#'
#' @section Where the recruits come from (`loss_recruit`):
#' `loss_recruit = "simulated"` (the default) lets the regeneration process
#' recruit, and its parameters are estimated with the rest. An inventory
#' samples small trees sparsely, so recruitment is the hardest process to
#' estimate jointly: a mis-fitted regeneration process feeds a wrong stand
#' into every later interval, and growth and mortality compensate for it.
#' `loss_recruit = "observed"` therefore drives the stand with the OBSERVED
#' recruitment: `reg` in `data` divided by `period_length` is the expected
#' recruitment per hectare and year in every simulated year of that
#' interval, drawn per patch by the recruit sampler as usual. In that call
#' the regeneration term is weighted 0, every parameter of the regeneration
#' process is frozen (and gets its previous state back when `fit()`
#' returns), and growth, mortality and competition are estimated on a stand
#' that recruits what was observed. The regeneration process is then fitted
#' in a call of its own, with the other processes' `optimizeSpecies` and
#' `optimizeEnv` off. Refused: an `NA` `reg` in an interval that needs it, a
#' simulated year in two intervals or before a site's first interval, and
#' `reanchor`/`teacher_forcing`. [FINN::simulateForest()] with
#' `recruits = "observed"` simulates the same way for checking such a model;
#' [recruitment] has the unit chain and the staged recipe.
#'
#' @section Reading the fit:
#' `model$history` holds the six (weighted) terms per epoch,
#' `model$loss_baseline` the intercept-only baselines, `model$loss_weights`
#' the weights in force, and `model$param_history` the parameters at every
#' `checkpoints` epochs. A fit has converged when the parameters stop moving,
#' not when the losses look flat.
#'
#' @references
#' Pichler, M., & Käber, Y. (2026). Inferring processes within dynamic forest
#' models using hybrid modelling. \emph{Methods in Ecology and Evolution}.
#' \doi{10.1111/2041-210x.70347}
#'
#' @seealso [FINN::fit()], [finn_data], [recruitment], [FINN::predictTrees()]
#' @name finn_loss
#' @aliases finn-loss
NULL
