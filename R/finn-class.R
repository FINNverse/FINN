#' Forest Informed Neural Network
#'
#' @description
#' Creates a Forest Informed Neural Network (FINN), a differentiable, cohort-based
#' dynamic forest (gap) model in the tradition of JABOWA/ForClim-style models, in
#' which any of the four demographic processes (growth, mortality, regeneration,
#' competition for light) can either be specified mechanistically or replaced by a
#' deep neural network (DNN). Mechanistic and DNN-based processes are calibrated
#' jointly, end-to-end, via gradient descent (see [FINN::fit]).
#'
#' @details
#' FINN represents the forest as cohorts of trees, grouped by site, patch and
#' cohort, each characterized by diameter at breast height (dbh), number of trees,
#' and species identity. Starting from an initial state, FINN simulates the forest
#' forward in discrete annual time steps by sequentially applying the four
#' demographic processes: competition (light availability, based on basal area
#' and species-specific shading), growth (diameter increment as a function of
#' light, size and environment), mortality (binomial death of trees as a function
#' of growth, light, size and environment) and regeneration (recruitment of new
#' cohorts as a function of light and environment, drawn from a negative binomial
#' distribution). Each process additionally depends on a species- and
#' process-specific environmental-response function that maps site-level
#' environmental predictors to a scalar effect on the process.
#'
#' Every process and its environmental-response function can be configured in one
#' of two ways: (1) mechanistically, with an explicit functional form and
#' interpretable parameters (e.g. light-response thresholds, allometric
#' coefficients), created via [FINN::createProcess]; or (2) as a hybrid process, in
#' which the environmental-response function or the entire process equation is
#' replaced by a DNN, created via [FINN::createHybrid]. The remaining mechanistic
#' processes constrain the DNN to ecologically plausible behaviour, while the DNN
#' absorbs misalignments and structural simplifications that would otherwise bias
#' the mechanistic processes; both are estimated jointly rather than calibrated in
#' isolation and plugged in afterwards. If a process argument is left at its
#' default `NULL`, the corresponding default mechanistic process (as described in
#' Pichler & Käber, 2026, Methods in Ecology and Evolution) is used.
#'
#' `finn()` only assembles the model architecture (analogous to instantiating a
#' `torch` `nn_module`); none of the process or environmental-response parameters
#' are estimated yet. Use [FINN::fit] to calibrate the returned object against
#' observed data, and [FINN::predict.finn_class]/[FINN::simulateForest] to
#' simulate forest dynamics from a (fitted) model.
#'
#' @param N_species (`integer(1)`)\cr Number of species.
#' @param mortality_process (`function`)\cr Mortality process, created by
#'   [FINN::createProcess] (mechanistic) or [FINN::createHybrid] (DNN-based). If
#'   `NULL` (default), the default mechanistic mortality process is used.
#' @param growth_process (`function`)\cr Growth process, created by
#'   [FINN::createProcess] (mechanistic) or [FINN::createHybrid] (DNN-based). If
#'   `NULL` (default), the default mechanistic growth process is used.
#' @param regeneration_process (`function`)\cr Regeneration process, created by
#'   [FINN::createProcess] (mechanistic) or [FINN::createHybrid] (DNN-based). If
#'   `NULL` (default), the default mechanistic regeneration process is used.
#' @param competition_process (`function`)\cr Competition (light availability)
#'   process, created by [FINN::createProcess]. If `NULL` (default), the default
#'   mechanistic competition process is used.
#' @param recruits_dbh (`numeric(1)`)\cr Starting dbh for recruits. Has value 1.0 as default.
#'
#' @references
#' Pichler, M., & Käber, Y. (2026). Inferring processes within dynamic forest
#' models using hybrid modelling. \emph{Methods in Ecology and Evolution}.
#' \doi{10.1111/2041-210x.70347}
#'
#' @return An object of class `finn_class` (a [torch::nn_module]): an assembled
#'   but un-fitted FINN model, ready to pass to [fit()] or [simulateForest()].
#' @export
#' @param reg_floor (`numeric(1)`)\cr Additive floor on the expected recruitment of
#'   every species in every patch, in stems ha^-1 per time step, used by
#'   [regeneration()] and [regeneration_saturation()]:
#'   `mean = regP * exp(env) + reg_floor`. It keeps the negative-binomial
#'   likelihood finite where the light gate closes. Default `1e-3`. Earlier versions
#'   hard-coded `0.2`, which with many species adds a large climate-independent
#'   recruitment (65 species x 0.2 = 13 stems ha^-1 yr^-1); models saved before this
#'   argument existed keep that value when reloaded.
#' @param recruit_obs_weight (`numeric(1)`)\cr Observation-operator inclusion
#'   weight applied to both predicted and observed regeneration before the loss
#'   (maps a per-hectare recruitment DENSITY into a sampling design's tallied-count
#'   space, e.g. `g(recruit_dbh)/BAF` for an angle-count inventory). `1.0` (default)
#'   is the per-hectare behaviour.
#' @param growth_period_scale (`logical(1)`)\cr If `TRUE`, growth is scored against
#'   the COMPOUNDED period increment (`prod(1+g)-1`) rather than the mean of the
#'   per-year rates, so a de-annualised period target aligns without collapsing
#'   growth's variance. Default `FALSE` (mean of annual rates).
#' @examplesIf torch::torch_is_installed()
#' # the three demographic processes, each with a linear temperature and
#' # precipitation response; ?createProcess for the hybrid and network variants
#' Nsp <- max(fia_obs_dt$species)
#' m <- finn(
#'   N_species            = Nsp, recruits_dbh = 12.9,
#'   growth_process       = createProcess(~ temp + prec, FINN::growth),
#'   mortality_process    = createProcess(~ temp + prec, FINN::mortality),
#'   regeneration_process = createProcess(~ temp + prec, FINN::regeneration))
#' m
finn = function(N_species,
                mortality_process = NULL,
                growth_process = NULL,
                regeneration_process = NULL,
                competition_process = NULL,
                recruits_dbh = 1.0,
                recruit_obs_weight = 1.0,
                growth_period_scale = FALSE,
                reg_floor = 1e-3) {
  if(!(is.numeric(reg_floor) && length(reg_floor) %in% c(1L, N_species) &&
       all(is.finite(reg_floor)) && all(reg_floor >= 0)))
    stop("`reg_floor` must be non-negative finite numbers, of length 1 or `N_species`.", call. = FALSE)
  finn_class(N_species = N_species,
             mortality_process = mortality_process,
             growth_process = growth_process,
             regeneration_process = regeneration_process,
             competition_process = competition_process,
             recruits_dbh = recruits_dbh,
             recruit_obs_weight = recruit_obs_weight,
             growth_period_scale = growth_period_scale,
             reg_floor = reg_floor)
}


#' Fit FINN
#'
#' @description
#' Fit (calibrate) a FINN model end-to-end by gradient-descent optimization. All
#' mechanistic process parameters, environmental-response coefficients and, if
#' present, the weights of any hybrid (DNN-based) processes are estimated jointly
#' in a single optimization, rather than calibrating each process in isolation.
#'
#' @details
#' Calibration relies on the fact that FINN is implemented in `torch` for R and is
#' therefore fully differentiable: at each simulated time step the predicted stand
#' variables and demographic rates are compared to the observed data through a
#' joint loss function, and gradients of this loss with respect to every model
#' parameter are obtained via automatic differentiation (backpropagation) and used
#' to update the parameters with a `torch` optimizer (Adam, [torch::optim_ignite_adam],
#' by default).
#'
#' The joint loss is the sum of per-variable losses (akin to negative
#' log-likelihoods), one for each of `dbh`, `ba` (basal area), `trees` (number of
#' trees), `growth`, `mortality` and `regeneration`. Following Pichler & Käber
#' (2026), reasonable choices are mean squared error (`"mse"`, equivalent to a
#' Gaussian likelihood) for `dbh` and `ba`, Poisson likelihood for `trees`,
#' negative binomial (`"nbinom"`) for `regeneration`, and `"mse"` for `growth` as
#' a continuous rate. Both rate responses are per YEAR: over a `period_length`
#' window they are compared against the model's mean over that window, so an
#' inventory interval has to be annualised before it becomes a target. `mortality`
#' is an observed annual *proportion* of trees that died
#' and the model predicts it through a sigmoid, so it defaults to `"binomial"` —
#' a Bernoulli/binomial likelihood (binary cross-entropy, which admits fractional
#' targets). This respects the \[0, 1\] support and the mean-variance link of a
#' proportion, both of which `"mse"` ignores (the model also supports
#' `"gaussian"` and `"poisson"` as alternatives via the `loss` argument; see
#' Appendix B of the paper for details). Each loss can be weighted individually via `weights`, and
#' missing values in the observed data are masked out of the corresponding loss
#' term. The model is trained for `epochs` iterations over the (optionally
#' batched and shuffled) data using `optimizer` with learning rate `lr`.
#'
#' Backpropagating gradients through a long simulated time series is prone to
#' vanishing gradients and is computationally expensive. FINN therefore uses
#' truncated backpropagation through time: the computational graph is detached
#' every `update_step` simulated years, gradients are accumulated and the
#' parameters are updated, before the simulation continues. Smaller `update_step`
#' values are faster and avoid vanishing gradients but provide a more short-sighted
#' learning signal; larger values let the loss integrate over a longer trajectory
#' at increased computational cost. `start_time` allows discarding an initial
#' burn-in period of the simulation from the loss, and `checkpoints`/`folder`
#' allow periodically saving the model state during training.
#'
#' @param model (`finn_class`)\cr Object of class `finn_class` created by [FINN::finn].
#' @param data (`data.table|data.frame`)\cr Data about demographic rates and stand
#'   variables must be passed as `data.table` or `data.frame`.
#'
#'   Optional columns that change how the responses are scored:
#'   \itemize{
#'     \item `period_length` — how many simulated timesteps each observation's
#'       inventory interval spans. The rate responses (`growth`, `mortality`) are
#'       then compared against the model's **mean** over that interval, and
#'       `regeneration` against its **sum** (recruits accumulate). Absent or `NA`
#'       means each observation is compared against a single timestep. NOTE: this
#'       must currently be the SAME for every site — a per-site interval is not
#'       supported yet, so inventories with varying remeasurement gaps need
#'       filtering to a constant interval first.
#'     \item `n_at_risk` — trees behind each observed `mortality` proportion; used
#'       to weight the binomial term (see `weights`).
#'     \item `growth_n` — trees behind each observed `growth` mean; used to weight
#'       the squared-error term, since the variance of a mean is \eqn{\sigma^2/n}.
#'   }
#'
#'   On the per-site path (`period_length` differing between sites, or sites
#'   observed in different years) each response's loss of a batch is ONE mean
#'   over all units it scores in that batch — observed trees for
#'   `loss_aggregation = "tree"`, `n_at_risk` / `growth_n` for the weighted
#'   cell terms, cells otherwise — across all observation years: the mean of
#'   each observation year enters with weight \eqn{n_{year} / n_{batch}}. A
#'   batch that scores a single observation year is unaffected.
#' @param env (`data.table|data.frame`)\cr Data with environmental covariates must be passed as `data.table` or `data.frame`.
#' @param disturbance (`data.table|data.frame`)\cr Data with disturbance rates must be passed as `data.table` or `data.frame`.
#' @param patches (`integer(1)`)\cr Number of patches.
#' @param patch_size (`numeric(1)`)\cr Patch size.
#' @param init_cohort (`CohortMat`)\cr Initial cohort matrix of class `CohortMat`, created by [FINN::CohortMat]
#' @param epochs (`integer(1)`)\cr Number of iteration steps. `epochs = 0`
#'   runs the set-up only (data tensors, `loss_baseline`, `loss_weights` for
#'   `weights = "auto"`) and trains nothing; the parameters stay as they were.
#' @param lr (`numeric(1)`)\cr Learning rate of the optimizer.
#' @param lr_scheduler (`character(1)|function`)\cr Learning-rate schedule applied
#'   on top of the constant `lr` above. `"none"` (default) keeps `lr` constant,
#'   reproducing the previous behavior exactly. Built-ins: `"step"` (decay by
#'   `gamma` every `step_size` epochs), `"exponential"` (multiply by `gamma` every
#'   epoch), `"cosine"` (cosine-anneal to `eta_min` over `T_max` epochs),
#'   `"plateau"` (reduce by `factor` after `patience` epochs without improvement
#'   in the total loss). Advanced: pass a `function(optimizer)` returning a
#'   `torch` lr scheduler object (must implement `$step()`).
#' @param lr_scheduler_params (`list()`)\cr Tuning overrides for the scheduler
#'   chosen via `lr_scheduler`, e.g. `list(step_size = 50, gamma = 0.5)` for
#'   `"step"`, `list(T_max = 200)` for `"cosine"`, `list(factor = 0.5, patience =
#'   20)` for `"plateau"`. Unset entries fall back to defaults scaled off
#'   `epochs`/`lr`.
#' @param loss (deprecated)\cr Former name of `loss_family`. Still accepted, with
#'   a warning; it will be removed in a future release.
#' @param loss_family (`character(6)`)\cr Named vector giving the **likelihood**
#'   used for each response. Names should be `dbh`, `ba`, `trees`, `growth`,
#'   `mortality`, and `regeneration`. Supported families are `mse`, `poisson`,
#'   `nbinom`, `gaussian`, `binomial`, and `none`.
#'
#'   `binomial` is a Bernoulli/binomial negative log-likelihood intended for
#'   `mortality`; it expects both the prediction and the observation to be
#'   proportions in \[0, 1\]. `nbinom` is implemented for `trees` and
#'   `regeneration` only. `none` switches a response off entirely — it
#'   contributes a constant zero and no gradient. Use it to retire a response
#'   whose information now enters through another term, for example setting
#'   `dbh = "none"` and `ba = "none"` when the size distribution is fitted
#'   through size-classed stem counts (see `loss_aggregation`).
#' @param loss_aggregation (`character`|`list`)\cr How the model's per-cohort
#'   quantities are reduced before they meet the observation, given per response
#'   as a named character vector or, when a type needs arguments, a named list.
#'   Anything not named keeps the default, so the argument only ever states what
#'   departs from it. Supported types:
#'   \itemize{
#'     \item `"mean"` — one value per site, year and species, the abundance-weighted
#'       mean over cohorts. The default for `dbh`, `growth` and `mortality`.
#'     \item `"sum"` — summed over cohorts. The default for `ba`, `trees` and
#'       `regeneration`.
#'     \item `"size_classes"` — stem counts binned by diameter, requiring
#'       `breaks` (the class boundaries in cm) and optionally `tau` (the
#'       softness of the bin edges in cm, default 1). `K` breaks give `K + 1`
#'       classes, and the observation is read from columns `n_class_1` …
#'       `n_class_K+1` of `data`. Only meaningful for `trees`.
#'
#'       This is the size-distribution counterpart of the mean diameter: with a
#'       `poisson` family it scores total abundance and size composition at
#'       once, since independent Poisson class counts factorise into
#'       `Poisson(N | lambda) x Multinomial(n | N, p)`. Class membership uses a
#'       difference of sigmoids rather than a hard cut so the term stays
#'       differentiable with respect to diameter; `tau` controls that softness
#'       and small values approach hard binning.
#'     \item `"quantiles"` — the abundance-weighted diameter quantiles of each
#'       site and species instead of their mean, requiring `probs` (e.g.
#'       `c(0.1, 0.5, 0.9)`) and optionally `tau` (the softness of the quantile
#'       on the cumulative-abundance scale, default 0.05). The observation is
#'       read from columns `dbh_q10`, `dbh_q50`, `dbh_q90` of `data`
#'       (`<response>_q<100*p>`). Only meaningful for `dbh`.
#'
#'       Where `"size_classes"` asks how many stems fall in a diameter
#'       interval, this asks how thick the stem at a given position of the
#'       distribution is, which keeps the statement on the diameter scale: the
#'       upper quantile speaks directly about the big trees that carry the
#'       basal area, where a count likelihood weighs them the same as saplings.
#'       The quantile is a soft selection along the cumulative abundance, so it
#'       stays differentiable with respect to diameter; `tau` controls that
#'       softness and small values approach the exact weighted quantile.
#'     \item `"tree"` — no aggregation: the growth and mortality kernels are
#'       evaluated on the **observed** trees at their observed diameters, with
#'       light computed from the observed stand, and compared tree by tree. Only
#'       available for `growth` and `mortality`, and requires `tree_data`. This
#'       removes the compositional confound in a stand mean (relative growth
#'       falls several-fold from the smallest to the largest trees, so a mean
#'       moves when the size distribution moves), at the cost of teacher-forcing
#'       those two responses onto observed stands rather than simulated ones.
#'   }
#'
#'   Example — score the size distribution through stem counts and the two rates
#'   per tree:
#'
#'   ```r
#'   fit(m, data = obs_dt, env = env, tree_data = tree_obs_dt,
#'       loss_family      = c(dbh = "none", ba = "none", trees = "poisson",
#'                            growth = "mse", mortality = "binomial",
#'                            regeneration = "nbinom"),
#'       loss_aggregation = list(trees  = list(type = "size_classes",
#'                                            breaks = seq(20, 60, by = 10)),
#'                               growth = "tree", mortality = "tree"))
#'   ```
#' @param tree_data (`data.table|data.frame`)\cr Individual observed trees, one
#'   row per tree and census, required by both ways of scoring single trees
#'   (`loss_aggregation = "tree"` and `loss_resolution = "tree"`, see the section
#'   *Scoring individual trees*). `year` is the census that CLOSES the tree's
#'   interval, i.e. the `year` the same interval carries in `data`. Columns:
#'   \itemize{
#'     \item always: `siteID`, `year`, `patchID`, `species`, `growth` (the
#'       observed relative growth rate per timestep, `NA` if unobserved) and
#'       `died` (0/1 over the interval, `NA` if unobserved);
#'     \item `loss_aggregation = "tree"`: additionally `dbh` (diameter at the
#'       START of the interval, cm), and optionally `trees` (stems the row
#'       stands for, default 1). The observed stand is rebuilt from these rows,
#'       so any `patchID` coding works.
#'     \item `loss_resolution = "tree"`: additionally `cohortID`, the id
#'       [makeInitCohorts()] / [obsDF2arrays()] assigned per (`siteID`,
#'       `patchID`) to the cohort the tree started as — its position along the
#'       cohort dimension of `init_cohort` (`makeInitCohorts(treeID_table =
#'       TRUE)` returns the table that carries it), or of the anchor at the
#'       previous census under `reanchor`/`teacher_forcing`. `patchID` must be
#'       the patch index of that layout. `species` (and `dbh`, if present) are
#'       checked against the cohort, and a row that names a missing or already
#'       used cohort, a site or year that is not scored, or a tree that does not
#'       match its cohort stops the fit: the targets are keyed by position, so a
#'       misaligned table would silently score every tree against another one.
#'   }
#' @details
#' `growth` is compared as a **relative** rate (`dbh/dbh_before - 1`), which is
#' the model's native parameter (`dbh_new = dbh * (1 + g)`). Training against the
#' absolute diameter increment (`dbh * g`) instead was tested over three seeds and
#' was worse on every one — even when scored on the absolute scale — and about
#' four times more seed-variable, because it couples the growth parameter to
#' whichever trees happen to be present.
#'
#' @param weights (`"auto"` or `numeric(6)`)\cr Weights of the six losses, in the
#'   order `dbh`, `ba`, `trees`, `growth`, `mortality`, `regeneration`.
#'
#'   The six terms are summed, and their raw magnitudes differ by orders of
#'   magnitude: on the bundled FIA data `dbh` is an MSE in cm^2 (~430 at the
#'   start of training) while `growth` is an MSE on a ratio (~0.03), a factor of
#'   ~1e4. On the raw scale, equal weights give `dbh` ~87% of the FIA objective
#'   and `growth` ~1%, and `growth` then never improves.
#'
#'   `"auto"` (the default) fixes this without needing to know the units. Each
#'   loss is divided by its **intercept-only baseline** — the loss you would get
#'   by predicting the single best constant (the null model). Every term then
#'   measures the same thing on the same scale: the fraction of its own null
#'   deviance. For a squared-error term the baseline is exactly \eqn{var(y)}, so
#'   this reduces to scaling by \eqn{1/sd(y)^2}; but it generalises to the
#'   Poisson, negative-binomial and binomial terms, where a standard deviation is
#'   not the right scale. It is the same idea as `env_autoscale = TRUE`, applied
#'   to the responses rather than the predictors.
#'
#'   A useful side effect: the per-response values in `m$history` become directly
#'   interpretable as **"how much better than the mean"** — `1` is no better than
#'   the intercept, below `1` is better, above `1` is worse. The baselines are
#'   stored in `m$loss_baseline`, the resulting weights in `m$loss_weights`, and
#'   both are reported once when fitting starts.
#'
#'   A `numeric(6)` multiplies the baseline-scaled losses (since 0.2.0), so each
#'   number reads directly as that response's relative share: `c(1, 1, 1, 2, 1, 1)`
#'   gives `growth` twice the weight of every other response. Every fit reports
#'   each response's resulting share of the objective at the intercept-only
#'   baseline. On the FIA example, balancing the objective is worth ~+0.24
#'   Spearman on held-out `growth`, far more than tuning `lr` or `epochs`. See
#'   the "Fitting FINN to forest inventory data" vignette.
#' @param weight_scale (`character(1)`)\cr How a numeric `weights` vector acts.
#'   `"auto"` (default): it multiplies the baseline-scaled losses, as described
#'   under `weights`. `"raw"`: it multiplies the raw losses and the baseline
#'   scaling is off, which is the meaning `weights` had before 0.2.0 and
#'   reproduces analyses fitted with earlier versions. The shares are reported
#'   in both modes.
#' @param optimizer (`torch_optimizer_generator`)\cr Optimizer from the `torch` package.
#' @param batchsize (`integer(1)`)\cr Batch size, model will be trained in random batch sizes of the data to preserve memory and improve convergence.
#' @section What makes a fit faster, and what does not:
#' A FINN fit is bound by the *number of tensor operations it dispatches*, not by the
#' arithmetic in them: a simulated year issues thousands of small kernels, and a tiny
#' operation costs tens of microseconds of dispatch whatever hardware runs it. Three
#' consequences, each measured on a 32,190 site-period panel with 25 patches:
#'
#' * **A GPU need not help.** The per-tree route ran at 84 s/epoch on an RTX A5000
#'   against 85-100 s/epoch on CPU under the same load, with the card at **3 %
#'   utilisation and 653 MB**. Measure one configuration before committing a long fit
#'   to a device.
#' * **Intra-op threads do not help and eventually cost.** On a 64-core machine,
#'   1/2/4/8/16 threads gave 4.59/4.77/4.68/4.73/**6.06** s/epoch on a comparable
#'   problem; on an 8-core laptop one thread was fastest by 30 %. The knee moves with
#'   the core count, the direction does not. `torch::torch_set_num_threads(1)` is the
#'   considerate setting on a shared machine and costs nothing.
#' * **Concurrency is what pays.** Seven fits side by side cost 13-33 % per epoch each
#'   and roughly doubled total throughput, because the machine was never short of cores.
#'   Run several fits at once rather than giving one fit more of the machine.
#'
#' What does shorten a fit is fewer epochs: read a convergence criterion off
#' `$param_history` and stop when the parameters settle, rather than fixing a large
#' `epochs` in advance.
#'
#' @section Lower bound on simulated mortality:
#' The simulation clamps the per-cohort mortality probability to `[0.01, 0.99]`
#' before the death draw (`mortality_draw`), so **however low the fitted mortality
#' is, at least 1 % of a cohort is at risk per model step**. The bound is per
#' *step*, not per year, so its annual equivalent is roughly `0.01 / dt`: 1 %/yr on
#' annual data, 0.2 %/yr on five-year steps, 0.1 %/yr on decadal ones. It exists
#' because `binomial_from_gamma()`'s shape parameter `n*p/(1-p)` degenerates as
#' `p -> 0`, and it is the reason a simulated low-mortality stand -- a mature
#' unthinned stand, a slow-growing species, an equilibrium run -- loses stems the
#' observations do not. The per-tree mortality loss is **not** affected; it clamps
#' at `1e-6`. So a fit can be correct while the trajectory it simulates is not.
#'
#' @param device (`character(1)`)\cr `"cpu"` or `"gpu"`; `"cuda"` is accepted as a
#'   synonym for `"gpu"`. Only NVIDIA cards are supported. A card does not always
#'   pay: the per-tree route (`loss_resolution = "tree"`, `loss_aggregation = "tree"`)
#'   issues thousands of small kernels per simulated year and is bound by dispatch
#'   rather than arithmetic, so it can run at the same speed on either device with
#'   the card almost idle. Measure one configuration before committing a long fit to
#'   a GPU.
#' @param update_step (`integer(1)`)\cr Number of steps for which the gradient should be calculated. Automatic differentation becomes slow for larger update steps and the risk of vanishing gradients increases.
#' @param start_time (`integer(1)`)\cr Starting from which year should the model be fitted. Can be used to use on burn-in.
#' @param plot_progress (`logical(1)`)\cr Plot fitting progress (losses) or not.
#' @param folder (`character(1)`)\cr Path to folder for saving checkpoint models. If `NULL`, no models will be saved during the training.
#' @param checkpoints (`integer(1)`)\cr Interval size in epochs for saving checkpoint models.
#' @param shuffle (`logical(1)`)\cr Shuffle data or not.
#' @param record_gradients (`logical(1)`)\cr Record the gradients of all parameters or not. Can get large for many epochs.
#' @param env_autoscale How to scale the environmental predictors in `env`.
#'   The learned transformation is stored on the model and re-applied
#'   automatically at every `predict()`/`simulate()` call, so raw (untransformed)
#'   `env` can be passed for both calibration and prediction. Accepts:
#'   * `TRUE` (default) -- z-standardise every predictor (mean / sd);
#'   * `FALSE` -- use `env` exactly as supplied (no scaling);
#'   * a **per-predictor spec** for more nuanced scaling: a length-1 string
#'     applied to all, or a vector/list with one entry per predictor (matched by
#'     name if named, else by position). Each entry is one of `"auto"`
#'     (z-standardise), `"identity"`/`"none"` (leave unchanged), `"0to1"`
#'     (min-max to `[0,1]`), or a function `f(x)` returning `list(center=, scale=)`
#'     for a custom affine scaler. Every mode is affine, so ALE stays invertible.
#'   Example: `env_autoscale = c(mat_c = "auto", management = "identity")` keeps
#'   an ordinal management forcing on its native 0..3 scale while z-scaling
#'   climate. The learned constants are available as `model$env_scaling`.
#' @param clip_norm (`numeric(1)|list()`)\cr Gradient-norm budget passed to
#'   `torch::nn_utils_clip_grad_norm_()`, applied separately to each of three
#'   parameter groups (mechanistic per-species rates, env-effect networks,
#'   loss-distribution nuisance parameters) rather than once globally across all
#'   parameters. A single number (default `2.0`) applies the same budget to
#'   every group; a named list/vector keyed by `"mechanistic"`/`"nn"`/`"loss"`
#'   overrides individual groups, e.g. `clip_norm = list(loss = 5, nn = 1)`.
#' @param nonfinite (`character(1)`)\cr What happens when one of the six loss
#'   terms is not finite at a scored timestep (an overflowing likelihood, a
#'   target outside the family's support). `"drop_term"` (default) leaves only
#'   that term out of the timestep's gradient and reported loss; the other
#'   responses are still backpropagated there, and the occurrences are counted
#'   per response in `model$nonfinite_counts` and summarised in a warning when
#'   the fit ends. `"skip_step"` discards the WHOLE timestep — none of the six
#'   terms is backpropagated or reported there — which is what every version
#'   before `drop_term` existed did; use it to reproduce such a fit. The
#'   difference only shows once a term overflows: with finite losses the two
#'   are identical.
#' @param reanchor (`logical(1)`)\cr Fit by **multiple shooting** instead of single
#'   shooting. With the default `FALSE` the cohort state is initialised once from
#'   `init_cohort` and then runs free over the whole horizon, so simulation error
#'   compounds across every observation. With `TRUE` the simulated state is
#'   discarded at every timestep for which `anchor_cohorts` supplies an observed
#'   tree list and rebuilt from it, after that timestep's loss has been
#'   backpropagated — drift is then bounded by one inter-census segment. This is
#'   not the same as the gradient truncation controlled by `update_step`, which
#'   cuts the gradient but keeps the state; set `update_step` to the segment length
#'   to let the gradient span the whole segment.
#' @param anchor_cohorts (`list()`)\cr Named list of [FINN::CohortMat] objects, the
#'   observed state to restart from. Names are the simulation timesteps to
#'   re-anchor at (integers as characters, and entries of the model's internal
#'   `year_sequence`, i.e. timesteps that carry an observation). Every anchor must
#'   have the same number of species, sites and patches as `init_cohort`; the
#'   number of cohorts may differ between anchors and from `init_cohort`. The last
#'   observation year needs no anchor — nothing is simulated after it — and an
#'   anchor may not fall inside another observation's `period_length` window,
#'   because the aggregated rates would then mix pre- and post-anchor state.
#'   `fit()` copies the anchors onto its own device and never modifies them, so one
#'   list can be reused across fits. Ignored unless `reanchor = TRUE`.
#'
#'   **Build the arrays directly** —
#'   `CohortMat(dbh = , trees = , species = , dims = , sp = )` — rather than via
#'   `CohortMat(obs_df = )`. The `obs_df` path runs `obsDF2arrays()`, which rounds
#'   `trees` to whole stems and silently renumbers any `siteID`/`patchID` that is
#'   not a complete `1:n` sequence; on a simulated or partly empty stand that
#'   corrupts the anchor, and a renumbered site no longer aligns with `init_cohort`.
#' @param anchor_mode (`character(1)`)\cr What an anchor does to the simulated state.
#'   `"replace"` (default) overwrites it wholesale. `"merge"` keeps the simulated
#'   cohorts with `dbh < anchor_min_dbh` and concatenates them with the full anchor.
#'   Use `"merge"` whenever the inventory has a minimum-dbh threshold (FIA: 12.7 cm):
#'   an anchor built from it contains no small trees, so `"replace"` deletes FINN's
#'   entire regeneration pool at every census — on FIA that inflates 200-year
#'   equilibrium basal area by 41%. Nothing is double counted, because the
#'   observation is complete at and above the threshold. Ignored unless
#'   `reanchor = TRUE`.
#' @param anchor_min_dbh (`numeric(1)`)\cr The inventory's minimum-dbh threshold, in
#'   the units of the `dbh` column. Required by `anchor_mode = "merge"`, ignored by
#'   `"replace"`.
#' @param teacher_forcing (`logical(1)`)\cr Fit by **direct estimation** instead of
#'   by simulation, and never run a trajectory at all. The cohort state is reset at
#'   **every** timestep to the observed stand at the census that STARTS that
#'   timestep's observation interval — `init_cohort` for the first interval, then
#'   `anchor_cohorts` — so each year's demographic rates are evaluated at an
#'   observed state and the simulated state is never carried forward. The interval
#'   rates are then aggregated over the window and scored exactly as `fit()`
#'   already scores them, so growth, mortality and regeneration keep their usual
#'   meaning. The three STATE responses (`dbh`, `ba`, `trees`) **leave the loss**:
#'   predicting them requires a trajectory, and there is none here.
#'
#'   Mutually exclusive with `reanchor`, which is a different estimator:
#'   re-anchoring still free-runs the whole segment between two censuses and scores
#'   where it lands. `update_step` is inert under teacher forcing — there is no
#'   multi-year trajectory to truncate — and the anchor validation refuses any
#'   `update_step` that would change which observations are scored.
#'
#'   `anchor_cohorts` must supply an anchor for every observed timestep except the
#'   last (nothing is conditioned on the final census, because nothing follows it).
#'   Only `fit()` teacher-forces; `predict()`/`simulateForest()` always simulate,
#'   so a teacher-forced model is still validated by single shooting.
#' @param loss_resolution (`character(1)`)\cr The resolution the `growth` and
#'   `mortality` terms are scored at. `"site"` (the default) is the validated path
#'   and is unchanged in every respect: both are per (site, census, species) cell.
#'   `"tree"` scores the same two responses per individual SIMULATED tree and
#'   **replaces** the cell-resolution terms rather than adding to them — scoring
#'   one response twice under two independent weights is a redundancy the
#'   optimiser can trade against itself. Needs `tree_data` with a `cohortID`
#'   column and an `init_cohort` with one cohort per tree (see the section
#'   *Scoring individual trees* for how it differs from `loss_aggregation =
#'   "tree"`; the two cannot be combined). The loss vector, `weights` and
#'   `$history` keep their six positions and their meaning in both modes; only
#'   the resolution of positions 4 and 5 changes, and the likelihood family is
#'   the one `loss_family` names for them. `dbh`, `ba`, `trees` and
#'   `regeneration` are untouched — a recruit has no per-tree counterpart at all.
#'
#'   The correspondence is carried by the cohort ids: a cohort keeps the id of
#'   the slot it was created in, so cohort *j* of a (site, patch) in a generation
#'   is tree *j* of the list that generation started from. Recruits get ids above
#'   that block and can never enter a residual. A generation starts at
#'   `init_cohort` and restarts at every anchor under `reanchor`/`teacher_forcing`;
#'   `anchor_mode = "merge"` keeps the carried sub-threshold cohorts out of it,
#'   which the anchor never saw. Growth is scored over the interval exactly as the
#'   aggregate term is, mortality as the compounded probability `1 - prod(1 - m)`
#'   against the tree's 0/1 fate — never the sampled death, which at `trees = 1`
#'   is a coin flip.
#'
#'   The term does not depend on the fitting regime but its residual does: under
#'   `teacher_forcing` the kernel runs on the observed tree, so the residual is
#'   parameter error alone; under `reanchor` the cohort has drifted for at most one
#'   census interval; under single shooting for up to the whole horizon.
#'
#'   A cohort the simulation kills inside a window is **kept** with `trees = 0` for
#'   as long as its generation lasts — to the next anchor under
#'   `reanchor`/`teacher_forcing`, which restart the generation there, and
#'   to the end of the horizon under single shooting, where `tree_live0` is rebuilt
#'   only at anchors. It shades nothing and contributes nothing to `dbh`, `ba` or
#'   `trees`, but its diameter and its hazard keep accumulating, because the
#'   observation it is scored against is the mean annual rate and the 0/1 fate over
#'   the whole interval. Both terms therefore score every tree, and the scored
#'   sample does not depend on the parameters. `$tree_deaths` reports per census and
#'   epoch how often the simulation killed a tree the inventory found alive
#'   (`frac_killed`) and in which year (`mean_kill_year`) — a diagnostic, not a
#'   correction. A generation outlives a window under single shooting, so a slot
#'   killed in an earlier window enters the next one already dead; `frac_ghost`
#'   is that share and it is excluded from `frac_killed`, which counts this
#'   window's kills alone. Under `reanchor`/`teacher_forcing` the generation
#'   restarts at every anchor, so `frac_ghost` is 0 by construction only where the
#'   anchors cover every census; a sparser anchor list carries ghosts like single
#'   shooting does.
#' @param rb_gradient (`logical(1)`)\cr Take the gradient of the drawn death and
#'   recruitment counts through their expectation instead of through the draw. The
#'   forward pass is unchanged to the bit -- the stand still dies and recruits
#'   discretely -- but the reparameterised draw and the rounding's surrogate no
#'   longer carry the sampler's noise into the gradient. On the German inventory
#'   97 % of the mortality gradient was draw noise without it. Default `TRUE` since
#'   0.2.0; set `FALSE` to reproduce a fit made before that. See FINNverse/FINN#53.
#' @param trees_offset (`numeric(1)`)\cr Offset added to the predicted stem count
#'   before the `trees` count likelihood, which is `-Inf` (poisson) or arbitrarily
#'   steep (nbinom) where a species' simulated count is exactly zero and the
#'   observed one is not. On the scale of the predicted count, so per patch. The
#'   default `0.01` keeps that cell's gradient comparable with every other term;
#'   the historical `1e-5` (poisson) and `0.001` (nbinom) left it at `y/offset`,
#'   four orders above the rest. `NULL` restores those two historical constants
#'   and so reproduces a fit made before 0.2.0.
#' @param ... Additional arguments passed to `optimizer`.
#'
#' @section Scoring individual trees:
#' By default, `growth` and `mortality` are scored per site, species and census,
#' as means over cohorts. Two options score them per observed tree instead; both
#' take the tree table in `tree_data`, but they ask different questions.
#'
#' `loss_aggregation = c(growth = "tree", mortality = "tree")` predicts each
#' OBSERVED tree: the kernels run on its measured dbh at the interval start with
#' light from the observed stand, and only the yearly environment changes across
#' the interval. It measures the error of the growth and mortality functions
#' alone. `loss_resolution = "tree"` predicts each SIMULATED cohort that started
#' as an observed tree, along the model's own trajectory (simulated growth,
#' competitors, light). It measures the error of the functions plus the drift of
#' the simulated stand.
#'
#' | | `loss_aggregation = "tree"` | `loss_resolution = "tree"` |
#' |---|---|---|
#' | the tree is | observed: measured dbh, observed competitors, yearly environment | simulated: the cohort it started as, along the model's trajectory |
#' | measures | the growth and mortality functions alone | the functions plus the drift of the simulated stand |
#' | responses | chosen separately, e.g. only `mortality = "tree"` | growth and mortality together |
#' | `init_cohort` | any, incl. binned (`makeInitCohorts(dbh_binsize = )`) | one cohort per tree (`makeInitCohorts()` without `dbh_binsize`) |
#' | `tree_data` needs | `dbh` at the interval start | `cohortID` of the initial cohort |
#' | combines with | the simulated `dbh`/`ba`/`trees`/`regeneration` terms | `reanchor` and `teacher_forcing` (under teacher forcing the trajectory is reset to the observed stand every year, which comes close to the first option) |
#' | diagnostics | [predictTrees()] | `$tree_deaths`, `$tree_pred` |
#'
#' Use the first to estimate the process functions from the trees as observed,
#' the second to test whether simulated trajectories follow the real trees. The
#' two cannot be combined in one fit.
#'
#' @return The fitted `model`, invisibly. `fit()` trains the model in place, so
#'   the returned object is the one passed in, now carrying the training results
#'   (e.g. `$history`, `$loss_weights`, `$loss_baseline`).
#' @export
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # three Oregon FIA sites over two epochs; a real calibration runs hundreds
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(
#'   N_species            = Nsp, recruits_dbh = 12.9,
#'   growth_process       = createProcess(~ temp + prec, FINN::growth),
#'   mortality_process    = createProcess(~ temp + prec, FINN::mortality),
#'   regeneration_process = createProcess(~ temp + prec, FINN::regeneration))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
#'     plot_progress = FALSE)
#'
#' # the per-response losses over the epochs, and the baselines they were
#' # scaled by to make `weights = "auto"` comparable across responses
#' tail(m$history, 2)
#' m$loss_baseline
fit = function(model,
               data = NULL,
               env,
               disturbance = NULL,
               patches = 100L,
               patch_size = 0.1,
               init_cohort = NULL,
               epochs = 20L,
               lr = 0.01,
               lr_scheduler = "none",
               lr_scheduler_params = list(),
               loss = NULL,
               loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
               loss_aggregation = NULL,
               tree_data = NULL,
               weights = "auto",
               weight_scale = c("auto", "raw"),
               optimizer = optim_ignite_adam,
               batchsize = NULL,
               device = c("cpu", "gpu"),
               update_step = 1L,
               start_time = 1L,
               plot_progress = TRUE,
               folder = NULL,
               checkpoints = 100L,
               shuffle = TRUE,
               record_gradients = FALSE,
               env_autoscale = TRUE,
               clip_norm = 2.0,
               nonfinite = c("drop_term", "skip_step"),
               reanchor = FALSE,
               anchor_cohorts = NULL,
               anchor_mode = "replace",
               anchor_min_dbh = NULL,
               teacher_forcing = FALSE,
               loss_resolution = c("site", "tree"),
               rb_gradient = TRUE,
               trees_offset = 0.01,
               ...) {
  invisible(model$fit(data = data,
                      env = env,
                      disturbance = disturbance,
                      patches = patches,
                      patch_size = patch_size,
                      init_cohort = init_cohort,
                      epochs = epochs,
                      lr = lr,
                      lr_scheduler = lr_scheduler,
                      lr_scheduler_params = lr_scheduler_params,
                      loss = loss,
                      loss_family = loss_family,
                      loss_aggregation = loss_aggregation,
                      tree_data = tree_data,
                      weights = weights,
                      weight_scale = weight_scale,
                      optimizer = optimizer,
                      batchsize = batchsize,
                      device = device,
                      update_step = update_step,
                      start_time = start_time,
                      plot_progress = plot_progress,
                      folder = folder,
                      checkpoints = checkpoints,
                      shuffle = shuffle,
                      record_gradients = record_gradients,
                      env_autoscale = env_autoscale,
                      clip_norm = clip_norm,
                      nonfinite = match.arg(nonfinite),
                      reanchor = reanchor,
                      anchor_cohorts = anchor_cohorts,
                      anchor_mode = anchor_mode,
                      anchor_min_dbh = anchor_min_dbh,
                      teacher_forcing = teacher_forcing,
                      rb_gradient = rb_gradient,
                      trees_offset = trees_offset,
                      loss_resolution = match.arg(loss_resolution),
                      ...))
}


#' Predict growth and mortality for individual observed trees
#'
#' @details
#' The tree-level counterpart of [FINN::predict.finn_class]. Where `predict()`
#' simulates a stand forward and reports site-level states, this evaluates the
#' growth and mortality kernels on trees you observed: light is computed from
#' the observed stand, the kernels run on those trees at their observed
#' diameters, and the annual rates are accumulated over each site's
#' remeasurement interval — growth as the mean of the annual rates, mortality as
#' the complement of the product of annual survivals.
#'
#' It is the same code path `fit(..., loss_aggregation = c(growth = "tree",
#' mortality = "tree"))` scores, so a validation report cannot drift from the
#' objective it is validating.
#'
#' @param model (`finn_class`)\cr A (fitted) model.
#' @param tree_data (`data.table|data.frame`)\cr Observed trees: `siteID`,
#'   `year`, `patchID`, `species`, `dbh` (at the interval start) and
#'   `period_length` (the interval, in timesteps).
#' @param env (`data.table|data.frame`)\cr Environment, as for [FINN::fit].
#' @param patches (`integer(1)`)\cr Not used: the patch dimension of the
#'   observed-tree tensors is the number of distinct `patchID`s in `tree_data`
#'   (it used to be `patches`, and the extra patches were empty). Kept so
#'   existing calls run unchanged.
#' @param patch_size (`numeric(1)`)\cr Patch size. Defaults to the one the model
#'   was fitted with (`model$patch_size_ha`); pass a value only to override it.
#' @param device (`character(1)`)\cr `"cpu"` or `"gpu"`; `"cuda"` is accepted as a
#'   synonym for `"gpu"`.
#' @return One row per row of `tree_data`, in the same order (`.row` records the
#'   original position), with `growth_pred` (mean annual relative growth over
#'   the interval) and `mort_pred` (probability of dying during the interval).
#'   Bind these onto `tree_data` directly rather than joining on
#'   siteID/year/patchID, which would be many-to-many.
#' @export
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # the growth and mortality kernels applied to the trees as observed, so the
#' # error of the two functions alone, with no simulated stand in between
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(
#'   N_species            = Nsp, recruits_dbh = 12.9,
#'   growth_process       = createProcess(~ temp + prec, FINN::growth),
#'   mortality_process    = createProcess(~ temp + prec, FINN::mortality),
#'   regeneration_process = createProcess(~ temp + prec, FINN::regeneration))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
#'     plot_progress = FALSE)
#'
#' # `fia_init_trees` is a pure initial state, so name the census the trees are
#' # scored against and the interval, in timesteps, that reaches it
#' trees <- fia_init_trees[siteID %in% sites,
#'                         .(siteID, patchID, species, dbh, year = 1L,
#'                           period_length = 1)]
#' pred  <- predictTrees(m, tree_data = trees,
#'                       env = fia_env_dt[siteID %in% sites], patch_size = 0.06)
#' head(pred)
predictTrees = function(model, tree_data, env, patches = 100L, patch_size = NULL, device = c("cpu", "gpu")) {
  model$predict_trees(tree_data = tree_data, env = env, patches = patches,
                      patch_size = patch_size, device = resolve_device(device))
}


#' Predict from a FINN model
#'
#' @details
#' Simulate from a (fitted) FINN model. This is an S3 method for the
#' [stats::predict] generic, so it is dispatched as `predict(model, ...)`.
#'
#' @param object (`finn_class`)\cr Object of class `finn_class` created by [FINN::finn].
#' @param env (`data.table|data.frame`)\cr Data with environmental covariates must be passed as `data.table` or `data.frame`.
#' @param disturbance (`data.table|data.frame`)\cr Data with disturbance rates must be passed as `data.table` or `data.frame`.
#' @param patches (`integer(1)`)\cr Number of patches.
#' @param patch_size (`numeric(1)`)\cr Patch size.
#' @param init_cohort (`CohortMat`)\cr Initial cohort matrix of class `CohortMat`, created by [FINN::CohortMat].
#' @param device (`character(1)`)\cr `"cpu"` or `"gpu"`; `"cuda"` is accepted as a
#'   synonym for `"gpu"`. Only NVIDIA cards are supported.
#' @param return_cohorts Controls whether the raw per-cohort state is returned in
#'   addition to the aggregated site output. Storing cohorts every timestep is
#'   expensive, so the default is `FALSE` (none). Use `TRUE` (or `"all"`) for
#'   every timestep, `"last"` for the final timestep only, or an integer vector
#'   of timesteps to store just those. Recorded cohorts appear as `$long$cohort`
#'   / `$wide$cohort`.
#' @param debug (`logical(1)`)\cr Debug modus or not. If `TRUE`, individual tree states are stored.
#' @param ... Advanced options forwarded to the internal simulator, chiefly
#'   `batchsize` (split the sites into batches of this many to cap memory for very
#'   large runs). The default processes all sites in one batch and is what almost
#'   all users want.
#'
#' @return A named list of predictions. `$long$site` (and `$wide$site`) give the
#'   site-level results (columns `siteID`, `year`, `species`, `variable`,
#'   `value`). When `return_cohorts` is set, `$long$cohort` / `$wide$cohort` add
#'   the per-cohort state (`dbh`, `trees`, `species`, growth `g`, mortality `m`,
#'   ...) for the requested timesteps.
#' @method predict finn_class
#' @export
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # simulate three Oregon FIA sites over the years in `env`, starting from the
#' # observed stands. The model is unfitted here, so these are its priors.
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' m <- finn(
#'   N_species            = Nsp, recruits_dbh = 12.9,
#'   growth_process       = createProcess(~ temp + prec, FINN::growth),
#'   mortality_process    = createProcess(~ temp + prec, FINN::mortality),
#'   regeneration_process = createProcess(~ temp + prec, FINN::regeneration))
#' sim <- predict(m, env = fia_env_dt[siteID %in% sites], init_cohort = init,
#'                patches = 4, patch_size = 0.06, return_cohorts = "last")
#' head(sim$long$site)
#'
#' # the raw per-cohort state of the final timestep, kept because of
#' # `return_cohorts`; storing every timestep is expensive
#' head(sim$long$cohort)
predict.finn_class = function(object,
                              env,
                              disturbance = NULL,
                              patches = 100L,
                              patch_size = 0.1,
                              init_cohort = NULL,
                              device = c("cpu", "gpu"),
                              return_cohorts = FALSE,
                              debug = FALSE,
                              ...) {
  object$simulate(env = env,
                  disturbance = disturbance,
                  patches = patches,
                  patch_size = patch_size,
                  init_cohort = init_cohort,
                  device = device,
                  return_cohorts = return_cohorts,
                  debug = debug,
                  ...)
}


#' Simulate
#'
#' @details
#' Simulate from a fitted FINN model. This is a thin, backwards-compatible
#' alias for [FINN::predict.finn_class]; new code should call
#' `predict(model, ...)` directly.
#'
#' @param model (`finn_class`)\cr Object of class `finn_class` created by [FINN::finn].
#' @param env (`data.table|data.frame`)\cr Data with environmental covariates must be passed as `data.table` or `data.frame`.
#' @param disturbance (`data.table|data.frame`)\cr Data with disturbance rates must be passed as `data.table` or `data.frame`.
#' @param patches (`integer(1)`)\cr Number of patches.
#' @param patch_size (`numeric(1)`)\cr Patch size.
#' @param init_cohort (`CohortMat`)\cr Initial cohort matrix of class `CohortMat`, created by [FINN::CohortMat].
#' @param device (`character(1)`)\cr `"cpu"` or `"gpu"`; `"cuda"` is accepted as a
#'   synonym for `"gpu"`. Only NVIDIA cards are supported.
#' @param return_cohorts Controls whether the raw per-cohort state is returned in
#'   addition to the aggregated site output. Storing cohorts every timestep is
#'   expensive, so the default is `FALSE` (none). Use `TRUE` (or `"all"`) for
#'   every timestep, `"last"` for the final timestep only, or an integer vector
#'   of timesteps to store just those. Recorded cohorts appear as `$long$cohort`
#'   / `$wide$cohort`.
#' @param debug (`logical(1)`)\cr Debug modus or not. If `TRUE`, individual tree states are stored.
#' @param ... Advanced options forwarded to the internal simulator, chiefly
#'   `batchsize` (split the sites into batches of this many to cap memory for very
#'   large runs). The default processes all sites in one batch and is what almost
#'   all users want.
#'
#' @return A named list of simulation results. The patch-averaged, stand-level
#'   state variables and demographic rates are in `$long$site` / `$wide$site`
#'   (long format: `siteID`, `year`, `species`, `variable`, `value`). When
#'   `return_cohorts` is set, `$long$cohort` / `$wide$cohort` hold the raw
#'   per-cohort state for the requested timesteps.
#' @export
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' m <- finn(
#'   N_species            = Nsp, recruits_dbh = 12.9,
#'   growth_process       = createProcess(~ temp + prec, FINN::growth),
#'   mortality_process    = createProcess(~ temp + prec, FINN::mortality),
#'   regeneration_process = createProcess(~ temp + prec, FINN::regeneration))
#' # identical to predict(m, ...), which is what new code should call
#' sim <- simulateForest(m, env = fia_env_dt[siteID %in% sites],
#'                       init_cohort = init, patches = 4, patch_size = 0.06)
#' head(sim$long$site)
simulateForest = function(model,
                          env,
                          disturbance = NULL,
                          patches = 100L,
                          patch_size = 0.1,
                          init_cohort = NULL,
                          device = c("cpu", "gpu"),
                          return_cohorts = FALSE,
                          debug = FALSE,
                          ...) {
  predict(model,
          env = env,
          disturbance = disturbance,
          patches = patches,
          patch_size = patch_size,
          init_cohort = init_cohort,
          device = device,
          return_cohorts = return_cohorts,
          debug = debug,
          ...)
}

# Translate the user-facing `return_cohorts` argument into the concrete vector of
# timesteps whose cohort states should be stored. Keeping this separate makes the
# rule (FALSE = none, TRUE / "all" = every step, "last" = final step, or an
# explicit integer vector) testable in isolation. `debug = TRUE` keeps its legacy
# meaning and records every timestep regardless of `return_cohorts`.
resolve_record_years = function(return_cohorts, time, debug = FALSE) {
  time = as.integer(time)
  if (isTRUE(debug)) return(seq_len(time))
  if (is.null(return_cohorts) || isFALSE(return_cohorts)) return(integer(0))
  if (isTRUE(return_cohorts)) return(seq_len(time))
  if (is.character(return_cohorts)) {
    if (length(return_cohorts) == 1L && return_cohorts %in% c("last", "all")) {
      return(if (return_cohorts == "last") time else seq_len(time))
    }
    stop('`return_cohorts` must be FALSE, TRUE, "last", or a vector of timesteps.',
         call. = FALSE)
  }
  if (is.numeric(return_cohorts)) {
    yrs = as.integer(round(return_cohorts))
    if (anyNA(yrs) || any(yrs < 1L) || any(yrs > time)) {
      stop(sprintf("`return_cohorts` timesteps must all be between 1 and %d.", time),
           call. = FALSE)
    }
    return(sort(unique(yrs)))
  }
  stop('`return_cohorts` must be FALSE, TRUE, "last", or a vector of timesteps.',
       call. = FALSE)
}

# The site sets of the tables fit() and predict() pair. Every array is indexed by
# the POSITION of a site among its table's sorted unique siteIDs, not by the id,
# so the ids need not be 1:N -- but `data`, `env` and `init_cohort` must hold the
# same sites. Otherwise a site is silently paired with another site's
# environment, or torch stops with an index error that names nothing (#33).
# `init_cohort` built from a table carries its user ids ($site_ids); one built
# from arrays is sites 1..N and is checked by count. Returns the sorted user ids:
# the lookup from array position to siteID.
check_site_sets = function(env, data = NULL, init_cohort = NULL, caller = "fit()") {
  ids_of = function(tab, what) {
    id = tab[["siteID"]]
    if (is.null(id)) stop(sprintf("%s: `%s` has no `siteID` column.", caller, what), call. = FALSE)
    if (!is.numeric(id) || anyNA(id) || any(id != round(id)))
      stop(sprintf("%s: `siteID` in `%s` must be integer-valued and not missing; resolveSiteIDs() turns site names into integer ids.", caller, what), call. = FALSE)
    sort(unique(id))
  }
  show = function(x) {
    x = sort(x)
    paste0(paste(utils::head(x, 10L), collapse = ", "),
           if (length(x) > 10L) sprintf(", ... (%d in all)", length(x)) else "")
  }
  env_ids = ids_of(env, "env")
  msg = character(0)
  if (!is.null(data)) {
    data_ids = ids_of(data, "data")
    a = setdiff(data_ids, env_ids); b = setdiff(env_ids, data_ids)
    if (length(a)) msg = c(msg, sprintf("in `data` but not in `env`: %s", show(a)))
    if (length(b)) msg = c(msg, sprintf("in `env` but not in `data`: %s", show(b)))
  }
  if (!is.null(init_cohort)) {
    init_ids = init_cohort$site_ids
    if (!is.null(init_ids)) {
      a = setdiff(init_ids, env_ids); b = setdiff(env_ids, init_ids)
      if (length(a)) msg = c(msg, sprintf("in `init_cohort` but not in `env`: %s", show(a)))
      if (length(b)) msg = c(msg, sprintf("in `env` but not in `init_cohort`: %s", show(b)))
    } else {
      n_init = init_cohort$species$shape[1]
      if (n_init != length(env_ids))
        msg = c(msg, sprintf("`init_cohort` (built from arrays, sites 1..%d) has %d sites, `env` has %d",
                             n_init, n_init, length(env_ids)))
    }
  }
  if (length(msg))
    stop(sprintf("%s: `data`, `env` and `init_cohort` must hold the same sites (the siteIDs need not be 1:N). Sites %s.",
                 caller, paste(msg, collapse = "; ")), call. = FALSE)
  env_ids
}

# Model fields FINN sets itself, during construction or in fit() (#42). A custom
# parameter of the same name would be registered first, and the later assignment
# would then be a silent no-op, because an nn_module drops an assignment to a
# name already registered as a parameter (e.g. `light_mask`, which the
# competition process sets after mortality registered its custom parameters).
# Checked against this fixed list, not against the fields that happen to exist
# when a process is registered. Prefixes and suffixes cover the per-process
# fields (par_<type>, nn_<type>, <type>_func, ...).
reserved_model_fields = c(
  "N_species", "ale", "anchor_cohorts", "anchor_min_dbh", "anchor_mode", "class_breaks",
  "class_obs", "class_tau", "conditional_effects", "continuous", "device", "dgtl_embedder",
  "dtype", "emb_dim", "env_embedder", "env_scaling", "fit_id", "gradients", "growth_period_scale",
  "history", "init_cohort", "light_mask", "loss_agg", "loss_baseline", "loss_raw",
  "loss_resolution", "loss_weights", "mortality_ceiling", "mortality_draw", "mortality_floor", "n_quantiles", "nonfinite",
  "nonfinite_counts", "num_species", "obs_rec", "optimizer", "optimizer_lr", "param_history",
  "patch_size_ha", "per_site_period", "period_mode", "perm_importance", "positional_encoding",
  "pred_rec", "quant_obs", "quant_probs", "quant_tau", "reanchor", "record_raws",
  "recruit_obs_weight", "recruits_dbh", "reg_floor", "regeneration_sampler",
  "rb_gradient", "trees_offset", "reg_cap", "reg_prior", "reg_prior_lambda", "reg_prior_value", "sample_regeneration", "site_lookup", "species_embedder", "teacher_forcing", "train_env",
  "train_init_cohort", "transformer_encoder", "tree_data", "tree_deaths", "tree_obs", "tree_pred")
is_reserved_model_field = function(name) {
  name %in% reserved_model_fields ||
    grepl("^(par|env|nn|loss|process)_|_(formula|func|optimized)$", name)
}

#' finn class (internal representation of FINN); use [finn()] to construct.
#' @noRd
finn_class = nn_module(
  "finn_class",
  initialize = function(
    N_species,
    mortality_process = NULL,
    growth_process = NULL,
    regeneration_process = NULL,
    competition_process = NULL,
    recruits_dbh = 1.0,
    recruit_obs_weight = 1.0,
    growth_period_scale = FALSE,
    reg_floor = 1e-3
  ) {
    self$N_species = N_species
    self$recruits_dbh = recruits_dbh
    self$recruit_obs_weight = recruit_obs_weight
    ## additive recruitment floor (stems ha^-1 step^-1 per species) used by the
    ## mechanistic regeneration functions; see finn() docs
    self$reg_floor = reg_floor
    self$growth_period_scale = growth_period_scale
    ## How each response is reduced from cohorts to the observed quantity. Filled
    ## by fit() from `loss_aggregation`; every entry defaults to the historical
    ## rule, so a model built and fitted without the argument behaves as before.
    self$loss_agg = NULL
    ## Observation tensors for the non-scalar aggregations, built once in fit():
    ##   class_obs  [sites, obs_years, species, K]  stem counts per diameter class
    ##   tree_obs   list of [sites, obs_years, patches, max_trees] tensors
    ## They are kept whole and sliced by the batch index, like init_cohort.
    ##   quant_obs  [sites, obs_years, species, n_probs]  observed dbh quantiles
    self$class_obs = NULL
    self$tree_obs = NULL
    self$quant_obs = NULL
    self$class_breaks = NULL
    self$class_tau = 1.0
    self$quant_probs = NULL
    self$quant_tau = 0.05
    self$record_raws = FALSE
    self$env_scaling = NULL
    self$train_env           = NULL
    self$train_init_cohort   = NULL
    self$fit_id              = NULL
    self$conditional_effects = NULL
    self$ale                 = NULL
    self$perm_importance     = NULL
    # Set by fit(). With weights = "auto", loss_baseline holds each response's
    # intercept-only (null) loss and loss_weights is 1/that - so every term in
    # `history` reads as a fraction of its own baseline.
    self$loss_weights        = NULL
    self$loss_baseline       = NULL
    # How often each response produced a non-finite loss during the last fit()
    # and was therefore left out of that timestep's gradient (see finite_terms).
    self$nonfinite_counts    = NULL
    private$add_process(mortality_process, "mortality")
    private$add_process(growth_process, "growth")
    private$add_process(regeneration_process, "regeneration")
    private$add_process(competition_process, "competition")

  },

  #' @description
  #' Predicts the growth, mortality, and regeneration of trees based on the given inputs.
  #'
  #' The `forward` method generates predictions for tree growth, mortality, and regeneration for the specified species across different environmental conditions. It uses the initialized model parameters and can handle optional input tensors like diameter at breast height (dbh), number of trees, and species. If these are not provided, they will be initialized internally.
  #'
  #' @param dbh torch.Tensor (Optional). Diameter at breast height of the trees.
  #' @param trees torch.Tensor (Optional). Number of trees.
  #' @param species torch.Tensor (Optional). Species of the trees.
  #' @param env torch.Tensor. Environmental data.
  #' @param y torch.Tensor. Response tensor for target data.
  #' @param disturbance torch.Tensor. Disturbance rates.
  #' @param start_time integer. Time at which to start recording the results.
  #' @param pred_growth torch.Tensor (Optional). Predicted growth values.
  #' @param pred_mort torch.Tensor (Optional). Predicted mortality values.
  #' @param pred_reg torch.Tensor (Optional). Predicted regeneration values.
  #' @param patches numeric. Number of patches.
  #' @param debug logical. Run in debug mode if TRUE.
  #' @param return_cohorts Which timesteps' cohort states to store: `FALSE` (none),
  #'   `TRUE`/`"all"` (every timestep), `"last"`, or an integer vector of timesteps.
  #' @param update_step integer. Backpropagation step length.
  #' @param verbose logical. Print progress if TRUE.
  #' @param year_sequence at which year indices should the predictions compared with the observed values
  #' @param class_obs torch.Tensor (Optional). Observed size-class stem counts,
  #'   `[sites, observed years, species, K + 1]`, used when `trees` is scored
  #'   with `loss_aggregation = "size_classes"`. Built by `fit()`.
  #' @param tree_obs list (Optional). Observed individual trees as padded tensors
  #'   (`dbh`, `species`, `trees`, `growth`, `died`), each
  #'   `[sites, observed years, patches, max trees]`, used when `growth` or
  #'   `mortality` is scored with `loss_aggregation = "tree"`. Built by `fit()`.
  #' @param nonfinite character. `"drop_term"` backpropagates the finite loss terms
  #'   of a timestep and drops the others; `"skip_step"` discards the whole
  #'   timestep as soon as one term is non-finite. See [FINN::fit].
  #' @param reanchor logical. Multiple shooting: replace the simulated cohort state
  #'   by the observed one at every timestep that has an entry in `anchor_cohorts`.
  #' @param anchor_cohorts named list. Names are timesteps (as characters), entries
  #'   carry `$dbh`, `$trees` and `$species` (a [CohortMat] or a plain list of the
  #'   three tensors, already subset to the sites of this call).
  #' @param anchor_mode character. `"replace"` (default) overwrites the simulated
  #'   state with the anchor; `"merge"` keeps the simulated cohorts below
  #'   `anchor_min_dbh` and appends the anchor to them.
  #' @param anchor_min_dbh numeric. The dbh below which the anchors are blind
  #'   (same units as `dbh`). Required by `anchor_mode = "merge"`.
  #' @param teacher_forcing logical. Direct estimation: reset the cohort state at
  #'   every timestep to the last observed stand (the incoming cohorts, then
  #'   `anchor_cohorts`), so no trajectory is ever simulated. The `dbh`/`ba`/`trees`
  #'   loss terms are dropped, because scoring them would need one.
  #' @param loss_resolution character. Score growth and mortality per (site, year,
  #'   species) cell ("site") or per individual tree ("tree"). See [FINN::fit].
  #' @param tree_data named list. Per-tree targets keyed by simulation timestep,
  #'   each a list of two tensors `growth` and `mort` with the dimensions of the
  #'   cohort state at the start of the interval, already subset to this call's sites.
  #'   Only used with `loss_resolution = "tree"`.
  #' @return list. A list of predicted values for dbh, number of trees, and other recorded time points. If `debug` is TRUE, raw results and cohorts are also returned.
  forward = function(dbh = NULL,
                     trees = NULL,
                     species = NULL,
                     env = NULL,
                     y = NULL,
                     disturbance = NULL,
                     start_time = 1L,
                     pred_growth = NULL,
                     pred_mort = NULL,
                     pred_reg = NULL,
                     patches = 100L,
                     debug = FALSE,
                     return_cohorts = FALSE,
                     update_step = 1L,
                     verbose = TRUE,
                     year_sequence = NULL,
                     class_obs = NULL,
                     quant_obs = NULL,
                     tree_obs = NULL,
                     nonfinite = c("drop_term", "skip_step"),
                     reanchor = FALSE,
                     anchor_cohorts = NULL,
                     anchor_mode = "replace",
                     anchor_min_dbh = NULL,
                     teacher_forcing = FALSE,
                     loss_resolution = c("site", "tree"),
                     tree_data = NULL){

    nonfinite = match.arg(nonfinite)
    anchor_mode = match.arg(anchor_mode, c("replace", "merge"))
    loss_resolution = match.arg(loss_resolution)
    if(isTRUE(reanchor) && isTRUE(teacher_forcing)) stop("`reanchor` and `teacher_forcing` are different estimators and cannot be combined: re-anchoring free-runs each segment from an observed state, teacher forcing never simulates one.")

    # if no cohorts exist initialize empty cohort array
    if(is.null(dbh)){
      cohorts = CohortMat(dims = c(env$shape[1], patches, self$sp), sp = self$sp)
      cohorts$to(device = self$device) # TODO: device must be set
      trees = cohorts$trees
      species = cohorts$species
      dbh = cohorts$dbh
    }

    # # repeat env if same env for the three processes
    if(is.list(env)) {
      env = lapply(env, function(e) torch_tensor(e, dtype=self$dtype, device=self$device))
    } else {
      # processes 1, 2, 3 are  mortality, growth, and regeneration
      env = lapply(1:3, function(i) torch_tensor(env, dtype=self$dtype, device=self$device))
    }
    names(env) <- c("mort", "growth", "reg")

    # get dimensions
    sites = env[[1]]$shape[1]
    time =  env[[1]]$shape[2]
    patches = dbh$shape[2]
    sp = self$N_species

    # Which timesteps' cohorts to store (see resolve_record_years()). `debug`
    # keeps recording every step; `return_cohorts` is the lightweight selector.
    record_years = resolve_record_years(return_cohorts, time, debug)
    record_cohorts = length(record_years) > 0

    # check dtype and device of disturbance
    if(!is.null(disturbance)) {
      disturbance = disturbance$to(dtype=self$dtype, device=self$device)
      disturbances_tens = torch::distr_bernoulli(probs = disturbance$squeeze(3L))$sample(patches)$permute(c(2, 3, 1))
      disturbances_tens = 1*(disturbances_tens==0)
    }

    # if y (response) is null -> simulation mode, gradients are not required (not neccessary to set them to zero but it is cleaner)
    if(is.null(y)) {
      lapply(self$parameters, function(p) p$requires_grad_(FALSE) )
      # make sure parameters are always re-enabled for gradient tracking, even if
      # an error is thrown later in this function (e.g. the "Light values > 1" guard
      # below) - otherwise the model is left permanently frozen and a later fit()
      # call will silently fail to learn.
      on.exit(lapply(self$parameters, function(p) p$requires_grad_(TRUE)), add = TRUE)
    }

    # send data to correct devices (and dtypes) TODO: unncessary?
    dbh = torch_tensor(dbh, dtype=self$dtype, device=self$device)
    trees = torch_tensor(trees, dtype=self$dtype, device=self$device)
    species = torch_tensor(species, dtype=torch_int64(), device=self$device)

    # Cohort ids, unique within a (site, patch): a cohort is traced by the triple
    # (siteID, patchID, cohortID).
    #
    # At t = 0 a cohort's id IS its position along the cohort dimension, which is
    # the `cohortID` obsDF2arrays() placed it at. So an initial cohort keeps the
    # identity `makeInitCohorts()` gave it -- and since makeInitCohorts() called
    # WITHOUT `dbh_binsize` does not bin (one row in, one cohort out), that is
    # the identity of the observed tree the row came from. The trace back to a
    # tagged inventory tree is then just a join on (siteID, patchID, cohortID),
    # with no assumption about how the arrays happen to be laid out.
    cohort_ids = torch_tensor(array(
      rep(seq_len(species$shape[3]), each = prod(species$shape[1:2])),
      dim = species$shape), dtype=torch_int32(), device = self$device
    )
    # Ids of new cohorts come from a counter that only ever goes up. Taking
    # max(cohort_ids) instead would re-issue the id of a cohort that has died and
    # been pruned, and a trace would then silently follow a recruit that inherited
    # a dead tree's id.
    next_cohort_id = as.integer(max(c(0L, species$shape[3])))

    # observed cohorts -> a state the time loop can run on: right device and dtypes, and ids
    # in the scheme above - cohort j of every (site, patch) gets `base_id + j` - so the caller
    # takes `base_id` from `next_cohort_id` and advances the counter past the block
    observed_state = function(anchor, base_id) list(
      dbh     = anchor$dbh$to(dtype = self$dtype, device = self$device),
      trees   = anchor$trees$to(dtype = self$dtype, device = self$device),
      species = anchor$species$to(dtype = torch_int64(), device = self$device),
      ids     = torch_tensor(array(rep(base_id + seq_len(anchor$species$shape[3]),
                                       each = prod(anchor$species$shape[1:2])),
                                   dim = anchor$species$shape),
                             dtype = torch_int32(), device = self$device))

    # Teacher forcing: the rates observed over an interval were generated by the stand at the
    # census that STARTS it, so the state is reset to that observation at every timestep and
    # the simulated state is never carried forward. The incoming cohorts are the year-0 census.
    forced = if(isTRUE(teacher_forcing)) list(dbh = dbh, trees = trees, species = species, ids = cohort_ids) else NULL

    ## Tree-level loss: setup ####
    # Within a generation - the block handed out at t0 or at an anchor - cohort j of a
    # (site, patch) carries id `base + j` and is observed tree j of that (site, patch), so its
    # target sits at the column-major flat index of slot (s, p, j). Recruits and later anchors
    # get ids above the block, so the tracked set is `base < id <= base + n_cohorts` and a
    # recruit cannot enter a residual by construction.
    tree_mode = identical(loss_resolution, "tree") && !is.null(tree_data) && !is.null(y)
    tree_gen_base = 0L
    tree_gen_n = as.integer(prod(species$shape))
    tree_gen_nc = as.integer(species$shape[3])
    # 0-based offset of (s, p) in the flat index, and the stride of one cohort position
    tree_sp_n = as.integer(sites * patches)
    tree_sp_off = torch_tensor(array((seq_len(sites) - 1L) + rep((seq_len(patches) - 1L) * sites, each = sites),
                                     dim = c(sites, patches, 1L)),
                               dtype = torch_long(), device = self$device)
    tree_acc = NULL
    tree_cov = list()
    tree_flat = function(a) a$to(dtype = self$dtype, device = self$device)$permute(c(3,2,1))$reshape(c(-1))
    tree_new = function(n) lapply(1:5, function(k) torch::torch_zeros(n + 1L, device = self$device))
    # Slot-wise zeroing keeps the graph of the sites still inside their window (the per-site
    # path backwards once, at the end); elsewhere that graph is already freed, so start fresh.
    # A generation's slots are the cohort layout it was built from, padding included.
    # Only the slots that held a real tree are tracked; the empty ones would otherwise
    # accumulate a rate against a NaN target and pollute `$tree_deaths`.
    tree_live0 = NULL
    tree_slot_idx = function(ids = cohort_ids) {
      loc = ids$to(dtype = torch_long()) - tree_gen_base
      inr = loc$ge(1)$bitwise_and(loc$le(tree_gen_nc))
      sl  = tree_sp_off + (loc - 1L)*tree_sp_n + 1L
      torch_where(inr, sl, torch_full_like(sl, tree_gen_n + 1L))
    }
    tree_mark_live = function(tr = trees, ids = cohort_ids) {
      idx = tree_slot_idx(ids)$flatten()
      v = tr$gt(0.5)$to(dtype = torch_float32())$flatten()
      torch::torch_zeros(tree_gen_n + 1L, device = self$device)$scatter_add(1, idx, v)$gt(0.5)
    }
    tree_tracked = function() {
      idx = tree_slot_idx()
      tree_live0$to(dtype = torch_float32())$gather(1, idx$flatten())$reshape(idx$shape)$gt(0.5)
    }
    tree_zero = function(acc, reset) {
      if(!any(reset > 0.5)) return(acc)
      if(all(reset > 0.5)) return(tree_new(tree_gen_n))
      k = torch_tensor(c(rep(1 - reset, times = tree_gen_n %/% sites), 1), dtype = self$dtype, device = self$device)
      lapply(acc, function(a) a*k)
    }
    if(tree_mode) {
      tree_acc = tree_new(tree_gen_n)
      tree_live0 = tree_mark_live()
      tree_data = lapply(tree_data, function(td) lapply(td, tree_flat))
      # The accumulators cover (i - period_length + 1):i, the window the aggregate rate terms
      # use, so they are zeroed at each window's first year - zeroing at each census instead
      # would hand FIA's first interval an extra year of spin-up.
      tree_starts = matrix(0, nrow = sites, ncol = time)
      pl_obs = as.matrix(y[,,1,7]$cpu())
      for(t in seq_along(year_sequence)) {
        w = pl_obs[,t]
        w[!is.finite(w)] = 1
        tree_starts[cbind(seq_len(sites), pmax(1, year_sequence[t] - w + 1))] = 1
      }
    }

    # init Result tensors
    Result = lapply(1:7,function(tmp) torch::torch_zeros(list(sites, time, self$N_species), device=self$device))
    names(Result) =  c("dbh","ba", "trees", "growth", "mort", "reg", "r_mean_ha")

    # running aggregate of the per-update_step loss, used so that the loss reported
    # back to fit() reflects all observation points hit during this forward() call,
    # not just whichever one happened to be evaluated last in the time loop.
    loss_total = torch::torch_zeros(7L, device = self$device)
    loss_count = 0L

    # Per-site aggregation windows (set by fit() when `period_length` differs
    # between sites observed in the same data-year, or when sites are observed
    # in different years). Each site's growth/mortality/regeneration prediction
    # is then aggregated over ITS OWN window (i - period_s + 1):i, and because
    # the windows of different sites overlap across steps the backward is
    # deferred to the end of the time loop instead of running at every
    # observation step (which would free the graph other sites still need).
    per_site_mode = !is.null(y) && isTRUE(self$per_site_period)
    loss_deferred = torch::torch_zeros(1L, device = self$device)
    # Per-site path: every response is ONE mean over all units (trees, cells,
    # tree counts) scored in this call, not a sum of per-step means. Each loss
    # function returns the mean over the units of its step, so the steps are
    # kept with their unit counts (`loss_n`) and weighted n_step / n_call at the
    # end of the time loop. Summing per-step means gave a step with ten trees the
    # weight of a step with thousands whenever sites are observed in different
    # years (FINN-bwi variant B: 2.9 % of the mortality objective for the
    # drought half of the split corners, 74.5 % for the pre-drought half). With
    # one scored step per call every weight is n/n = 1 and nothing changes.
    step_losses = list()
    scored_n = function(true, n = NULL) if (per_site_mode) private$scored_units(true, n) else 0

    # storage for the raw per-timestep state. Cohort snapshots are kept only for
    # the requested timesteps (`record_cohorts`); patch snapshots stay tied to the
    # legacy `debug` path. Entries are indexed by their true timestep, so skipped
    # steps stay NULL and pred2DF() labels every recorded step with its real year.
    if(record_cohorts) {
      Raw_cohort_results = list()
      Raw_cohort_ids = list()
    }
    if(debug) {
      Raw_patch_results = list()
    }

    # if simulation mode, environmental predictions can be made a priori (it would interrupt the gradients in inference mode)
    if(is.null(y)) {
      if(!inherits(self$process_growth, "hybrid")) predGrowthGlobal = self$nn_growth(env[["growth"]])
      if(!inherits(self$process_mortality, "hybrid")) predMortGlobal = self$nn_mortality(env[["mort"]])
      if(!inherits(self$process_regeneration, "hybrid")) predRegGlobal = self$nn_regeneration(env[["reg"]])
    }

    # create process bar
    if(verbose) cli::cli_progress_bar(format = "Year: {cli::pb_current}/{cli::pb_total} {cli::pb_bar} ETA: {cli::pb_eta} ", total = time, clear = FALSE)
    for(i in 1:time){

      if(isTRUE(teacher_forcing)) {
        dbh = forced$dbh
        trees = forced$trees
        species = forced$species
        cohort_ids = forced$ids
      }

      if(tree_mode) tree_acc = tree_zero(tree_acc, tree_starts[,i])

      # store this timestep's cohorts?
      rec_i = record_cohorts && (i %in% record_years)

      # In inference mode, make env predictions in each time step (to get the gradients)
      # otherwise, just take the i-th prediction
      if(!is.null(y)) {
        if(!inherits(self$process_growth, "hybrid")) pred_growth = self$nn_growth(env[["growth"]][,i,])
        if(!inherits(self$process_mortality, "hybrid")) pred_mort = self$nn_mortality(env[["mort"]][,i,])
        if(!inherits(self$process_regeneration, "hybrid")) pred_reg = self$nn_regeneration(env[["reg"]][,i,])
      } else {
        if(!inherits(self$process_growth, "hybrid")) pred_growth = predGrowthGlobal[,i,]
        if(!inherits(self$process_mortality, "hybrid")) pred_mort = predMortGlobal[,i,]
        if(!inherits(self$process_regeneration, "hybrid")) pred_reg = predRegGlobal[,i,]
      }

      # empty rate objects/tensors, [sites, patches, cohorts] like the state they
      # stand in for when no cohort is alive (they used to be sized with `time`
      # in place of `patches`, which broke the cohort record of an empty stand)
      light = torch_zeros(list(sites, patches, dbh$shape[3]), device=self$device)
      g = torch_zeros(list(sites, patches, dbh$shape[3]), device=self$device)
      m = torch_zeros(list(sites, patches, dbh$shape[3]), device=self$device)
      r = torch_zeros(list(sites, patches, dbh$shape[3]), device=self$device)
      if(rec_i) trees_before = torch::torch_zeros_like(g)

      # detach previous cohort objects (to interrupt the gradients)

      # if(!is.null(y)) {
      #   #for(j in 1:3) Result[[j]] = Result[[j]]$detach()
      #   # if period_length = NA on all sites, it is assumed that we need only yearly gradients
      #   if(as.numeric(y[,,,7]$isnan()$bitwise_not()$sum()) < 0.5) {
      #     dbh=dbh$detach()
      #     trees=trees$detach()
      #     species=species$detach()
      #     cohort_ids=cohort_ids$detach()
      #   }
      #
      # }

      # dbh=dbh$detach()
      # trees=trees$detach()
      # species=species$detach()
      # cohort_ids=cohort_ids$detach()

      # Apply disturbance
      if(!is.null(disturbance)) {
        trees = trees*disturbances_tens[,i,]$unsqueeze(3L)
      }

      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      ## Demographic processes ####
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      # check if there is at least one cohort alive (otherwise just jump directly to regeneration?)
      if(dbh$shape[3] > 0.5){
        # calculate available light for each cohort
        ### Competition 1 ####
        light = self$competition_func(
          dbh = dbh,
          species = species,
          trees = trees,
          # parHeight = parHeight,
          parComp = self$par_competition,
          h = NULL,
          # minLight = self$minLight,
          patch_size_ha = self$patch_size_ha,
          ba = NULL,
          cohortHeights = NULL,
          n_quantiles = self$n_quantiles,
          continuous = self$continuous,
          light_mask = isTRUE(self$light_mask)
        )

        if (light$gt(1.0)$sum()$item() > 0) {
          stop("Light values > 1")
        }

        ## Growth ####
        if(!inherits(self$process_growth, "hybrid")) pred = index_species(pred_growth, species)
        else pred = env[[1]][,i,]

        g = self$growth_func(
          dbh = dbh,
          species = species,
          parGrowth = self$par_growth,
          pred = pred,
          light = light,
          trees = trees
        )
        self$g = g
        dbh_growth = dbh*g
        dbh = dbh + dbh_growth

        ## Competition 2 ####
        light = self$competition_func(
          dbh = dbh,
          species = species,
          trees = trees,
          parComp = self$par_competition,
          h = NULL,
          # minLight = self$minLight,
          patch_size_ha = self$patch_size_ha,
          ba = NULL,
          cohortHeights = NULL,
          n_quantiles = self$n_quantiles,
          continuous = self$continuous,
          light_mask = isTRUE(self$light_mask)
        )
        # cat("Second section  B2\n")

        ## Mortality ####
        if(!inherits(self$process_mortality, "hybrid")) pred = index_species(pred_mort, species)
        else pred = env[[2]][,i,]

        m = self$mortality_func(
          dbh = dbh,
          species = species,
          trees = trees + 0.001,
          parMort = self$par_mortality,
          pred = pred,
          light = light,
          growth = g
        )

        n_at_risk = torch::torch_clamp(trees+trees$le(0.5)$float()+0.01, min = 1.0)
        # The floor guards the Gamma draw, whose shape n*p/(1-p) degenerates as p -> 0.
        # NULL (a model saved before mortality_floor existed) keeps the published 0.01.
        m_draw = torch::torch_clamp(m,
                                    if(is.null(self$mortality_floor)) 0.01 else self$mortality_floor,
                                    if(is.null(self$mortality_ceiling)) 0.99 else self$mortality_ceiling)
        # NULL (a model saved before mortality_draw existed, or a hybrid built elsewhere) is the default
        alive_mask = trees$ge(0.5)$float()
        # "binomial" resolves by whether the draw has to carry a gradient at all.
        # With rb_gradient the draw is detached and the derivative comes from the
        # expectation, so torch's exact binomial is free to be non-differentiable;
        # without it the gradient reaches the mortality parameters THROUGH the state
        # (deaths -> trees -> ba), and the draw must be differentiable, where the
        # moment-matched Gamma gives -39.36 against the expectation's -39.32 while
        # the straight-through Bernoulli gives -22.59. See FINNverse/FINN#62.
        trees_dead = (switch(if(is.null(self$mortality_draw)) "gamma" else self$mortality_draw,
                             binomial  = if(isTRUE(self$rb_gradient)) torch_binomial_exact(n_at_risk, m_draw)
                                         else binomial_from_gamma(n_at_risk, m_draw),
                             bernoulli = binomial_from_bernoulli(n_at_risk, m_draw),
                             poisson   = binomial_from_poisson(n_at_risk, m_draw),
                             gamma     = binomial_from_gamma(n_at_risk, m_draw)))*alive_mask
        trees_dead = trees_dead + trees_dead$round()$detach() - trees_dead$detach()
        # `rb_gradient`: keep the drawn, rounded count in the FORWARD pass -- the stand must
        # die discretely -- but take the derivative through the expectation n_at_risk * m.
        # The reparameterised draw and the rounding's surrogate both carry the sampler's
        # noise into the gradient; measured on BWI, 97 % of the mortality gradient is draw
        # noise (SNR 0.52) against 43 without it. The forward value is unchanged to the bit,
        # so simulations and any fit with rb_gradient = FALSE reproduce exactly.
        if (isTRUE(self$rb_gradient)) {
          dead_exp   = n_at_risk * m_draw * alive_mask
          trees_dead = dead_exp + (trees_dead - dead_exp)$detach()
        }
        trees_before = trees
        trees = torch_clamp(trees - trees_dead, min = 0) #### TODO if trees = 0 then NA...prevent!

        ## Tree-level loss: accumulate this year's rates per tracked cohort ####
        # A tracked cohort accumulates for the WHOLE window whether or not the simulation
        # killed it: the observation it is scored against is the mean annual rate and the
        # 0/1 fate over the whole interval, so a truncated accumulation would compare two
        # different quantities. A killed slot keeps `trees = 0`, shades nothing and
        # contributes nothing to the stand terms, but its diameter and hazard stay defined.
        # Accumulators 4 and 5 are the two liveness counts `$tree_deaths` needs, and one of
        # them is not enough. 4 counts the years the cohort ENTERED alive, so a cohort killed
        # in year k of an L-year window gets k and a cohort that was already dead at the
        # window's start gets 0 -- but a kill in year L is then indistinguishable from a
        # survivor. 5 counts the years it ENDED alive, which separates those two (k-1 against
        # L) but gives 0 both for a year-1 kill and for a ghost. Survivor: 5 == n_ok.
        # Ghost: 4 == 0. Killed in this window: 4 > 0 and 5 < n_ok, in year 4.
        if(tree_mode) {
          keep = tree_tracked()$to(dtype = torch_float32())
          idx = torch_where(keep$gt(0.5), tree_slot_idx(),
                            torch_full_like(cohort_ids$to(dtype = torch_long()), tree_gen_n + 1L))$flatten()
          gi = if(isTRUE(self$growth_period_scale)) torch_log1p(g) else g
          tree_acc[[1]] = tree_acc[[1]]$scatter_add(1, idx, (gi*keep)$flatten())
          tree_acc[[2]] = tree_acc[[2]]$scatter_add(1, idx, (torch_log((1-m)$clamp(1e-6, 1))$negative()*keep)$flatten())
          tree_acc[[3]] = tree_acc[[3]]$scatter_add(1, idx, keep$flatten())
          tree_acc[[4]] = tree_acc[[4]]$scatter_add(1, idx, (keep*trees_before$gt(0.5)$to(dtype = torch_float32()))$flatten())
          tree_acc[[5]] = tree_acc[[5]]$scatter_add(1, idx, (keep*trees$gt(0.5)$to(dtype = torch_float32()))$flatten())
        }
      }

      ### Competition 3 ####
      # start reg
      AL_reg = self$competition_func( # must have dimension = n species in last dim
        dbh = dbh,
        species = species,
        trees = trees,
        parComp = self$par_competition,
        h = 1,
        # minLight = self$minLight,
        patch_size_ha = self$patch_size_ha,
        ba = NULL,
        cohortHeights = NULL,
        n_quantiles = self$n_quantiles,
        continuous = self$continuous,
        light_mask = isTRUE(self$light_mask)
      )

      # browser()
      ### Regeneration ####
      if(!inherits(self$process_regeneration, "hybrid"))pred = pred_reg
      else pred = env[["reg"]][,i,]

      ## Density dependence (e.g. a Beverton-Holt cap) is now the responsibility
      ## of the regeneration process function itself -- see FINN::regeneration
      ## (unbounded) vs FINN::regeneration_saturation (capped, with its own
      ## reg_logK parameter declared via createProcess(custom_parameters=)).
      ## The stand itself (dbh, trees) goes only to functions that ask for it,
      ## e.g. FINN::regeneration_adult (seed source = conspecific adults).
      reg_args = list(species = species,
                      parReg = self$par_regeneration[,1],
                      pred = pred,
                      light = AL_reg)
      if(all(c("dbh", "trees") %in% names(formals(self$regeneration_func))))
        reg_args = c(reg_args, list(dbh = dbh, trees = trees))
      r_mean_ha = do.call(self$regeneration_func, reg_args)

      r_mean_patch = r_mean_ha*self$patch_size_ha

      if(self$sample_regeneration){
        r = if(identical(self$regeneration_sampler, "refined")) rnbinom_torch_refined(r_mean_patch, self$par_theta_recruits)
            else rnbinom_torch(r_mean_patch, self$par_theta_recruits)
      }else if(!self$sample_regeneration){
        r = r_mean_patch
      }

      r = r + r$round()$detach() - r$detach()
      # the same split for recruitment: drawn forward, expectation backward (see the
      # death draw above). The regeneration loss already scores r_mean_ha, so the
      # recruitment parameters keep an analytic gradient of their own either way.
      if (isTRUE(self$rb_gradient)) r = r_mean_patch + (r - r_mean_patch)$detach()

      # New recruits
      dbh_new = ((r-1+0.1)/1e-3)$sigmoid() + (self$recruits_dbh - 1.0) # TODO: check!!! --> when r 0 dann dbh = 0, ansonsten dbh = 1 dbh[r==0] = 0
      trees_new = r
      species_new = torch_arange(1, sp, dtype=torch_int64(), device = self$device)$unsqueeze(1)$`repeat`(c(r$shape[1], r$shape[2], 1))

      # assign cohortIDs from the monotone counter (see where it is initialised)
      n_new = prod(r$shape)
      new_cohort_id = torch_tensor(array(
        next_cohort_id + seq_len(n_new),
        dim = r$shape), dtype=torch_int32(), device = self$device
      )
      next_cohort_id = next_cohort_id + n_new

      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      ## Aggregation of rates####
      # from previous cohort
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      if(dbh$shape[3] != 0){
        samples = vector("list", 3)
        names(samples) = c(
          "g*trees_before",
          "m*trees_before",
          "trees_before"
        )
        # Result[[4]] is the tree-weighted mean of g, matching an observed
        # mean(dbh/dbh_before - 1). `g` is the model's native parameter
        # (dbh_new = dbh*(1+g)), which is why the relative rate is the right
        # target here: training on the absolute increment (dbh*g) instead was
        # tested over 3 seeds and came out WORSE on every one, even when scored
        # on the absolute scale, and ~4x more seed-variable - it couples the
        # growth parameter to whichever trees happen to be present.
        samples[[1]] = g*trees_before
        samples[[2]] = m*trees_before
        samples[[3]] = trees_before # original number of trees

        Results_tmp = replicate(length(samples), torch_zeros_like(Result[[1]][,i,]))
        # calculate the sum of all elements in samples over all patches for each species, patch and site
        tmp_res = aggregate_results(species, samples, Results_tmp, aggregation = "sum")

        # Aggregation [sites, patches, cohorts] --> [sites, species]:
        if(as.numeric(trees$sum() > 0)) {
          mask = tmp_res[[3]]$gt(0.5)
          # Growth
          Result[[4]][,i,][mask] = Result[[4]][,i,][mask]+tmp_res[[1]][mask]/tmp_res[[3]][mask] # summe dbh_growth / summe trees  #/cohort_counts[alive_species])
          # Mort
          Result[[5]][,i,][mask] = Result[[5]][,i,][mask]+tmp_res[[2]][mask]/tmp_res[[3]][mask] # summe dbh_growth / summe trees  #/cohort_counts[alive_species])
        }
      }
      # reg extra
      ## Regeneration count
      tmp_res = aggregate_results(species_new, list(r), list(torch::torch_zeros(Result[[1]][,i,]$shape[1], sp, device = self$device )))
      Result[[6]][,i,] = Result[[6]][,i,]$add(tmp_res[[1]])/patches

      ## Regeneration rate mean
      tmp_res = aggregate_results(species_new, list(r_mean_ha), list(torch::torch_zeros(Result[[1]][,i,]$shape[1], sp, device = self$device )))
      r_mean_ha = tmp_res[[1]]/patches
      Result[[7]][,i,] = r_mean_ha
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      ## Update arrays ####
      # 1. Combine old cohorts and recruit cohorts
      # 2. Find dead cohorts and remove them, if possible (by finding the minimal required cohort dimension),
      #    a few dead cohorts are preserved for padding.
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

      dbh = torch::torch_cat(list(dbh, dbh_new), 3)
      trees = torch::torch_cat(list(trees, trees_new), 3)
      species = torch::torch_cat(list(species, species_new), 3)
      cohort_ids = torch::torch_cat(list(cohort_ids, new_cohort_id), 3)

      if (rec_i) {
        tmp = torch::torch_zeros_like(dbh_new)
        tmp[] = NaN
        g = torch::torch_cat(list(g, tmp), 3)
        m = torch::torch_cat(list(m, tmp), 3)
        trees_before = torch::torch_cat(list(trees_before, tmp), 3)
      }
      rm(dbh_new,trees_new,species_new,new_cohort_id )

      # Pad tensors, expensive, currently each timestep
      if(i %% 1 == 0){

        # Gradient shouldn't be required, also expensive for backpropagation because of reshape/view operations!
        #torch::with_no_grad({
        # Masks to find alive cohorts. A tracked cohort is kept even when the simulation
        # killed it, so its diameter and hazard stay defined to the end of the interval.
          alive_now = trees > 0.5
          if(tree_mode) alive_now = alive_now$bitwise_or(tree_tracked())
          mask = alive_now$flatten(start_dim = 1, end_dim = 2)
          org_dim = species$shape[1:2]
          org_dim_t = torch::torch_tensor(org_dim, dtype = torch_long(), device = "cpu")

          # Minimal number of cohort dimension
          non_zero_counts = mask$sum(dim=2)
          max_non_zeros = non_zero_counts$max()$item()
          sorted_tensor = torch::torch_sort(mask$float(), dim=2, descending=TRUE)[[2]]

          dbh = dbh$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
          trees = trees$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
          cohort_ids = cohort_ids$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
          species = species$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)

          if(rec_i) {
            g = g$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
            m = m$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
            trees_before = trees_before$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
          }

        #})
      }

      # aggregate results
      # TODO: Position of the first block?
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      ## Aggregation of stand variables ####
      # from updated cohorts
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

      if(i > 0){
        if(dbh$shape[3] != 0){
          #dead_trees_mask = trees == 0
          if(tree_mode) dbh = dbh*trees$gt(0.5)$bitwise_or(tree_tracked())$float()
          else dbh = dbh*trees$gt(0.5)$float()
          BA_stem_values = BA_stem(dbh = dbh)*trees
          species = species
          samples = vector("list", 3)
          samples[[1]] = trees * dbh #* mask
          samples[[2]] = BA_stem_values #* mask# torch_sigmoid((trees - 0.5)/1e-3) # better to just use greater? (Masking!) Gradients shouldn't be needed! (I think?)
          samples[[3]] = trees # torch_sigmoid((trees - 0.5)/1e-3)
          Results_tmp = replicate(length(samples), torch_zeros_like(Result[[1]][,i,]))

          tmp_res = aggregate_results(species, samples, Results_tmp)
          # BA and number of trees Result[[1]] and Result[[2]]
          if(as.numeric(trees$sum() > 0)) {
            # mask only alive trees
            mask = tmp_res[[3]]$gt(0.5)
            Result[[1]][,i,][mask] = Result[[1]][,i,][mask]+ (tmp_res[[1]][mask]/tmp_res[[3]][mask]) # /cohort_counts[alive_species]
            # BA
            Result[[2]][,i,]= Result[[2]][,i,]$add(tmp_res[[2]]$div_(patches))
            # Trees
            Result[[3]][,i,]= Result[[3]][,i,]$add(tmp_res[[3]]$div_(patches))
          }
          rm(BA_stem_values)
        }
      }
      # end aggregation

      if (rec_i) {
        # index by the true timestep `i`; unrecorded steps stay NULL
        Raw_cohort_results[[i]] = list(
          "species" = torch::as_array(species$cpu()),
          "trees" = torch::as_array(trees$cpu()),
          "dbh" = torch::as_array(dbh$cpu()),
          "m" = torch::as_array(m$cpu()),
          "g" = torch::as_array(g$cpu()),
          "trees_before" = torch::as_array(trees_before$cpu()) # noob
        )
        Raw_cohort_ids[[i]] = torch::as_array(cohort_ids$cpu())
      }
      if (debug) {
        Raw_patch_results[[i]] = list(
          "r" = torch::as_array(r$cpu())
        )
      }

      loss = torch_zeros(7L, device = self$device)
      loss_n = numeric(7L)   # units behind each term's mean at this step (per-site path)

      ##### Rework #######
      # NOTE: truncation of the recurrent cohort state (below) is intentionally
      # decoupled from `i %in% year_sequence` - see comment further down.
      update_boundary = i > 0 && dbh$shape[3] != 0 && !is.null(y) && (i %% update_step == 0)
      if(update_boundary) {
        if(i %in% year_sequence) {
          tmp_index = which(year_sequence %in% i, arr.ind = TRUE)
          # browser()
          # Under teacher forcing the state at step i is the LAST OBSERVED stand carried one
          # year, not a trajectory that reached year i, so the three state responses are not
          # predictions of this census at all and leave the objective.
          if(!isTRUE(teacher_forcing)) {
            # #dbh
            loss[1] = self$loss_dbh_func(y[, tmp_index,,1], Result[[1]][,i,] )
            loss_n[1] = scored_n(y[, tmp_index,,1])

            # ---- diameter quantiles replace the scalar mean diameter ----------
            # A mean cannot separate one exceptional tree from a shift of the whole
            # distribution; three quantiles state where the thin, the typical and
            # the thick stems of each species actually are. Unlike a count in a
            # size class this stays on the diameter scale, so the upper quantile
            # is a direct statement about the big trees.
            if (!is.null(quant_obs) && identical(self$loss_agg$dbh$type, "quantiles")) {
              q_pred = private$quantile_values(dbh, trees, species, self$quant_probs, self$quant_tau)
              loss[1] = self$loss_dbh_func(quant_obs[, tmp_index, , ], q_pred)
              loss_n[1] = scored_n(quant_obs[, tmp_index, , ])
            }
            # ba
            loss[2] = self$loss_ba_func(y[, tmp_index,,2], Result[[2]][,i,] )
            loss_n[2] = scored_n(y[, tmp_index,,2])
            # counts
            loss[3] = self$loss_trees_func(y[,tmp_index,,3], Result[[3]][,i,])
            loss_n[3] = scored_n(y[,tmp_index,,3])
          }

          # ---- size-classed stem counts replace the scalar `trees` term ------
          # The prediction is a PLOT total per species and diameter class, which
          # is what the observed integer counts are, so the Poisson applies
          # directly. Contains both abundance and size composition, because
          # independent Poisson class counts factorise into Poisson(total) x
          # Multinomial(shape).
          if (!isTRUE(teacher_forcing) && !is.null(class_obs) && identical(self$loss_agg$trees$type, "size_classes")) {
            cls_pred = private$class_counts(dbh, trees, species, self$class_breaks, self$class_tau)
            loss[3] = self$loss_trees_func(class_obs[, tmp_index, , ], cls_pred)
            loss_n[3] = scored_n(class_obs[, tmp_index, , ])
          }

          # growth rates - check for NA in period_length, if not, then accumulate gradients?
          # constant-period path below assumes the same period on every site
          accumulate_gradients = y[,tmp_index,,7] |> as.matrix()
          period = unique(accumulate_gradients[,1])

          ## Tree resolution: score the closed window, in place of the cell terms ####
          # Every tracked cohort covers the whole window by construction (see the accumulator
          # above), so neither term selects its own sample any more. Mortality scores the
          # compounded per-step probability, never the sampled fate, which at trees = 1 is a
          # coin flip. `$tree_deaths` reports how often the simulation killed a tree the
          # inventory found alive, and in which year - a diagnostic, not a correction.
          if(tree_mode && !is.null(tree_data[[as.character(i)]])) {
            td = tree_data[[as.character(i)]]
            w = accumulate_gradients[,1]
            w[!is.finite(w)] = 1
            w_flat = torch_tensor(rep(w, times = tree_gen_n %/% sites), dtype = self$dtype, device = self$device)
            n_ok = tree_acc[[3]]$narrow(1, 1, tree_gen_n)
            n_entered = tree_acc[[4]]$narrow(1, 1, tree_gen_n)
            n_alive = tree_acc[[5]]$narrow(1, 1, tree_gen_n)
            g_acc = tree_acc[[1]]$narrow(1, 1, tree_gen_n)
            g_pred = if(isTRUE(self$growth_period_scale)) torch_expm1(g_acc) else g_acc/n_ok$clamp(min = 1)
            p_death = torch_expm1(tree_acc[[2]]$narrow(1, 1, tree_gen_n)$negative())$negative()
            # fast-setup (9c03a90) keeps a killed tracked cohort to the end of its
            # window, so no target is hidden any more; F1 (5badcda) counts the units
            # behind each term's mean for the per-site pooling after the time loop.
            loss[4] = self$loss_growth_func(td$growth, g_pred)
            loss[5] = self$loss_mortality_func(td$mort, p_death$clamp(1e-6, 1 - 1e-6))
            loss_n[4] = scored_n(td$growth)
            loss_n[5] = scored_n(td$mort)
            self$tree_pred = list(growth = as.numeric(g_pred$detach()$cpu()),
                                  mort = as.numeric(p_death$detach()$cpu()),
                                  n = as.numeric(n_ok$detach()$cpu()),
                                  n_alive = as.numeric(n_alive$detach()$cpu()))
            # of the slots carrying a fate: how many the simulation killed before the window
            # closed, and how long those lived. A survivor ends every tracked year alive,
            # `n_alive == n_ok`; a kill in the window's LAST year is a survivor by the
            # entering count alone, which is why both counts are kept.
            # A generation outlives a window under single shooting, so a slot killed in an
            # earlier window enters this one already dead: it accumulates nothing and its
            # `n_entered` is 0. Counting it with the trees killed INSIDE this window would
            # mix two different statements and pull `mean_kill_year` toward zero, which is
            # what it did. `ghost` is that share; `killed` is this window's kills alone, and
            # the year it died is the number of years it entered alive.
            has = td$mort$isnan()$bitwise_not()$to(dtype = self$dtype)
            ghost = has*n_entered$lt(0.5)$to(dtype = self$dtype)
            killed = has*n_alive$lt(n_ok - 0.5)$to(dtype = self$dtype)*n_entered$gt(0.5)$to(dtype = self$dtype)
            tree_cov[[length(tree_cov) + 1L]] = c(i, as.numeric(torch_stack(list(
              has$sum(), killed$sum(), (killed*n_entered)$sum(),
              (has*n_ok)$sum(), ghost$sum()))$detach()$cpu()))
          }
          if(per_site_mode) {
            # ---- per-site windows: mask M[s, t] = 1 for t in (i - p_s, i] ----
            p_site = accumulate_gradients[,1]
            p_t = torch_tensor(ifelse(is.na(p_site), 0, p_site), dtype = self$dtype, device = self$device)
            if(any(!is.na(p_site))) {
              tt = torch_arange(1, i, dtype = self$dtype, device = self$device)
              M  = tt$unsqueeze(1)$gt((i - p_t)$unsqueeze(2))$to(dtype = self$dtype) *
                   p_t$gt(0.5)$to(dtype = self$dtype)$unsqueeze(2)
              M3 = M$unsqueeze(3)
              denom = p_t$clamp(min = 1)$unsqueeze(2)
              G = Result[[4]][,1:i,,drop = FALSE]
              g_pred = if (isTRUE(self$growth_period_scale)) (G*M3 + 1)$prod(2) - 1 else (G*M3)$sum(2)/denom
              m_pred = (Result[[5]][,1:i,,drop = FALSE]*M3)$sum(2)/denom
              r_pred = (Result[[7]][,1:i,,drop = FALSE]*M3)$sum(2)
              if(!tree_mode) {
                loss[4] = self$loss_growth_func(y[,tmp_index,,4], g_pred, y[,tmp_index,,9])
                loss[5] = self$loss_mortality_func(y[,tmp_index,,5], m_pred, y[,tmp_index,,8])
                loss_n[4] = scored_n(y[,tmp_index,,4], y[,tmp_index,,9])
                loss_n[5] = scored_n(y[,tmp_index,,5], y[,tmp_index,,8])
              }
              loss[6] = self$loss_regeneration_func(y[,tmp_index,,6], r_pred)
              loss_n[6] = scored_n(y[,tmp_index,,6])

              # ---- per-tree growth / mortality on the OBSERVED stand ---------
              # No correspondence is sought between model cohorts and observed
              # trees; the observed stand replaces the simulated one for these
              # terms, so light and the kernels are evaluated on the real trees
              # at their real diameters. The mse and binomial closures already
              # mask NaN, so padding slots and non-survivors drop out.
              if (!is.null(tree_obs)) {
                # tp$growth is the window MEAN of the annual rates (matching the
                # annualised observation) and tp$mortality the interval death
                # probability (matching the observed 0/1 over the interval).
                tp = private$tree_predictions(tree_obs, tmp_index, env, i, M, p_t)
                if (identical(self$loss_agg$growth$type, "tree")) {
                  loss[4] = self$loss_growth_func(tree_obs$growth[, tmp_index, , ], tp$growth)
                  loss_n[4] = scored_n(tree_obs$growth[, tmp_index, , ])
                }
                if (identical(self$loss_agg$mortality$type, "tree")) {
                  loss[5] = self$loss_mortality_func(tree_obs$died[, tmp_index, , ], tp$mortality)
                  loss_n[5] = scored_n(tree_obs$died[, tmp_index, , ])
                }
              }
              self$obs_rec = y[,tmp_index,,6] |> as.matrix()
              self$pred_rec = r_pred |> as.matrix()
              self$loss_raw = as.numeric(loss)
              keep = private$scored_terms(loss, nonfinite)
              # kept with its unit counts; weighted and summed after the time loop
              if (length(keep)) step_losses[[length(step_losses) + 1L]] = list(loss = loss, n = loss_n, keep = keep)
            }
          } else if(!is.na(period)) {
            # loss[1] = self$loss_dbh_func(y[, tmp_index,,1], Result[[1]][,(i-period+1):(i),]$mean(2) )
            # # ba
            # loss[2] = self$loss_ba_func(y[, tmp_index,,2], Result[[2]][,(i-period+1):(i),]$mean(2) )
            # # counts
            # loss[3] = self$loss_trees_func(y[,tmp_index,,3], Result[[3]][,(i-period+1):(i),]$mean(2))

            #browser()
            ## Growth aggregation over the period. Default: MEAN of the per-year
            ## rates, matched by an annualised per-year target. With
            ## `growth_period_scale = TRUE` the per-year rates are COMPOUNDED over
            ## the window (prod(1+g)-1) to give the cumulative PERIOD increment,
            ## matched by a de-annualised period target. This keeps env annual (so
            ## a per-year drought forcing still lands on the right year) while
            ## scoring growth at full period magnitude -- annualising the target
            ## instead collapses growth's variance and flattens its env response.
            if(!tree_mode) {
              g_pred = if (isTRUE(self$growth_period_scale))
                         (Result[[4]][,(i-period+1):(i),] + 1)$prod(2) - 1
                       else Result[[4]][,(i-period+1):(i),]$mean(2)
              loss[4] = self$loss_growth_func(y[,tmp_index,,4], g_pred, y[,tmp_index,,9])
              # mort rates
              loss[5] = self$loss_mortality_func(y[,tmp_index,,5], Result[[5]][,(i-period+1):(i),]$mean(2), y[,tmp_index,,8])
            }
            # reg rates ha
            loss[6] = self$loss_regeneration_func(y[,tmp_index,,6], Result[[7]][,(i-period+1):(i),]$sum(2) )#(Result[[7]][,(i-period+1),] - Result[[7]][,i,])$clamp(min = 0.0)  )
            self$obs_rec = y[,tmp_index,,6] |> as.matrix()
            self$pred_rec = Result[[7]][,(i-period+1):(i),]$sum(2) |> as.matrix()
            self$loss_raw = as.numeric(loss)
            # ---- Check loss before backward ----
            # (this branch used to backward() unconditionally, unlike the per-step
            # branch below; a non-finite period-averaged loss would otherwise inject
            # NaN gradients into every parameter with no diagnostic at all)
            keep = private$scored_terms(loss, nonfinite)
            if (length(keep) && isTRUE(loss$requires_grad)) {
              # (a loss with no graph - every response NA at this step - has
              # nothing to backpropagate; calling backward() on it is a hard error)
              l = private$loss_at_step(loss, keep, nonfinite)
              l$sum$backward()
              for(j in 1:7) Result[[j]] = Result[[j]]$detach()
              loss_total = loss_total + l$value
              loss_count = loss_count + 1L
            }
          } else {

            # loss[1] = self$loss_dbh_func(y[, tmp_index,,1], Result[[1]][,i,] )
            # # ba
            # loss[2] = self$loss_ba_func(y[, tmp_index,,2], Result[[2]][,i,] )
            # # counts
            # loss[3] = self$loss_trees_func(y[,tmp_index,,3], Result[[3]][,i,])

            if(!tree_mode) {
              loss[4] = self$loss_growth_func(y[,tmp_index,,4], Result[[4]][,i,], y[,tmp_index,,9])
              # mort rates
              loss[5] = self$loss_mortality_func(y[,tmp_index,,5], Result[[5]][,i,], y[,tmp_index,,8])
            }
            # reg rates ha
            loss[6] = self$loss_regeneration_func(y[,tmp_index,,6], Result[[7]][,i,])
            self$loss_raw = as.numeric(loss)
            # ---- Check loss before backward ----
            keep = private$scored_terms(loss, nonfinite)
            if (length(keep) && isTRUE(loss$requires_grad)) {
              l = private$loss_at_step(loss, keep, nonfinite)
              l$sum$backward()
              for(j in 1:7) Result[[j]] = Result[[j]]$detach()
              loss_total = loss_total + l$value
              loss_count = loss_count + 1L
            }
          }
          # Also check each component tensor you pass into distribution losses

        }

        # Truncate the recurrent cohort state every `update_step` years, independent of
        # whether this particular year has an observation to score a loss against. This
        # detach used to live inside the `i %in% year_sequence` block above, which meant
        # that for irregularly/sparsely observed data (e.g. multi-year FIA remeasurement
        # intervals that rarely line up with update_step) the state was effectively never
        # detached - i.e. full, un-truncated BPTT over the whole simulated horizon instead
        # of the truncated window `update_step` is supposed to bound. Detaching here does
        # not affect the gradient of losses already computed above: `Result` has already
        # recorded its dependency on the (still-attached) state at step i by this point;
        # this detach only prevents *future* steps from continuing to backprop through the
        # state as it existed at/before step i.
        anchor = if(isTRUE(reanchor)) anchor_cohorts[[as.character(i)]] else NULL
        if(is.null(anchor)) {
          dbh=dbh$detach()
          trees=trees$detach()
          species=species$detach()
          cohort_ids=cohort_ids$detach()
        } else {
          # Multiple shooting: the segment ends here, so the simulated state is REPLACED by
          # the observed cohorts rather than carried forward - a freshly built array carries
          # no graph, and ids continue from the counter so none is reused.
          base_id = next_cohort_id
          a = observed_state(anchor, base_id)
          next_cohort_id = base_id + a$species$shape[3]
          # a new segment is a new generation of tree identities; anything "merge" carries
          # across keeps its old id and stays out of the residual - the anchor never saw it
          if(tree_mode) {
            tree_gen_base = as.integer(base_id)
            tree_gen_n = as.integer(prod(a$species$shape))
            tree_gen_nc = as.integer(a$species$shape[3])
            tree_acc = tree_new(tree_gen_n)
            tree_live0 = tree_mark_live(a$trees, a$ids)
          }
          a_dbh = a$dbh
          a_trees = a$trees
          a_species = a$species
          a_ids = a$ids
          if(identical(anchor_mode, "merge")) {
            if(is.null(anchor_min_dbh)) stop("`anchor_mode = \"merge\"` requires `anchor_min_dbh`, the dbh below which the anchors are blind.")
            # The anchor is complete at and above `anchor_min_dbh`, so carrying the simulated
            # cohorts below it double counts nothing - replacing them instead deletes the whole
            # regeneration pool at every census. Zeroed trees is the dead-cohort convention, so
            # the above-threshold simulated slots are dropped by the next padding pass.
            keep = dbh$lt(anchor_min_dbh)$bitwise_and(trees$gt(0.5))$to(dtype = self$dtype)
            dbh = torch_cat(list(dbh$detach()*keep, a_dbh), 3)
            trees = torch_cat(list(trees$detach()*keep, a_trees), 3)
            species = torch_cat(list(species$detach(), a_species), 3)
            cohort_ids = torch_cat(list(cohort_ids$detach(), a_ids), 3)
          } else {
            dbh = a_dbh
            trees = a_trees
            species = a_species
            cohort_ids = a_ids
          }
        }
      }

      # this timestep's census closes the interval it observed, so everything after it is
      # conditioned on THIS observed stand instead. Outside `update_boundary` on purpose:
      # there is no trajectory here, so nothing depends on where the gradient is truncated.
      if(isTRUE(teacher_forcing) && !is.null(anchor_cohorts[[as.character(i)]])) {
        base_id = next_cohort_id
        forced = observed_state(anchor_cohorts[[as.character(i)]], base_id)
        next_cohort_id = base_id + forced$species$shape[3]
        if(tree_mode) {
          tree_gen_base = as.integer(base_id)
          tree_gen_n = as.integer(prod(forced$species$shape))
          tree_gen_nc = as.integer(forced$species$shape[3])
          tree_acc = tree_new(tree_gen_n)
          tree_live0 = tree_mark_live(forced$trees, forced$ids)
        }
      }

      if(!is.null(y)) {
        #for(j in 1:3) Result[[j]] = Result[[j]]$detach()
        # if period_length = NA on all sites, it is assumed that we need only yearly gradients
        if(as.numeric(y[,,,7]$isnan()$bitwise_not()$sum()) < 0.5) for(j in 1:7) Result[[j]] = Result[[j]]$detach()

      }
      if(!per_site_mode) loss$detach_()

      if(verbose) cli::cli_progress_update()

    }

    # browser()

    # per-site windows: each term becomes one mean over all units scored in this
    # call - step s enters with weight n_s / sum(n), the units of the steps where
    # the term was kept - and one backward runs over the result. A term no step
    # scored a unit of (no observations; its losses are 0) keeps weight 1.
    if(per_site_mode && length(step_losses)) {
      n_call = numeric(7L)
      for(s in step_losses) n_call[s$keep] = n_call[s$keep] + s$n[s$keep]
      for(s in step_losses) {
        w = ifelse(n_call > 0, s$n / n_call, 1)
        l = private$loss_at_step(s$loss * torch_tensor(w, dtype = s$loss$dtype, device = self$device),
                                 s$keep, nonfinite)
        loss_deferred = loss_deferred + l$sum
        loss_total = loss_total + l$value
      }
      # loss_total already is the pooled mean; reported as is
      loss_count = 1L
    }
    if(per_site_mode && isTRUE(loss_deferred$requires_grad)) {
      loss_deferred$backward()
      for(j in 1:7) Result[[j]] = Result[[j]]$detach()
    }

    # report the mean loss over every update_step boundary that was actually scored
    # during this forward() call, rather than whichever timestep happened to be last
    # in the loop - the latter is an arbitrary, often misleading slice of one call's
    # worth of observations and made training-progress plots hard to interpret.
    loss_out = if (loss_count > 0L) loss_total / loss_count else loss_total

    names(Result) =  c("dbh","ba", "trees", "growth", "mort", "reg", "r_mean_ha")
    Predictions = list(Site = lapply(Result, function(x) torch::as_array(x)))
    # `debug` additionally exposes the per-timestep patch state (unchanged)
    if(debug) Predictions$Patch = Raw_patch_results
    # cohorts are attached whenever any timestep was recorded (debug or the
    # user-facing `return_cohorts`); pred2DF() turns them into a tidy table
    if(record_cohorts) {
      Predictions$Cohort = list(
        cohortStates = Raw_cohort_results,
        cohortID = Raw_cohort_ids
      )
    }
    Result_out = list(Predictions = Predictions, loss = loss_out)
    if(length(tree_cov)) Result_out$tree_deaths = do.call(rbind, tree_cov)

    # (parameters' requires_grad is restored via on.exit() above, even if this
    # function exits early due to an error)

    return(Result_out)
  },

  #' @description
  #' Predict growth and mortality for INDIVIDUAL observed trees, the way
  #' `loss_aggregation = "tree"` scores them during fitting.
  #'
  #' The observed stand replaces the simulated one: light is computed from the
  #' observed trees at their observed diameters, the growth and mortality
  #' kernels are evaluated on those trees, and the annual rates are accumulated
  #' over each site's remeasurement interval (growth as the mean, mortality as
  #' the complement of the product of annual survivals). This is what validating
  #' a tree-level fit needs, and it is deliberately the same code path the loss
  #' uses, so a report cannot drift from the objective.
  #'
  #' @param tree_data (`data.table|data.frame`)\cr Observed trees, as for
  #'   [FINN::fit]'s `tree_data`: `siteID`, `year`, `patchID`, `species`, `dbh`
  #'   (at the interval start) and `period_length` (the interval, in timesteps).
  #' @param env (`data.table|data.frame`)\cr Environment, as for [FINN::fit].
  #' @param patches (`integer(1)`)\cr Not used: the patch dimension is the
  #'   number of distinct `patchID`s in `tree_data`. Kept for compatibility.
  #' @param patch_size (`numeric(1)`)\cr Patch size. Defaults to the one the model
  #'   was fitted with (`$patch_size_ha`); pass a value only to override it.
  #' @param device (`character(1)`)\cr `"cpu"` or `"gpu"`; `"cuda"` is accepted as a
#'   synonym for `"gpu"`.
  #' @return One row per row of `tree_data`, IN THE SAME ORDER (`.row` records
  #'   the original position), with `growth_pred` (mean annual relative growth
  #'   over the interval) and `mort_pred` (probability of dying during the
  #'   interval). Bind these onto `tree_data` directly; joining on
  #'   siteID/year/patchID would be many-to-many, since a patch holds many trees
  #'   in the same year.
  predict_trees = function(tree_data, env, patches = 100L, patch_size = NULL, device = c("cpu", "gpu")) {
    device = resolve_device(device)
    self$device = if (identical(device, "gpu")) "cuda" else "cpu"
    # light is f(BA / patch_size_ha), so scoring at another patch size than the fit
    # silently rescales competition - the one thing this route exists to rule out
    if (is.null(patch_size)) patch_size = self$patch_size_ha
    self$patch_size_ha = patch_size
    d = data.table::as.data.table(tree_data)
    if (!"period_length" %in% colnames(d))
      stop("`tree_data` needs a `period_length` column: the rates are interval quantities.", call. = FALSE)
    obs_years = sort(unique(d$year)); site_ids = sort(unique(d$siteID))
    tobs = private$build_tree_obs(d, site_ids, obs_years)
    envs = private$extract_env_method(env)
    # the environment networks are built lazily, on the first fit() or predict();
    # build them here too so this works on a model that has never been run
    private$create_nn(self$process_growth, "growth", dim(envs$growth_env)[3])
    private$create_nn(self$process_mortality, "mortality", dim(envs$mortality_env)[3])
    envt = list(growth = torch::torch_tensor(envs$growth_env, dtype = self$dtype, device = self$device),
                mort   = torch::torch_tensor(envs$mortality_env, dtype = self$dtype, device = self$device))
    year_sequence = which(levels(as.factor(env$year)) %in% levels(as.factor(d$year)))
    out = data.table::copy(tobs$rows)[, `:=`(growth_pred = NA_real_, mort_pred = NA_real_)]
    torch::with_no_grad({
      for (k in seq_along(obs_years)) {
        i = year_sequence[k]
        pl = merge(data.table::data.table(siteID = site_ids),
                   unique(d[year == obs_years[k], .(siteID, period_length)]), by = "siteID", all.x = TRUE)
        p_t = torch::torch_tensor(ifelse(is.na(pl$period_length), 0, pl$period_length),
                                  dtype = self$dtype, device = self$device)
        tt = torch::torch_arange(1, i, dtype = self$dtype, device = self$device)
        M  = tt$unsqueeze(1)$gt((i - p_t)$unsqueeze(2))$to(dtype = self$dtype) *
             p_t$gt(0.5)$to(dtype = self$dtype)$unsqueeze(2)
        tp = private$tree_predictions(tobs, k, envt, i, M, p_t)
        g  = as.array(tp$growth$cpu()); m = as.array(tp$mortality$cpu())
        sel = out$y == k
        ix  = cbind(out$s[sel], out$p[sel], out$slot[sel])
        out[sel, `:=`(growth_pred = g[ix], mort_pred = m[ix])]
      }
    })
    # back in the caller's own row order, so the predictions can simply be
    # cbind()ed onto `tree_data` (a join on siteID/year/patchID would be
    # many-to-many: a patch holds many trees in the same year)
    data.table::setorder(out, .row)
    out[, c("s", "y", "p", "slot") := NULL][]
  },

  simulate = function(env,
                      disturbance = NULL,
                      patches = 100L,
                      patch_size = 0.1,
                      init_cohort = NULL,
                      batchsize = NULL,
                      device = c("cpu", "gpu"),
                      return_cohorts = FALSE,
                      debug = FALSE
                      ) {
    device = resolve_device(device)
    if(device == "gpu") device="cuda:0"
    self$device = device
    self$patch_size_ha = patch_size
    if (is.null(self$recruit_obs_weight)) self$recruit_obs_weight = 1.0  # inventory-simulator: observation-operator inclusion weight (1 = per-ha, unchanged)
    if (is.null(self$growth_period_scale)) self$growth_period_scale = FALSE  # FALSE = mean of annual rates (default); TRUE = compounded period increment
    # `env` and `init_cohort` are paired by site position; the sorted user ids map
    # the positions back, so the output carries the caller's own siteIDs.
    site_ids_user = check_site_sets(env, init_cohort = init_cohort, caller = "predict()")
    envs = private$extract_env_method(env)

    if(!is.null(disturbance)) {
      disturbance = extract_env(~0+intensity, disturbance)
    }

    # init networks if not yet done
    private$create_nn(self$process_mortality, "mortality", dim(envs$mortality_env)[3])
    private$create_nn(self$process_growth, "growth", dim(envs$growth_env)[3])
    private$create_nn(self$process_regeneration, "regeneration", dim(envs$regeneration_env)[3])

    if(is.null(batchsize)) batchsize = nrow(envs[[1]])
    sites =  nrow(envs[[1]])

    sp <- self$N_species
    local_init <- if (is.null(init_cohort)) {
      CohortMat(
        dims  = c(sites, patches, sp),
        dbh   = array(1, dim = c(sites, patches, sp)),
        trees = array(1, dim = c(sites, patches, sp)),
        sp    = sp
      )
    } else {
      init_cohort
    }
    self$to(device = device)
    if(is.null(disturbance)) {
      data = tensor_dataset(torch_tensor(envs[[1]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[2]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[3]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_arange(1,dim(envs[[1]])[1])$to(dtype = torch_int64() , device=torch_device('cpu'))
      )
    } else {
      data = tensor_dataset(torch_tensor(envs[[1]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[2]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[3]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_arange(1, dim(envs[[1]])[1] )$to(dtype = torch_int64() , device=torch_device('cpu')),
                            torch_tensor(disturbance, dtype=self$dtype, device=torch_device('cpu'))
      )
    }

    self$eval()

    DataLoader = torch::dataloader(data, batch_size=batchsize, shuffle=FALSE, num_workers=0, pin_memory=TRUE, drop_last=FALSE)

    predictions_batch = list()
    coro::loop(for (b in DataLoader) {
      x_mort = b[[1]]$to(device = self$device, non_blocking=TRUE)
      x_growth = b[[2]]$to(device = self$device, non_blocking=TRUE)
      x_reg = b[[3]]$to(device = self$device, non_blocking=TRUE)
      ind = b[[4]]$to(device = self$device, non_blocking=TRUE)
      dist = NULL
      if(!is.null(disturbance)) dist = b[[5]]$to(device = self$device, non_blocking=TRUE)

      trees   = local_init$trees$to(device = self$device, non_blocking = TRUE)[ind,]
      species = local_init$species$to(device = self$device, non_blocking = TRUE)[ind,]
      dbh     = local_init$dbh$to(device = self$device, non_blocking = TRUE)[ind,]
      pred_tmp = self$forward(dbh = dbh,
                              trees = trees,
                              species = species,
                              env = list(x_mort, x_growth, x_reg),
                              disturbance = dist,
                              verbose = TRUE,
                              return_cohorts = return_cohorts,
                              debug = debug)
      pred = pred_tmp[[1]]
      long_b = pred2DF(list(Predictions = pred), "long")
      wide_b = pred2DF(list(Predictions = pred), "wide")
      # pred2DF() labels sites 1..batch_size within this batch; remap those to the
      # batch's global site positions (`ind`), and those to the caller's own
      # siteIDs (`site_ids_user`, sorted as the env array). For ids that are 1:N
      # in a single batch both steps are the identity.
      ind_r = as.integer(torch::as_array(ind$cpu()))
      for (nm in intersect(names(long_b), c("site", "patch", "cohort"))) {
        long_b[[nm]][, siteID := site_ids_user[ind_r[siteID]]]
        wide_b[[nm]][, siteID := site_ids_user[ind_r[siteID]]]
      }
      predictions_batch = append(predictions_batch, list(list(long = long_b, wide = wide_b)))

      })
    predictions =
      list(long = list(site = rbindlist( lapply(predictions_batch, function(X) X$long$site ))),
           wide = list(site = rbindlist( lapply(predictions_batch, function(X) X$wide$site )))
      )
    # If cohorts were requested, forward() will have stored them and pred2DF()
    # will have produced a `$cohort` table per batch; surface them alongside
    # the site output in the same nested structure.
    if(!is.null(predictions_batch[[1]]$long$cohort)) {
      predictions$long$cohort = rbindlist( lapply(predictions_batch, function(X) X$long$cohort) )
      predictions$wide$cohort = rbindlist( lapply(predictions_batch, function(X) X$wide$cohort) )
    }
    if(debug){
      predictions = list(wide = pred2DF(pred_tmp, format = "wide"), long = pred2DF(pred_tmp, format = "long"))
    }
    return(predictions)

  },


  fit = function(data = NULL,
                 env,
                 disturbance = NULL,
                 patches = 100L,
                 patch_size = 0.1,
                 init_cohort = NULL,
                 epochs = 20L,
                 lr = 0.01,
                 # learning-rate schedule applied on top of the constant `lr` above.
                 # "none" (default) reproduces the previous constant-lr behavior exactly.
                 # Built-ins: "step" (lr_step), "exponential" (lr_multiplicative),
                 # "cosine" (lr_cosine_annealing), "plateau" (lr_reduce_on_plateau, driven
                 # by the epoch's total loss). Advanced: pass a function(optimizer) that
                 # returns a torch lr_scheduler object (must implement $step()).
                 lr_scheduler = "none",
                 # tuning overrides for the scheduler chosen above, e.g.
                 # list(step_size = 50, gamma = 0.5) for "step", list(T_max = 200) for
                 # "cosine", list(factor = 0.5, patience = 20) for "plateau". Unset
                 # entries fall back to defaults scaled off `epochs`/`lr`.
                 lr_scheduler_params = list(),
                 # # 1 -> dbh, 2 -> ba, 3 -> trees, 4 -> growth rates, 5 -> mort rates, 6 -> reg rates
                 # `loss` is the old name for `loss_family`; kept as a deprecated alias.
                 loss = NULL,
                 loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
                 # How the model's per-cohort quantities are reduced to the observed
                 # quantity, one entry per response. Character when no parameters are
                 # needed, a list when they are. See the argument docs.
                 loss_aggregation = NULL,
                 # Per-tree observations, one row per tree and census: siteID, year,
                 # patchID, species, growth, died, plus dbh for aggregation "tree"
                 # (growth / mortality) and cohortID for loss_resolution = "tree".
                 tree_data = NULL,
                 weights = "auto",
                 weight_scale = c("auto", "raw"),
                   optimizer = optim_ignite_adam,
                 batchsize = NULL,
                 device = c("cpu", "gpu"),
                 update_step = 1L,
                 start_time = 1L,
                 plot_progress = TRUE,
                 folder = NULL,
                 checkpoints = 100L,
                 shuffle = TRUE,
                 record_gradients = FALSE,
                 env_autoscale = TRUE,
                 # gradient-norm budget passed to nn_utils_clip_grad_norm_(), applied
                 # separately per parameter group (see grouping below) rather than once
                 # globally. A single number applies the same budget to every group;
                 # a named list/vector (any of "mechanistic", "nn", "loss") overrides
                 # individual groups, e.g. clip_norm = list(loss = 5, nn = 1).
                 clip_norm = 2.0,
                 # a non-finite loss term at a timestep: "drop_term" backpropagates the
                 # other terms of that timestep, "skip_step" discards the whole timestep
                 # (the behaviour before drop_term existed).
                 nonfinite = c("drop_term", "skip_step"),
                 # Keep the drawn count in the forward pass and take the derivative
                 # through its expectation (FINNverse/FINN#53); FALSE reproduces the
                 # previous behaviour exactly.
                 rb_gradient = TRUE,
                 # Offset keeping the count likelihoods finite at a simulated zero.
                 trees_offset = 0.01,
                 # multiple shooting: restart the simulated cohort state from the observed
                 # cohorts in `anchor_cohorts` (named by simulation timestep) instead of
                 # running free from t = 0. FALSE is single shooting, i.e. unchanged.
                 reanchor = FALSE,
                 anchor_cohorts = NULL,
                 # "replace" overwrites the simulated state with the anchor; "merge" keeps the
                 # simulated cohorts below `anchor_min_dbh`, which the anchor cannot observe.
                 anchor_mode = "replace",
                 anchor_min_dbh = NULL,
                 # direct estimation: reset the cohort state to the last observed stand at
                 # EVERY timestep, so the rates are evaluated at an observed state and no
                 # trajectory is simulated. Drops the dbh/ba/trees loss terms.
                 teacher_forcing = FALSE,
                 # the resolution growth and mortality are scored at: per (site, year, species)
                 # cell, or per individual tree. "site" is the validated, bit-for-bit path.
                 loss_resolution = c("site", "tree"),
                 ...) {

    # Only touch the global graphics state when we actually draw the training
    # curve. Capturing par(no.readonly = TRUE) and restoring it otherwise is a
    # no-op that can still emit a spurious "par(new=TRUE) with no plot" warning
    # (when a previous plot left par$new = TRUE and no device has an active
    # plot). Restore via on.exit() so an early error can't leak the caller's par.
    if (isTRUE(plot_progress)) {
      old_par = par(no.readonly = TRUE)
      on.exit(do.call(par, old_par), add = TRUE)
    }
    # (The families are validated per response in create_loss_functions(), which
    # knows the full supported set and which families a given response allows;
    # the blanket any()-check that used to sit here predates `binomial`/`none`.)

    # Internal z-standardization of environmental predictors. Learn the
    # centre/scale from the (raw) training env now and store them on the model;
    # extract_env_method() re-applies them here and at every predict/simulate
    # call, so raw env can be supplied throughout. See compute_env_scaling().
    # env_autoscale may be TRUE/FALSE, or a per-predictor spec (strings
    # "auto"/"identity"/"0to1" or functions) -- see compute_env_scaling(). Any
    # non-FALSE value computes a scaling table; FALSE leaves env untouched.
    # Backward-compatible dispatch: 0.1.0 scaled only when env_autoscale was TRUE
    # (isTRUE), so FALSE / NULL / NA all skipped scaling. Preserve that exactly,
    # and additionally accept the 0.2.0 per-predictor spec forms (a character
    # vector, a list, or a function).
    scale_env <- isTRUE(env_autoscale) || is.character(env_autoscale) ||
                 is.list(env_autoscale) || is.function(env_autoscale)
    if (scale_env) self$env_scaling = compute_env_scaling(env, env_autoscale)

    if(is.null(data)) {
      print("No data. Switching into simulation modus...")
      return(self$predict(env = env, patches = patches, patch_size = patch_size, batchsize=batchsize, device=device))
    }
    # `data`, `env` and `init_cohort` are paired by site POSITION, so they must
    # hold the same site set (the ids need not be 1:N). The sorted ids are kept
    # as the lookup table from array position to user id.
    site_ids_user = check_site_sets(env, data = data, init_cohort = init_cohort, caller = "fit()")
    self$site_lookup = data.table::data.table(siteID = site_ids_user, index = seq_along(site_ids_user))

    self$train_env         = data.table::copy(data.table::as.data.table(env))
    self$train_init_cohort = if(is.null(init_cohort)) NULL else list(cohort = init_cohort)
    self$fit_id            = as.numeric(Sys.time())
    self$conditional_effects = NULL
    self$ale                 = NULL
    self$perm_importance     = NULL
    # The observation tensors of the non-scalar aggregations are rebuilt below
    # for exactly the aggregations of THIS call. Left over from an earlier fit()
    # (or a torch_load()ed model) they would be stale, or dead tensors, and
    # forward() would still evaluate them (tree_obs runs tree_predictions()
    # whenever it is non-NULL).
    self$tree_obs  = NULL
    self$class_obs = NULL
    self$quant_obs = NULL

    # setup loss functions
    # weights = "auto" (the default): epoch 1 runs unweighted so each loss term
    # can be MEASURED, then the weights are set to the inverse of those
    # magnitudes (see the rescale block in the epoch loop). A numeric vector is
    # used as-is and disables the rescaling.
    # `loss` (the argument) names the loss FAMILY per response. The batch loop
    # below reuses the name `loss` for a torch tensor of the current losses,
    # which shadows this argument - so snapshot the spec now, and use the
    # snapshot anywhere the families are needed later (e.g. the "auto" rescale).
    if(!is.null(loss)) {
      warning("`loss` is deprecated; use `loss_family`. Using the value of `loss` for now.", call. = FALSE)
      loss_family = loss
    }
    loss_spec = loss_family
    nonfinite = match.arg(nonfinite)
    self$nonfinite = nonfinite
    # `loss_resolution` moves the growth and mortality terms between resolutions, it does not
    # add terms: the spec, the weights and the history keep their six positions in both modes
    loss_resolution = match.arg(loss_resolution)
    tree_res = identical(loss_resolution, "tree")
    if (tree_res && is.null(tree_data))
      stop("`loss_resolution = \"tree\"` requires `tree_data`: the observed trees, one row per tree and census, with a `cohortID` column naming the initial cohort each tree started as.", call. = FALSE)
    if (tree_res && is.null(init_cohort))
      stop("`loss_resolution = \"tree\"` needs `init_cohort` with one cohort per observed tree (makeInitCohorts() without `dbh_binsize`); the trees are scored on the cohorts they started as.", call. = FALSE)
    # `tree_data` is a TABLE for the user in both per-tree options; the arrays
    # loss_resolution = "tree" works on are built from it below (build_tree_targets).
    if (!is.null(tree_data) && !is.data.frame(tree_data))
      stop("`tree_data` must be a data.frame / data.table with one row per observed tree and census (columns siteID, year, patchID, species, growth, died, plus dbh or cohortID depending on the option); see ?fit, section 'Scoring individual trees'.", call. = FALSE)
    # Resolve how each response is reduced from cohorts to the observed quantity.
    # Stored on the model because it fixes tensor widths used inside forward().
    self$loss_agg = private$resolve_loss_aggregation(loss_aggregation, names(loss_spec))
    # Two ways to score single trees exist side by side and ask different
    # questions (observed vs simulated tree, see ?fit); the same response cannot
    # be scored both ways in one fit.
    # Retention keeps a killed tracked cohort in the array with trees = 0. Every stand
    # reduction weights by stem number and drops it -- except the height quantile of the
    # BINNED competition, which counts slots. Measured on a one-patch stand, the same live
    # cohorts with and without ghosts beside them: without the mask a live cohort's light
    # moves by up to 1.00 (fully shaded to open); with it, by 0 exactly, in every cell of a
    # 3 x 5 x 20 grid. So the tree term is only inert on the trajectory under a masked or
    # unbinned light.
    if (tree_res && !isTRUE(self$continuous) && !isTRUE(self$light_mask))
      warning("loss_resolution = \"tree\" keeps a cohort the simulation kills in the array until its interval ends, and the binned competition takes its height quantiles over every slot, so those cohorts move the light of the living ones -- by up to 1.00 in a one-patch test. Use createProcess(light_mask = TRUE) or continuous = TRUE for the competition process, or the per-tree term changes the trajectory it is scoring.",
              call. = FALSE)
    if (tree_res && any(vapply(self$loss_agg[c("growth", "mortality")], function(a) identical(a$type, "tree"), logical(1))))
      stop("`loss_resolution = \"tree\"` (the simulated cohort each tree started as) cannot be combined with `loss_aggregation = \"tree\"` for growth or mortality (the observed tree); choose one of the two ways of scoring individual trees.", call. = FALSE)
    weight_scale = match.arg(weight_scale, c("auto", "raw"))
    auto_weights = identical(weight_scale, "auto")
    if (is.character(weights) && length(weights) == 1L && identical(weights, "auto")) weights = rep(1, 6)
    if (!is.numeric(weights) || length(weights) != 6L || !all(is.finite(weights)) || any(weights < 0))
      stop("`weights` must be \"auto\" or a numeric vector of length 6 (dbh, ba, trees, growth, mortality, regeneration), non-negative and finite.",
           call. = FALSE)
    user_weights = weights
    if (!auto_weights) self$loss_weights = weights
    # NOTE: the loss closures are built further down, AFTER the response tensor Y
    # exists - "auto" needs the observed responses to compute the intercept-only
    # baseline it scales by.

    # setup data
    device = resolve_device(device)
    if(device == "gpu") device="cuda:0"
    self$device = device
    self$patch_size_ha = patch_size
    if (is.null(self$recruit_obs_weight)) self$recruit_obs_weight = 1.0  # inventory-simulator: observation-operator inclusion weight (1 = per-ha, unchanged)
    if (is.null(self$growth_period_scale)) self$growth_period_scale = FALSE  # FALSE = mean of annual rates (default); TRUE = compounded period increment

    # model.matrix() in extract_env() must keep NA rows so the response arrays
    # stay dimensionally aligned. Set na.pass, and restore the user's prior
    # na.action via on.exit() (immediately, per CRAN policy) so an error before
    # the reset below cannot leak the change into the user's session.
    oldopts <- options(na.action = 'na.pass')
    on.exit(options(oldopts), add = TRUE)
    sp = self$N_species
    if(!"period_length" %in% colnames(data)) data$period_length = NA_real_
    # Number of trees at risk behind each observed mortality proportion. Weighting
    # the binomial likelihood by it is what makes a 100-tree observation count for
    # more than a 1-tree one. Data built before makeObsData emitted this column
    # falls back to 1, which reduces the binomial exactly to an unweighted one.
    if(!"n_at_risk" %in% colnames(data)) data$n_at_risk = 1.0
    if(!"growth_n" %in% colnames(data)) data$growth_n = 1.0
    response = list(
      dbh = abind::abind(lapply(1:self$N_species, function(i) extract_env(~0+dbh, data[data$species==i,])), along = 3L),
      ba = abind::abind(lapply(1:sp, function(i) extract_env(~0+ba, data[data$species==i,])), along = 3L),
      trees = abind::abind(lapply(1:sp, function(i) extract_env(~0+trees, data[data$species==i,])), along = 3L),
      growth = abind::abind(lapply(1:sp, function(i) extract_env(~0+growth, data[data$species==i,])), along = 3L),
      mort = abind::abind(lapply(1:sp, function(i) extract_env(~0+mort, data[data$species==i,])), along = 3L),
      reg = abind::abind(lapply(1:sp, function(i) extract_env(~0+reg, data[data$species==i,])), along = 3L),
      period_length = abind::abind(lapply(1:sp, function(i) extract_env(~0+period_length, data[data$species==i,])), along = 3L),
      # slice 8 - appended last so the existing y[,,,1:7] indices stay valid
      n_at_risk = abind::abind(lapply(1:sp, function(i) extract_env(~0+n_at_risk, data[data$species==i,])), along = 3L),
      # slice 9 - trees behind each growth mean. 34% of FIA growth observations
      # rest on a SINGLE tree while others average 32, and the variance of a mean
      # is sigma^2/n - so the Gaussian term is weighted by n, exactly as the
      # binomial term is weighted by n_at_risk.
      growth_n = abind::abind(lapply(1:sp, function(i) extract_env(~0+growth_n, data[data$species==i,])), along = 3L)
    )
    # Restore the user's na.action now that the NA-preserving response arrays are
    # built, so the rest of fit() behaves exactly as in the submitted version
    # (the env extraction below runs under the user's na.action, not na.pass).
    # The on.exit() set above is the error-safety net for the block in between.
    options(oldopts)

    Y = torch::torch_cat(lapply(response, function(x) torch::torch_tensor(x, dtype=torch::torch_float32(), device="cpu")$unsqueeze(4)), 4)
    # The recruitment cap's prior needs the data, which finn() never sees. `reg`
    # is the SUM over an observation's interval while the cap bounds the mean of
    # ONE timestep, so the observations are annualised before the maximum is
    # taken - on a 10-year inventory the prior would otherwise sit 10x above
    # anything the cap can bind.
    private$resolve_cap_prior(if(is.null(data$reg)) numeric(0) else {
      pl_obs = as.numeric(data$period_length)
      pl_obs[!is.finite(pl_obs) | pl_obs <= 0] = 1
      as.numeric(data$reg) / pl_obs
    })
    # ---- constant vs per-site aggregation windows ----
    # The original path takes ONE period_length per data-year across all sites.
    # If sites observed in the same data-year carry different period_length
    # values, or if some sites have no observation in a year others do (NA
    # period next to a finite one), each site needs its own window: switch
    # forward() to the per-site path. `self$period_mode` ("auto" / "constant" /
    # "per_site") overrides the detection.
    pl = array(response$period_length[,,1], dim = dim(response$period_length)[1:2])
    agg_types = vapply(self$loss_agg, function(a) a$type, character(1))
    varies = FALSE
    for(k in seq_len(ncol(pl))) {
      v = pl[,k]; nn = v[!is.na(v)]
      if(length(nn) && (anyNA(v) || length(unique(nn)) > 1L)) { varies = TRUE; break }
    }
    # The per-tree terms are only evaluated on the per-site path, while
    # compute_baseline_losses() measures their null on the per-tree observation
    # regardless. On a constant period that combination does not fall back to
    # the cell term, it scores the cell term under a weight taken from the
    # per-tree null, so there is no correct constant path and an explicit
    # period_mode = "constant" is a contradiction rather than a preference.
    tree_agg = any(agg_types == "tree")
    if(tree_agg && identical(self$period_mode, "constant"))
      stop("period_mode = \"constant\" cannot be combined with loss_aggregation = \"tree\": the per-tree terms are only evaluated on the per-site path, and on the constant path the cell term would be scored under a weight taken from the per-tree null. Use period_mode = \"auto\" or \"per_site\", or drop the \"tree\" aggregation.", call. = FALSE)
    per_site = varies || tree_agg || identical(self$period_mode, "per_site")
    if(identical(self$period_mode, "constant")) per_site = FALSE
    self$per_site_period = per_site
    if(per_site) cli::cli_alert_info(paste0(
      "Using per-site aggregation windows (one backward per batch): ",
      if(varies) "period_length varies between sites"
      else if(tree_agg) "loss_aggregation = \"tree\" needs them"
      else "period_mode = \"per_site\""))

    # Observation tensors for the non-scalar aggregations. Built once, kept whole
    # on the model, and sliced by the batch index in the epoch loop (like
    # init_cohort), so the batching stays exactly as it was. This has to happen
    # before the baseline weights below, which need these observations for the
    # responses that no longer live in `Y`.
    obs_years = sort(unique(data$year))
    site_ids  = sort(unique(data$siteID))
    if (any(agg_types == "size_classes")) {
      r = names(agg_types)[agg_types == "size_classes"][1]
      self$class_breaks = self$loss_agg[[r]]$args$breaks
      self$class_tau    = if (!is.null(self$loss_agg[[r]]$args$tau)) self$loss_agg[[r]]$args$tau else 1.0
      self$class_obs    = private$build_class_obs(data, K = length(self$class_breaks) + 1L)
      message(sprintf("[loss] '%s' aggregated into %d diameter classes at %s (tau = %g)",
                      r, length(self$class_breaks) + 1L, paste(self$class_breaks, collapse = ", "), self$class_tau))
    }
    if (any(agg_types == "quantiles")) {
      r = names(agg_types)[agg_types == "quantiles"][1]
      self$quant_probs = self$loss_agg[[r]]$args$probs
      self$quant_tau   = if (!is.null(self$loss_agg[[r]]$args$tau)) self$loss_agg[[r]]$args$tau else 0.05
      self$quant_obs   = private$build_quantile_obs(data, self$quant_probs, r)
      message(sprintf("[loss] '%s' aggregated into the %s quantiles (tau = %g)",
                      r, paste(self$quant_probs, collapse = ", "), self$quant_tau))
    }
    if (any(agg_types == "tree")) {
      if (is.null(tree_data)) stop("loss_aggregation 'tree' needs `tree_data`.", call. = FALSE)
      self$tree_obs = private$build_tree_obs(tree_data, site_ids, obs_years)
    }
    # loss_resolution = "tree": the table becomes one target array per scored
    # census, aligned to the cohort slots of the layout that census's interval
    # starts from (init_cohort, or the previous anchor). Built here, before the
    # baselines, because the per-tree null model is measured on these arrays.
    tree_targets = NULL
    if (tree_res) {
      anchor_steps = if ((isTRUE(reanchor) || isTRUE(teacher_forcing)) && !is.null(anchor_cohorts))
        as.integer(names(anchor_cohorts)) else integer(0)
      tree_targets = private$build_tree_targets(
        tree_data, site_ids, obs_years,
        year_sequence = which(levels(as.factor(env$year)) %in% levels(as.factor(data$year))),
        init_cohort = init_cohort, anchor_cohorts = anchor_cohorts, anchor_steps = anchor_steps)
    }

    # ---- weights = "auto": scale each loss by its intercept-only baseline ----
    # The six losses are summed, and their raw magnitudes differ by orders of
    # magnitude (on FIA: dbh is an MSE in cm^2 ~ 1e2; growth an MSE on a ratio
    # ~ 1e-2). Summed as-is, dbh takes ~87% of the objective and growth ~1%, so
    # growth never learns - and no fixed weight vector fixes this in general,
    # because the imbalance follows the units of the user's data.
    #
    # The reference is the null model: what the loss would be if we predicted the
    # best CONSTANT (the intercept). Weighting by 1/baseline makes every term a
    # fraction of its own null deviance, so they are commensurable and each one
    # reads directly as "how much better than the mean": 1 = no better, < 1 =
    # better. For an MSE this is exactly Max's 1/sd(y)^2 - the loss of predicting
    # the mean IS the variance - but it generalises to the Poisson, negative
    # binomial and binomial terms, where sd() is not the right scale.
    # set before the loss closures are built: they capture `self$trees_offset`
    # A draw without a pathwise gradient, with the expectation gradient switched off,
    # leaves the mortality parameters with EXACTLY zero gradient -- measured -- and
    # nothing warns, because the other loss terms keep falling.
    if(!isTRUE(rb_gradient) && identical(self$mortality_draw, "poisson"))
      stop('mortality_draw = "poisson" has no pathwise gradient, so with rb_gradient = FALSE the ',
           'mortality parameters receive exactly zero gradient. Use rb_gradient = TRUE, or ',
           'mortality_draw = "gamma" (exact moments and the right gradient) or "bernoulli".',
           call. = FALSE)
    self$rb_gradient  = isTRUE(rb_gradient)
    self$trees_offset = trees_offset
    private$create_loss_functions(loss_spec, weights)   # unweighted for now
    # computed in both modes: under weight_scale = "raw" it is not used for scaling,
    # only to say what each response actually contributes (FINNverse/FINN#37)
    base = if (auto_weights)
      private$compute_baseline_losses(Y, loss_spec, self$class_obs, self$tree_obs, self$quant_obs)
    else
      # only used to report the shares, so a baseline we cannot compute is not worth
      # a warning here - the scaled path still warns, because there it changes the fit
      suppressWarnings(private$compute_baseline_losses(Y, loss_spec, self$class_obs, self$tree_obs, self$quant_obs))
    self$loss_baseline = base
    # A state term on trees or ba cannot train regeneration: recruits reach dbh through
    # a sigmoid of scale 1e-3, whose derivative is zero almost everywhere, and the stem
    # path runs through a rounding whose forward value is 0 below 0.5 recruits per patch
    # and species. Measured on BWI, a stems term left recruitment at 3.08x the target and
    # kept the recruitment loss from moving (FINNverse/FINN#43). The loss can lack an
    # observed recruitment in two ways: the term is switched off, or the family is set but
    # the scored data carry no finite value for it.
    if (isTRUE(self$par_regeneration_optimized) || isTRUE(self$env_regeneration_optimized)) {
      reg_off = identical(as.character(loss_spec[["regeneration"]]), "none")
      reg_empty = !reg_off && torch::with_no_grad({
        m = Y[,,,6]$isnan()$bitwise_not()
        !as.logical(m$any()$cpu())
      })
      if (reg_off || reg_empty)
        warning(sprintf("the regeneration process is being calibrated, but %s, so no observed recruitment enters the loss. No other response can stand in for it: recruits reach `dbh` through a near-step sigmoid and `trees` through a rounding that is 0 below 0.5 recruits per patch and species. Score regeneration against an observation, or set optimizeSpecies = FALSE and optimizeEnv = FALSE on the regeneration process. See FINNverse/FINN#43.",
                        if (reg_off) "`loss_family['regeneration']` is \"none\"" else "`data` carries no finite `reg` value in the scored cells"),
                call. = FALSE)
    }
    if (auto_weights) {
      # at tree resolution those two terms are scored over trees, so the null model they are
      # measured against is the intercept-only fit to the per-tree targets, not to the cells
      if (tree_res) base[c("growth", "mortality")] = sapply(c(growth = "growth", mortality = "mort"), function(r) {
        v = unlist(lapply(tree_targets, function(td) as.numeric(td[[r]])))
        v = v[is.finite(v)]
        # a target with one class only has no entropy, so 1/baseline explodes: leave it unscaled
        if (!length(v) || (r == "mort" && length(unique(v)) < 2L)) return(NA_real_)
        tt = torch_tensor(v, dtype = self$dtype, device = "cpu")
        f = if (r == "growth") self$loss_growth_func else self$loss_mortality_func
        as.numeric(f(tt, torch::torch_full_like(tt, mean(v))))
      })
      w = 1/base
      w[!is.finite(w) | w <= 0] = 1        # no data / degenerate -> leave alone
      # the user's numbers multiply the scaled objective, they do not replace it:
      # on the scaled objective a weight reads directly as that response's share,
      # which is what a hand-set vector is almost always meant to say
      w = w * user_weights
      names(w) = names(loss_spec)
      weights = w
      self$loss_baseline = base
      self$loss_weights  = w
      private$create_loss_functions(loss_spec, weights)
      cli::cli_alert_info(
        "weights: intercept-only baseline x your vector ({paste(sprintf('%s=%.3g', names(w), w), collapse = ', ')})")
    }
    # what each response actually contributes at the baseline, which is the thing a
    # hand-set vector is usually wrong about (see FINNverse/FINN#37)
    private$report_loss_shares(loss_spec, weights, base, auto_weights)

    envs = private$extract_env_method(env)

    if(!is.null(disturbance)) {
      disturbance = extract_env(~0+intensity, disturbance)
    }

    year_sequence = which(levels(as.factor(env$year)) %in% levels(as.factor(data$year)), arr.ind = TRUE)

    # A loss only fires where an observation year is ALSO an update_step boundary (see
    # forward()). If no observation year is, the objective is empty and the fit silently
    # trains on nothing - e.g. observations at 5/10/15/20 with update_step = 3.
    if(!any(year_sequence %% update_step == 0)) stop(paste("update_step", update_step, "scores no observation: none of the observation timesteps (", paste(year_sequence, collapse = ", "), ") is a multiple of it"))

    # init networks if not yet done
    private$create_nn(self$process_mortality, "mortality", dim(envs$mortality_env)[3])
    private$create_nn(self$process_growth, "growth", dim(envs$growth_env)[3])
    private$create_nn(self$process_regeneration, "regeneration", dim(envs$regeneration_env)[3])

    if(is.null(batchsize)) batchsize = nrow(envs[[1]])
    sites =  nrow(envs[[1]])

    if(is.null(init_cohort)) {
      sp = self$N_species
      self$init_cohort = CohortMat(dims = c(sites, patches, sp),
                                   dbh = array(1, dim = c(sites, patches, sp)),
                                   trees = array(1, dim = c(sites, patches, sp)),
                                   sp = sp)
    } else {
      self$init_cohort = init_cohort
      if(self$init_cohort$sp != sp) stop(paste("sp in cohort", self$init_cohort$sp, "does not match sp from data", sp))
      patches = dim(self$init_cohort$species_r)[2]
    }

    anchor_mode = match.arg(anchor_mode, c("replace", "merge"))
    # `||` coerces, the execution gates use isTRUE(): `reanchor = 1` would pass every
    # check below and then silently not re-anchor.
    for (nm_ in c("reanchor", "teacher_forcing")) {
      v_ = get(nm_)
      if (!is.logical(v_) || length(v_) != 1L || is.na(v_))
        stop(sprintf("`%s` must be TRUE or FALSE.", nm_), call. = FALSE)
    }
    if(reanchor && teacher_forcing) stop("`reanchor` and `teacher_forcing` are different estimators and cannot be combined: re-anchoring free-runs each segment from an observed state, teacher forcing never simulates one.")
    # `update_step` is inert under teacher forcing - there is no multi-year trajectory to
    # truncate - and this is the invariant that makes that true: every census is still scored.
    if(teacher_forcing && !all(year_sequence %% update_step == 0)) stop(paste("teacher_forcing = TRUE needs every observation to be scored, but update_step", update_step, "drops the observations at timestep(s)", paste(year_sequence[year_sequence %% update_step != 0], collapse = ", ")))
    # a missing anchor would silently condition an interval on a much older census. The first
    # interval is conditioned on init_cohort and nothing follows the last observation, so
    # those two are the only timesteps that need no anchor.
    if(teacher_forcing) {
      gaps = setdiff(utils::head(sort(year_sequence), -1L), as.integer(names(anchor_cohorts)))
      if(length(gaps)) stop(paste("`teacher_forcing = TRUE` needs an anchor at every observed timestep except the last; missing", paste(gaps, collapse = ", ")))
    }
    if(reanchor || (teacher_forcing && !is.null(anchor_cohorts))) {
      if(is.null(anchor_cohorts)) stop("`reanchor = TRUE` requires `anchor_cohorts`, a named list of CohortMat objects (names = simulation timesteps).")
      if(anchor_mode == "merge" && is.null(anchor_min_dbh)) stop("`anchor_mode = \"merge\"` requires `anchor_min_dbh`, the dbh below which the anchors are blind (same units as the `dbh` column).")
      if(length(anchor_cohorts) == 0 || is.null(names(anchor_cohorts)) || any(names(anchor_cohorts) == "")) stop("`anchor_cohorts` must be a non-empty list, named by the simulation timesteps to re-anchor at.")
      init_dims = self$init_cohort$dbh$shape
      for(k in names(anchor_cohorts)) {
        step = as.integer(k)
        anchor = anchor_cohorts[[k]]
        if(!(step %in% year_sequence)) stop(paste("anchor_cohorts: timestep", k, "is not in year_sequence (", paste(year_sequence, collapse = ", "), ") and would never be used"))
        # the state is only replaced where it is also truncated, i.e. at an update_step
        # boundary - an anchor anywhere else would silently never fire
        if(!teacher_forcing && step %% update_step != 0) stop(paste("anchor_cohorts: timestep", k, "is not a multiple of update_step", update_step, "and would never be used"))
        if(anchor$sp != self$init_cohort$sp) stop(paste("anchor_cohorts: timestep", k, "has sp", anchor$sp, "but init_cohort has sp", self$init_cohort$sp))
        if(anchor$dbh$shape[1] != init_dims[1]) stop(paste("anchor_cohorts: timestep", k, "has", anchor$dbh$shape[1], "sites but init_cohort has", init_dims[1]))
        if(anchor$dbh$shape[2] != init_dims[2]) stop(paste("anchor_cohorts: timestep", k, "has", anchor$dbh$shape[2], "patches but init_cohort has", init_dims[2]))
        # growth/mortality/regeneration are aggregated over (i - period + 1):i, so an anchor
        # anywhere in that window makes the aggregate mix pre- and post-anchor state
        for(t in seq_along(year_sequence)) {
          p = pl[,t]
          inside = !is.na(p) & step >= (year_sequence[t] - p + 1) & step <= (year_sequence[t] - 1)
          if(any(inside)) stop(paste("anchor_cohorts: timestep", k, "falls inside the rate-aggregation window of the observation at timestep", year_sequence[t], "(period_length", max(p[inside]), "); the aggregated growth/mortality/regeneration would mix pre- and post-anchor state"))
        }
      }
    }
    self$reanchor = reanchor
    self$teacher_forcing = teacher_forcing
    self$anchor_mode = anchor_mode
    self$anchor_min_dbh = anchor_min_dbh
    # Detached copies on the fit's device. CohortMat always builds its buffers on the CPU
    # (its `device` argument only sets device_r), and the caller's anchors must survive a
    # fit unmoved and unmodified so the same list can be reused across fits.
    self$anchor_cohorts = if(!reanchor && !teacher_forcing) NULL else lapply(anchor_cohorts, function(a)
      list(dbh     = a$dbh$detach()$to(dtype = self$dtype, device = self$device)$clone(),
           trees   = a$trees$detach()$to(dtype = self$dtype, device = self$device)$clone(),
           species = a$species$detach()$to(dtype = torch_int64(), device = self$device)$clone()))

    if(tree_res) {
      # One accumulator per cohort holds one window at a time, so a window may not reach back
      # past the census before it - a configuration in which the aggregate terms would score
      # the same year twice too, and which no inventory we know of produces.
      for(t in seq_along(year_sequence)) {
        prev = c(0L, year_sequence[year_sequence < year_sequence[t]])
        w = pl[,t]
        if(any(!is.na(w) & year_sequence[t] - w + 1 <= max(prev)))
          stop(paste("loss_resolution = \"tree\": the observation at timestep", year_sequence[t], "has period_length",
                     max(w[!is.na(w)]), "reaching back to or past the previous observation at timestep",
                     max(prev), "; the per-tree window would overlap the one before it"))
      }
      # the arrays were built from the table against the layout each interval starts from
      # (build_tree_targets), so their dimensions are right by construction; this only
      # guards the anchor validation above against a layout that changed in between
      for(k in names(tree_targets)) {
        prev = anchor_steps[anchor_steps < as.integer(k)]
        dims = if(length(prev)) anchor_cohorts[[as.character(max(prev))]]$dbh$shape else self$init_cohort$dbh$shape
        stopifnot(identical(as.integer(dim(tree_targets[[k]]$growth)), as.integer(dims)))
      }
    }
    self$loss_resolution = loss_resolution
    self$tree_data = if(!tree_res) NULL else lapply(tree_targets, function(td) lapply(td, function(a)
      torch_tensor(a, dtype = self$dtype, device = self$device)))

    if(device == "gpu") {
      device = "cuda:0"
      if(!torch::cuda_is_available()) {
        cli::cli_text("GPU device not available...")
        device = "cpu"
      }
    }

    self$to(device = device)
    if(is.null(disturbance)) {
      data = tensor_dataset(torch_tensor(envs[[1]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[2]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[3]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(Y[,,,], dtype=self$dtype, device=torch_device('cpu')),
                            torch_arange(1,dim(envs[[1]])[1])$to(dtype = torch_int64() , device=torch_device('cpu'))
      )
    } else {
      data = tensor_dataset(torch_tensor(envs[[1]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[2]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[3]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(Y[,,,], dtype=self$dtype, device=torch_device('cpu')),
                            torch_arange(1, dim(envs[[1]])[1] )$to(dtype = torch_int64() , device=torch_device('cpu')),
                            torch_tensor(disturbance, dtype=self$dtype, device=torch_device('cpu'))
      )
    }

    DataLoader = torch::dataloader(data, batch_size=batchsize, shuffle=shuffle, num_workers=0, pin_memory=TRUE, drop_last=TRUE)

    self$history = list()
    self$tree_deaths = NULL
    # written only by forward()'s tree branch, so a plain fit() after a tree one would
    # otherwise leave the first call's predictions in place with nothing marking them stale
    self$tree_pred = NULL
    self$nonfinite_counts = NULL
    self$gradients = list()

    # Rebuild the optimizer if this is the first fit() call, or if `lr` changed since
    # the optimizer was last (re)built. Previously the optimizer (and its lr) was only
    # ever constructed once per model object, so a second fit() call with a different
    # `lr` silently kept training at the old rate. Calling fit() again with the *same*
    # lr still reuses the existing optimizer instance (preserving e.g. Adam's momentum
    # state) so that resuming training for more epochs behaves as before.
    #
    # The existing optimizer is reused only if it still drives THIS model's
    # parameter tensors. After torch_save()/torch_load() it does not: the
    # optim_ignite_* optimizers come back as a dead external pointer (any access
    # errors), and the pure-R ones (optim_adam, ...) come back holding copies of
    # the old tensors, so step() would update orphans and the model would
    # silently not train. Adam's moment estimates are lost in both cases; a
    # resumed fit starts a fresh optimizer.
    if(is.null(self$optimizer) || !isTRUE(all.equal(self$optimizer_lr, lr)) ||
       !private$optimizer_owns_parameters(self$optimizer)) {
      self$optimizer = optimizer(self$parameters, lr = lr, ...)
      self$optimizer_lr = lr
    }

    # ---- Learning-rate scheduler (built fresh for this fit() call) ----
    # The scheduler tracks its own internal epoch counter starting at 0, so a
    # second fit() call restarts the chosen schedule rather than continuing a
    # previous one (consistent with epochs/checkpoints also being per-call).
    # lr_scheduler = "none" (the default) skips this entirely and reproduces
    # the constant-lr behavior of every previous version of fit().
    lr_scheduler_obj = NULL
    if (!(is.character(lr_scheduler) && length(lr_scheduler) == 1 && identical(lr_scheduler, "none"))) {
      lr_scheduler_obj = private$build_lr_scheduler(self$optimizer, lr_scheduler, lr, epochs, lr_scheduler_params)
    }

    cli::cli_progress_bar(format = "Epoch: {cli::pb_current}/{cli::pb_total} {cli::pb_bar} ETA: {cli::pb_eta} DBH: {dbh_l} BA: {ba_l} Trees: {trees_l} g: {g_l} m: {m_l} r: {r_l}", total = epochs, clear = FALSE)
    if(is.null(year_sequence)) year_sequence = 1:envs[[1]]$shape[2]

    self$train()

    for(epoch in seq_len(epochs)){
      counter = 1
      batch_loss = matrix(NA, nrow = 10000, ncol = 7L)
      tree_cov = list()
      coro::loop(for (b in DataLoader) {
        x_mort =   b[[1]]$to(device = self$device, non_blocking=TRUE)
        x_growth = b[[2]]$to(device = self$device, non_blocking=TRUE)
        x_reg =    b[[3]]$to(device = self$device, non_blocking=TRUE)
        y =        b[[4]]$to(device = self$device, non_blocking=TRUE)
        ind =      b[[5]]$to(device = self$device, non_blocking=TRUE)
        dist = NULL
        if(!is.null(disturbance))
          dist =  b[[6]]$to(device = self$device, non_blocking=TRUE)

        self$optimizer$zero_grad()

        if (is.null(self$init_cohort)) {
          cohorts = CohortMat(dims = c(x$shape[1], patches, self$sp),
                                  sp = sp,
                                  device=self$device)
          trees = cohorts$trees
          species = cohorts$species
          dbh = cohorts$dbh
        } else {
          trees = self$init_cohort$trees$to(device = self$device, non_blocking=TRUE)[ind,]
          species = self$init_cohort$species$to(device = self$device, non_blocking=TRUE)[ind,]
          dbh = self$init_cohort$dbh$to(device = self$device, non_blocking=TRUE)[ind,]
        }

        # Slice the extra observation tensors to this batch of sites, exactly as
        # init_cohort is sliced above.
        class_obs_b = if (!is.null(self$class_obs)) self$class_obs$to(device = self$device, non_blocking = TRUE)[ind, , , ] else NULL
        quant_obs_b = if (!is.null(self$quant_obs)) self$quant_obs$to(device = self$device, non_blocking = TRUE)[ind, , , ] else NULL
        tree_obs_b  = if (!is.null(self$tree_obs))
          lapply(self$tree_obs[private$TREE_OBS_TENSORS],
                 function(t) t$to(device = self$device, non_blocking = TRUE)[ind, , , ]) else NULL

        # the anchors are site-aligned state (already on device), so they must be subset by
        # the batch's site index exactly as init_cohort is
        anchors = NULL
        if(reanchor || teacher_forcing) anchors = lapply(self$anchor_cohorts, function(a)
          list(dbh = a$dbh[ind,], trees = a$trees[ind,], species = a$species[ind,]))

        # site-aligned like the anchors, flattened in forward() in the ids' column-major order
        tree_targets = NULL
        if(tree_res) tree_targets = lapply(self$tree_data, function(td) lapply(td, function(a) a[ind,]))

        pred_tmp = self$forward(dbh = dbh,
                                trees = trees,
                                species = species,
                                env = list(x_mort, x_growth, x_reg),
                                disturbance = dist,
                                start_time = start_time,
                                y = y,
                                class_obs = class_obs_b,
                                quant_obs = quant_obs_b,
                                tree_obs = tree_obs_b,
                                update_step = update_step,
                                verbose = FALSE,
                                year_sequence = year_sequence,
                                nonfinite = nonfinite,
                                reanchor = reanchor,
                                anchor_cohorts = anchors,
                                anchor_mode = anchor_mode,
                                anchor_min_dbh = anchor_min_dbh,
                                teacher_forcing = teacher_forcing,
                                loss_resolution = loss_resolution,
                                tree_data = tree_targets)

        # browser()
        pred = pred_tmp[[1]]
        loss = pred_tmp[[2]]
        # The cap's prior is NOT a loss term: forward() backwards the weighted
        # losses per step, so adding it there would put it through the
        # weights = "auto" baseline scaling and make `prior_lambda` mean something
        # different on every dataset. Backward it once per batch instead -- the
        # gradients accumulate into the same step.
        cap_pen = private$cap_prior_penalty()
        if (!is.null(cap_pen)) cap_pen$backward()
        if(!is.null(pred_tmp$tree_deaths)) tree_cov[[length(tree_cov) + 1L]] = pred_tmp$tree_deaths
        # ---- Guard before gradient clipping ----
        # Safe gradient checker: skips NULL or undefined grads
        check_grads <- function(params) {
          bad <- list()

          has_defined_grad <- function(p) {
            if (is.null(p$grad)) return(FALSE)
            ok <- tryCatch({
              # if this fails, grad is undefined
              g <- p$grad$detach()
              TRUE
            }, error = function(e) FALSE)
            ok
          }

          for (nm in names(params)) {
            p <- params[[nm]]
            if (!has_defined_grad(p)) next

            g <- p$grad$detach()

            all_finite <- as.logical(g$isfinite()$all()$item())
            if (!all_finite) {
              bad[[nm]] <- list(
                n_nan  = as.numeric(g$isnan()$sum()$item()),
                n_inf  = as.numeric(g$isinf()$sum()$item()),
                max_abs = tryCatch(as.numeric(g$abs()$max()$item()),
                                   error = function(e) NA_real_)
              )
            }
          }

          bad
        }

        bad <- check_grads(self$parameters)

        if (length(bad) > 0) {
          # Previously we only warned here and then fell through into clip_grad_norm_()
          # and optimizer$step() regardless. A NaN gradient makes the total grad norm
          # (and hence clip_coef) NaN too, which is exactly what produced the follow-on
          # crash "Error in if (clip_coef$item() < 1) : missing value where TRUE/FALSE
          # needed". Skip this batch's update entirely instead: zero the corrupted
          # gradients and move on, rather than applying (or crashing on) a bad step.
          warning("Non-finite gradients detected in epoch ", epoch,
                  " - skipping optimizer step for this batch. Affected parameters: ",
                  paste(names(bad), collapse = ", "), call. = FALSE)
          self$optimizer$zero_grad()
        } else {
          # ---- Per-parameter-group gradient clipping ----
          # A single global clip_grad_norm_() over ALL parameters mixes three
          # groups with very different natural gradient scales: the mechanistic
          # per-species/per-process rate parameters, the env-effect networks
          # (nn_growth/nn_mortality/nn_regeneration - a small NN even when not
          # "hybrid"), and the loss-distribution nuisance parameters
          # (par_loss_*_scale / par_loss_*_theta). Clipping them together means
          # whichever group happens to have the largest norm dominates the shared
          # budget and can over-clip the others. Clip each group against its own
          # budget (private$clip_grad_norm_grouped()) instead.
          private$clip_grad_norm_grouped(self$parameters, clip_norm)
          self$optimizer$step()
        }

        batch_loss[counter, ] =  as.numeric(loss$data()$cpu())
        counter <<- counter + 1
      })

      # too expensive with many parameters
      if(epoch %% checkpoints == 0) {
        self$param_history = c(self$param_history, list(lapply(self$parameters, function(p) as.matrix(p$cpu()))))
        if(record_gradients) {
          self$gradients = c(self$gradients, list(lapply(self$parameters, function(p) {
            if(prod(p$grad$shape) > 0.5) return(p$grad |> as.matrix())
          })))
        }
      }

      bl = colMeans(batch_loss, na.rm = TRUE)
      bl = round(bl, 5)
      dbh_l = bl[1]
      ba_l = bl[2]
      trees_l = bl[3]
      g_l = bl[4]
      m_l = bl[5]
      r_l = bl[6]

      #cat("Epoch: ", epoch, "Loss: ", bl, "\n")
      self$history[[epoch]] = colMeans(batch_loss, na.rm = TRUE)

      # How often the simulation killed a tree the inventory found alive, summed over this
      # epoch's batches. Every tracked cohort is scored over its whole window, so this is a
      # statement about the model's mortality, not about the sample the loss saw.
      if(length(tree_cov)) {
        cv = do.call(rbind, tree_cov)
        cv = rowsum(cv[, -1, drop = FALSE], cv[, 1])
        # a census can legitimately carry no fate at all (`died` all NA); every share below
        # would then be 0/0, which also breaks `frac_killed + frac_ghost <= 1`
        cv = cv[cv[,1] > 0, , drop = FALSE]
        self$tree_deaths = rbind(self$tree_deaths, data.frame(
          epoch = epoch, census = as.integer(rownames(cv)), n_target = cv[,1],
          frac_killed = cv[,2]/cv[,1], mean_kill_year = cv[,3]/pmax(cv[,2], 1),
          window = cv[,4]/cv[,1], frac_ghost = cv[,5]/cv[,1], row.names = NULL))
      }

      if (!is.null(lr_scheduler_obj)) {
        # "plateau" (torch::lr_reduce_on_plateau) is the one built-in that needs
        # a monitored metric instead of just an epoch count; drive it off this
        # epoch's total (summed) loss across all six response channels. Checked
        # via class rather than the raw `lr_scheduler` argument so that partial
        # string matches (match.arg() above) still resolve correctly. Every
        # other scheduler - including a user-supplied custom one - is stepped
        # with no argument.
        if (inherits(lr_scheduler_obj, "lr_reduce_on_plateau")) {
          lr_scheduler_obj$step(sum(bl, na.rm = TRUE))
        } else {
          lr_scheduler_obj$step()
        }
      }

      if(plot_progress) {
        if(epoch == 1) max_losses = self$history[[1]]*1.1
        losses = do.call(rbind, self$history)
        losses = losses/matrix(max_losses, nrow = nrow(losses), ncol = ncol(losses), byrow = TRUE)

        if(epoch == 1) {
          visualize.training(losses, epoch, epochs, "Loss", new = TRUE)
        } else {
          visualize.training(losses, epoch, epochs, "Loss", new = FALSE)
        }



      #   losses = apply(losses, 2, scales::rescale)
      #
      #   par(mfrow = c(1,1))
      #     e = tryCatch({
      #
      #        matplot(losses[,1:6], type = "l", lty = 1, col = c("#FA4D19", "#49777A", "#42DA5D","#FA8A19", "#000000", "#1992FA"),
      #                lwd = 1.5,
      #                las = 1, xlab = "Epochs", ylab = "Loss in %")
      #        legend("topright", pch = 15, bty = "n", col = c("#FA4D19", "#49777A", "#42DA5D","#FA8A19", "#000000", "#1992FA"),
      #               legend = c("dbh", "ba", "trees", "growth", "mort", "reg"))
      #
      #
      #     # plot(1:epochs, y = losses[,i],
      #     #      xlab = "Epochs",
      #     #      ylab = paste0("Loss ", c("dbh", "ba", "trees", "growth", "mort", "regeneration")[i]),
      #     #      type = "l")
      #     }, error = function(e) e)
      #
      }

      if(!is.null(folder)) {
        if(!dir.exists(folder)) dir.create(folder)
        # written to a temporary name and renamed: a kill or a full disk during the
        # save would otherwise leave a truncated .pt that torch_load() cannot read,
        # and a checkpoint stands for an hour of compute (FINNverse/FINN#36)
        if(epoch %% checkpoints == 0) {
          final = paste0(folder, "/", "epoch_", epoch, "model.pt")
          tmp   = paste0(folder, "/.epoch_", epoch, "model.pt.tmp", Sys.getpid())
          torch::torch_save(self$clone(deep=FALSE), path = tmp)
          if(!file.rename(tmp, final)) {
            unlink(tmp)
            warning(sprintf("could not write the checkpoint for epoch %d to '%s'.", epoch, final), call. = FALSE)
          }
        }
      }

      cli::cli_progress_update()

    }
    cli::cli_progress_done()


    if(torch::cuda_is_available()) torch::cuda_empty_cache()

    # Say it plainly if a response kept overflowing: under "drop_term" those
    # timesteps still trained every other response, under "skip_step" they
    # trained nothing at all.
    if (!is.null(self$nonfinite_counts) && any(self$nonfinite_counts > 0)) {
      nz = self$nonfinite_counts[seq_along(loss_spec)] > 0
      warning(sprintf("Non-finite loss %s for %s (%s).",
                      if (identical(nonfinite, "skip_step")) "discarded the whole timestep" else "skipped",
                      paste(sprintf("%s (%d timesteps)", names(loss_spec)[nz],
                                    self$nonfinite_counts[seq_along(loss_spec)][nz]), collapse = ", "),
                      if (identical(nonfinite, "skip_step")) "nonfinite = \"skip_step\": no response was scored at those timesteps"
                      else "the other responses at those timesteps were still used"),
              call. = FALSE)
    }

    # ignore debugging method; with epochs = 0 nothing was predicted
    if (epochs > 0) self$pred = list(long = pred2DF(list(Predictions = pred), "long"), wide = pred2DF(list(Predictions = pred), "wide"))
    # par is restored via on.exit() (set above, only when plot_progress = TRUE)

  },

  private = list(

    # The tensor entries of build_tree_obs()'s result (the rest is bookkeeping).
    TREE_OBS_TENSORS = c("dbh", "trees", "species", "growth", "died"),

    # Does `opt` still step on this model's own parameter tensors? FALSE for an
    # optimizer restored by torch_load(): an ignite optimizer is then a dead
    # external pointer (reading param_groups errors), a pure-R one holds copies
    # of the old tensors. identical() on two tensors compares their external
    # pointers, i.e. whether they are the same tensor.
    optimizer_owns_parameters = function(opt) {
      if (is.null(opt)) return(FALSE)
      tryCatch({
        held = unlist(lapply(opt$param_groups, function(g) g$params), recursive = FALSE)
        own = self$parameters
        length(held) == length(own) &&
          all(vapply(seq_along(own), function(i) identical(held[[i]], own[[i]]), logical(1)))
      }, error = function(e) FALSE)
    },

    # Which loss components can be backpropagated at this timestep.
    #
    # A single non-finite response used to discard the WHOLE timestep, which
    # meant one overflowing term removed every OTHER response's signal there
    # too: on the FIA fits the nbinom recruitment term overflowed at the
    # dominant 10-year interval and took most of the training signal with it.
    # Indexing the finite components (rather than masking or zeroing them) keeps
    # the offending term out of the graph entirely, so no NaN gradient can reach
    # a parameter through it. The occurrences are counted and reported by fit().
    finite_terms = function(loss) {
      ok = as.logical(as.numeric(loss$isfinite()$to(dtype = torch::torch_float32())$cpu()))
      ok[is.na(ok)] = FALSE
      if (!all(ok)) {
        cnt = as.numeric(!ok)
        self$nonfinite_counts = if (is.null(self$nonfinite_counts)) cnt else self$nonfinite_counts + cnt
      }
      which(ok)
    },

    # The components to score at this timestep under fit()'s `nonfinite` option:
    # "drop_term" keeps every finite term, "skip_step" keeps all or nothing. Both
    # count the non-finite responses (finite_terms), which fit() reports.
    scored_terms = function(loss, nonfinite) {
      keep = private$finite_terms(loss)
      if (identical(nonfinite, "skip_step") && length(keep) < length(loss)) return(integer(0))
      keep
    },

    # The summed loss to backpropagate and the detached values to report for the
    # components in `keep`. Under "skip_step" `keep` is the whole vector, and it
    # is taken exactly as it was before "drop_term" existed; under "drop_term"
    # the finite components are indexed out, so a non-finite one leaves no
    # NaN in the graph, and the report carries 0 in its place.
    loss_at_step = function(loss, keep, nonfinite) {
      if (identical(nonfinite, "skip_step"))
        return(list(sum = loss$sum(), value = loss$detach()))
      list(sum   = loss[keep]$sum(),
           value = loss$detach()$nan_to_num(nan = 0, posinf = 0, neginf = 0))
    },

    # The number of units behind a loss function's mean, with the same mask and
    # weights the loss functions use (create_loss_functions): the non-NaN
    # entries of `true`, or the sum of their clamped weights `n` (n_at_risk for
    # the binomial, growth_n for the weighted MSE). The per-site path weights
    # each step's mean by this count, which makes the call's loss the mean over
    # all its units - the same pooling the "auto" baselines use.
    scored_units = function(true, n = NULL) {
      mask = true$isnan()$bitwise_not()
      if (is.null(n)) return(as.numeric(mask$sum()$cpu()))
      as.numeric(n[mask]$clamp(min = 0)$sum()$cpu())
    },

    #=# Loss aggregation: how a response is reduced from cohorts to the observation #=#
    # Returns one entry per response: list(type =, args =). Unnamed responses keep
    # the historical rule, so omitting `loss_aggregation` changes nothing.
    resolve_loss_aggregation = function(spec, responses) {
      default_type = c(dbh = "mean", ba = "sum", trees = "sum",
                       growth = "mean", mortality = "mean", regeneration = "sum")
      supported = c("mean", "sum", "size_classes", "quantiles", "tree")
      out = lapply(responses, function(r) {
        d = if (r %in% names(default_type)) unname(default_type[r]) else "sum"
        list(type = d, args = list())
      })
      names(out) = responses
      if (is.null(spec)) return(out)
      if (is.character(spec)) spec = as.list(spec)
      if (!is.list(spec) || is.null(names(spec)))
        stop("`loss_aggregation` must be a NAMED character vector or list, one entry per response.", call. = FALSE)
      if (anyDuplicated(names(spec)))
        stop(sprintf("`loss_aggregation` names '%s' more than once; only the first entry would be used.",
                     paste(unique(names(spec)[duplicated(names(spec))]), collapse = "', '")), call. = FALSE)
      for (r in names(spec)) {
        if (!r %in% responses)
          stop(sprintf("`loss_aggregation` names '%s', which is not a response in `loss_family` (%s).",
                       r, paste(responses, collapse = ", ")), call. = FALSE)
        s = spec[[r]]
        if (is.character(s) && length(s) == 1L) {
          out[[r]] = list(type = s, args = list())
        } else if (is.list(s)) {
          # accepts list(type = "size_classes", breaks = ...) and list("size_classes", breaks = ...)
          nm = names(s)
          if (!is.null(nm) && "type" %in% nm) { type = s$type; args = s[setdiff(nm, "type")] }
          else                                { type = s[[1]];  args = s[-1] }
          out[[r]] = list(type = as.character(type), args = args)
        } else {
          stop(sprintf("`loss_aggregation$%s` must be a string or a list.", r), call. = FALSE)
        }
        if (!out[[r]]$type %in% supported)
          stop(sprintf("Unsupported loss_aggregation '%s' for '%s'. Supported: %s.",
                       out[[r]]$type, r, paste(supported, collapse = ", ")), call. = FALSE)
        if (out[[r]]$type == "size_classes" && is.null(out[[r]]$args$breaks))
          stop(sprintf("loss_aggregation '%s' = size_classes needs `breaks`, e.g. list(type = \"size_classes\", breaks = seq(20, 60, by = 10)).", r),
               call. = FALSE)
        if (out[[r]]$type == "quantiles") {
          if (!identical(r, "dbh"))
            stop("loss_aggregation 'quantiles' is implemented for 'dbh' only.", call. = FALSE)
          p = out[[r]]$args$probs
          if (is.null(p))
            stop("loss_aggregation 'dbh' = quantiles needs `probs`, e.g. list(type = \"quantiles\", probs = c(0.1, 0.5, 0.9)).",
                 call. = FALSE)
          if (any(p <= 0 | p >= 1))
            stop("`probs` must lie strictly between 0 and 1.", call. = FALSE)
        }
        if (out[[r]]$type == "tree" && !r %in% c("growth", "mortality"))
          stop("loss_aggregation 'tree' is implemented for 'growth' and 'mortality' only.", call. = FALSE)
      }
      out
    },

    # Model side of "size_classes": tree-weighted stem counts per diameter class,
    # summed over cohorts AND patches, so the result is a PLOT total per species,
    # matching an integer observed count. Class membership is a difference of two
    # sigmoids rather than a hard bucket, so a cohort straddling a break
    # contributes to both and gradients reach dbh. Reuses aggregate_results()
    # once per class, which already scatters by species and sums over patches.
    class_counts = function(dbh, trees, species, breaks, tau) {
      lower = c(-Inf, breaks); upper = c(breaks, Inf)
      per_class = lapply(seq_along(lower), function(k) {
        w = if (is.infinite(lower[k]))      ((upper[k] - dbh) / tau)$sigmoid()
            else if (is.infinite(upper[k])) 1 - ((lower[k] - dbh) / tau)$sigmoid()
            else ((upper[k] - dbh) / tau)$sigmoid() - ((lower[k] - dbh) / tau)$sigmoid()
        zero = torch::torch_zeros(list(dbh$shape[1], self$N_species), dtype = dbh$dtype, device = dbh$device)
        aggregate_results(species, list(trees * w), list(zero))[[1]]
      })
      torch::torch_stack(per_class, dim = 3)          # [sites, species, K+1]
    },

    # Model side of "quantiles": the abundance-weighted diameter quantiles of
    # each site x species, over cohorts AND patches.
    #
    # Where "size_classes" asks how MANY stems sit in a diameter interval, this
    # asks how THICK the stem at a given position of the distribution is. The
    # difference matters: a count likelihood weighs the smallest stem exactly as
    # much as the largest, so it says almost nothing about the rare big trees
    # that carry the basal area, whereas the 0.9 quantile is a statement about
    # those trees on the diameter scale itself.
    #
    # Differentiable by construction. Sorting by diameter passes gradients
    # through the sorted VALUES, and the quantile is then a soft selection along
    # that order: with the per-species cumulative abundance F running 0 -> 1, the
    # weight of cohort i for probability p is
    #     a_i = sigmoid((F_i - p)/tau) - sigmoid((F_{i-1} - p)/tau),
    # a smooth bump on whichever cohort the crossing falls in, which collapses on
    # the exact weighted quantile as tau -> 0. One sort serves every species,
    # because a per-species cumulative sum taken along a globally diameter-sorted
    # order is still ordered by diameter.
    quantile_values = function(dbh, trees, species, probs, tau) {
      d = dbh$flatten(start_dim = 2)                       # [sites, patches*cohorts]
      n = trees$flatten(start_dim = 2)
      s = species$flatten(start_dim = 2)
      # The ORDER carries no gradient, only the values do, so the sort is taken
      # without one and the values are then gathered: gather's backward scatters
      # each gradient back to the cohort it came from, which is exactly right,
      # and it keeps sort's saved indices out of the graph (R's gather shifts a
      # 1-based index tensor in place, which the sort backward then refuses).
      idx = torch::with_no_grad(torch::torch_sort(d$detach(), dim = 2)[[2]])
      d_s = d$gather(2, idx$clone())
      n_s = n$gather(2, idx$clone())
      s_s = s$gather(2, idx$clone())
      per_species = lapply(seq_len(self$N_species), function(j) {
        w   = n_s * s_s$eq(j)$to(dtype = d_s$dtype)
        tot = w$sum(dim = 2, keepdim = TRUE)$clamp(min = 1e-8)
        F   = w$cumsum(dim = 2) / tot
        Fp  = F - w / tot                                   # F_{i-1}
        q = lapply(probs, function(p) {
          a = ((F - p) / tau)$sigmoid() - ((Fp - p) / tau)$sigmoid()
          (a * d_s)$sum(dim = 2) / a$sum(dim = 2)$clamp(min = 1e-8)
        })
        torch::torch_stack(q, dim = 2)                      # [sites, n_probs]
      })
      torch::torch_stack(per_species, dim = 2)              # [sites, species, n_probs]
    },

    # Observed diameter quantiles -> [sites, obs_years, species, n_probs].
    # Columns are expected in `data` as <response>_q<100*p>, e.g. dbh_q10,
    # dbh_q50, dbh_q90 for probs = c(0.1, 0.5, 0.9).
    build_quantile_obs = function(data, probs, response = "dbh") {
      cols = sprintf("%s_q%g", response, 100 * probs)
      missing = setdiff(cols, colnames(data))
      if (length(missing))
        stop(sprintf("loss_aggregation quantiles needs columns %s in `data`; missing: %s.",
                     paste(cols, collapse = ", "), paste(missing, collapse = ", ")), call. = FALSE)
      per_prob = lapply(cols, function(cc) {
        f = stats::as.formula(paste0("~0+", cc))
        abind::abind(lapply(1:self$N_species, function(i) extract_env(f, data[data$species == i, ])), along = 3L)
      })
      torch::torch_tensor(abind::abind(per_prob, along = 4L), dtype = torch::torch_float32())
    },

    # Model side of "tree": evaluate the growth and mortality kernels on the
    # OBSERVED trees of one observation year. There is no correspondence between
    # model cohorts and observed trees, so none is attempted: the observed stand
    # replaces the simulated one for this term. Light comes from the observed
    # stand, and mortality uses the growth this same call predicts, exactly as
    # the simulation does. `obs` holds padded [sites, patches, max_trees] tensors.
    #
    # The rates are accumulated over the SAME per-site window the aggregated path
    # uses (mask `M`, one column per simulated year), because the observation is
    # an interval quantity: growth is the mean of the annual rates, mortality the
    # complement of the product of annual survivals. Evaluating a single year
    # instead would compare one year's climate - including its drought anomaly -
    # against a decade's mean. The stand is held at its state at the interval
    # start, so light and the diameters are computed once and only the
    # environment moves across the window.
    tree_predictions = function(obs, k, env, i, M, p_t) {
      dbh = obs$dbh[, k, , ]; sp = obs$species[, k, , ]; n = obs$trees[, k, , ]
      light = self$competition_func(dbh = dbh, species = sp, trees = n,
                                    parComp = self$par_competition, h = NULL,
                                    patch_size_ha = self$patch_size_ha, ba = NULL,
                                    cohortHeights = NULL,
                                    n_quantiles = self$n_quantiles, continuous = self$continuous,
                                    light_mask = isTRUE(self$light_mask))
      # only the years some site actually has in its window
      p_max  = as.numeric(p_t$max()$item())
      j_from = max(1L, i - as.integer(ceiling(p_max)) + 1L)
      g_sum = NULL; lsurv = NULL
      for (j in j_from:i) {
        w3 = M[, j]$reshape(c(-1, 1, 1))
        pred_g = if (inherits(self$process_growth, "hybrid")) env[["growth"]][, j, ]
                 else index_species(self$nn_growth(env[["growth"]][, j, ]), sp)
        g = self$growth_func(dbh = dbh, species = sp, parGrowth = self$par_growth,
                             pred = pred_g, light = light, trees = n)
        pred_m = if (inherits(self$process_mortality, "hybrid")) env[["mort"]][, j, ]
                 else index_species(self$nn_mortality(env[["mort"]][, j, ]), sp)
        m = self$mortality_func(dbh = dbh, species = sp, trees = n, parMort = self$par_mortality,
                                pred = pred_m, light = light, growth = g)
        # Out-of-window years contribute 0 to the growth sum and a survival of 1.
        # The survival product is accumulated in LOG space: the direct form,
        # prod (1 - m_j)^w_j, needs a tensor exponent, and torch's backward for
        # that is not defined for the mask's dtype.
        g_sum = if (is.null(g_sum)) g * w3 else g_sum + g * w3
        ls_j  = (1 - m$clamp(1e-6, 1 - 1e-6))$log() * w3
        lsurv = if (is.null(lsurv)) ls_j else lsurv + ls_j
      }
      denom = p_t$clamp(min = 1)$reshape(c(-1, 1, 1))
      list(growth = g_sum / denom, mortality = 1 - lsurv$exp())
    },

    # Observed stem counts per diameter class -> [sites, obs_years, species, K+1].
    # Columns are expected in `data` as <prefix>1 ... <prefix>K+1, one row per
    # site x year x species, holding PLOT totals so the counts are integers.
    build_class_obs = function(data, K, prefix = "n_class_") {
      cols = paste0(prefix, seq_len(K))
      missing = setdiff(cols, colnames(data))
      if (length(missing))
        stop(sprintf("loss_aggregation size_classes needs columns %s in `data`; missing: %s.",
                     paste(cols, collapse = ", "), paste(missing, collapse = ", ")), call. = FALSE)
      per_class = lapply(cols, function(cc) {
        f = stats::as.formula(paste0("~0+", cc))
        abind::abind(lapply(1:self$N_species, function(i) extract_env(f, data[data$species == i, ])), along = 3L)
      })
      torch::torch_tensor(abind::abind(per_class, along = 4L), dtype = torch::torch_float32())
    },

    # loss_resolution = "tree": the observed tree table -> one pair of target
    # arrays (`growth`, `mort`) per scored census, each with the dimensions of
    # the cohort layout that census's interval STARTS from: `init_cohort`, or the
    # anchor at the previous census under reanchor / teacher_forcing. Element
    # [s, p, c] is the tree that is cohort c of (site s, patch p) in that layout,
    # NaN where no tree was observed. `year` maps to a timestep exactly as in
    # build_tree_obs (the k-th observation year is year_sequence[k]).
    #
    # Strict by design: the targets are keyed by position, so a row that cannot
    # be placed is an error rather than a dropped row - a misaligned table would
    # otherwise score every tree against another one without any sign of it.
    build_tree_targets = function(tree_data, site_ids, obs_years, year_sequence,
                                  init_cohort, anchor_cohorts, anchor_steps) {
      d = data.table::as.data.table(tree_data)
      need = c("siteID", "year", "patchID", "cohortID", "species", "growth", "died")
      missing = setdiff(need, colnames(d))
      if (length(missing))
        stop(sprintf("loss_resolution = \"tree\": `tree_data` needs columns %s; missing: %s.",
                     paste(need, collapse = ", "), paste(missing, collapse = ", ")), call. = FALSE)
      bad_year = setdiff(unique(d$year), obs_years)
      if (length(bad_year))
        stop(sprintf("loss_resolution = \"tree\": `tree_data` has rows in year(s) %s, which carry no observation in `data` (scored years: %s); every row must be a tree of a scored census.",
                     paste(bad_year, collapse = ", "), paste(obs_years, collapse = ", ")), call. = FALSE)
      bad_site = setdiff(unique(d$siteID), site_ids)
      if (length(bad_site))
        stop(sprintf("loss_resolution = \"tree\": `tree_data` names site(s) %s that are not in `data`.",
                     paste(utils::head(bad_site, 10), collapse = ", ")), call. = FALSE)
      for (cc in c("patchID", "cohortID")) {
        v = d[[cc]]
        if (anyNA(v) || any(v != round(v)) || any(v < 1))
          stop(sprintf("loss_resolution = \"tree\": `tree_data$%s` must be positive integers (the %s index of the cohort layout).",
                       cc, if (cc == "patchID") "patch" else "cohort"), call. = FALSE)
      }
      # the layout each interval starts from
      layout = function(src, what) {
        if (is.null(src) || is.null(src$dbh) || is.null(src$trees) || is.null(src$species))
          stop(sprintf("loss_resolution = \"tree\": %s does not carry dbh, trees and species arrays.", what), call. = FALSE)
        list(dbh = torch::as_array(src$dbh$cpu()), trees = torch::as_array(src$trees$cpu()),
             species = torch::as_array(src$species$cpu()))
      }
      out = list()
      for (k in seq_along(obs_years)) {
        rows = d[d$year == obs_years[k], ]
        if (!nrow(rows)) next
        step = year_sequence[k]
        prev = anchor_steps[anchor_steps < step]
        L = if (length(prev)) layout(anchor_cohorts[[as.character(max(prev))]], sprintf("anchor_cohorts[[\"%d\"]]", max(prev)))
            else layout(init_cohort, "`init_cohort`")
        dims = dim(L$trees)
        s = match(rows$siteID, site_ids); p = as.integer(rows$patchID); cid = as.integer(rows$cohortID)
        where = sprintf("(siteID %s, patchID %s, cohortID %s, year %s)", rows$siteID, rows$patchID, rows$cohortID, rows$year)
        off = p > dims[2] | cid > dims[3]
        if (any(off))
          stop(sprintf("loss_resolution = \"tree\": %d tree(s) of census year %s name a cohort outside the %d x %d (patches x cohorts) layout the interval starts from, e.g. %s. `cohortID` is the position along the cohort dimension of %s.",
                       sum(off), obs_years[k], dims[2], dims[3], where[which(off)[1]],
                       if (length(prev)) sprintf("the anchor at timestep %d", max(prev)) else "`init_cohort`"), call. = FALSE)
        idx = cbind(s, p, cid)
        dup = duplicated(idx)
        if (any(dup))
          stop(sprintf("loss_resolution = \"tree\": %d tree(s) of census year %s name a cohort another row already uses, e.g. %s; one cohort is one tree, so `cohortID` must be unique within (siteID, patchID, year). Binned initial cohorts (makeInitCohorts(dbh_binsize = )) cannot be scored per tree.",
                       sum(dup), obs_years[k], where[which(dup)[1]]), call. = FALSE)
        empty = L$trees[idx] < 0.5
        if (any(empty))
          stop(sprintf("loss_resolution = \"tree\": %d tree(s) of census year %s name a cohort that does not exist (no trees) in the layout the interval starts from, e.g. %s.",
                       sum(empty), obs_years[k], where[which(empty)[1]]), call. = FALSE)
        sp_bad = L$species[idx] != rows$species
        if (any(sp_bad))
          stop(sprintf("loss_resolution = \"tree\": %d tree(s) of census year %s have a species other than their cohort's, e.g. %s is species %s but cohort %s is species %s. The table and %s do not describe the same trees - were the initial cohorts binned (makeInitCohorts(dbh_binsize = ))?",
                       sum(sp_bad), obs_years[k], where[which(sp_bad)[1]], rows$species[which(sp_bad)[1]],
                       cid[which(sp_bad)[1]], L$species[idx][which(sp_bad)[1]],
                       if (length(prev)) "the anchor" else "`init_cohort`"), call. = FALSE)
        if ("dbh" %in% colnames(rows)) {
          # obsDF2arrays() rounds to 4 decimals and the tensor is float32
          have = is.finite(rows$dbh)
          d_bad = have & abs(L$dbh[idx] - rows$dbh) > 0.01 + 1e-4 * abs(rows$dbh)
          if (any(d_bad))
            stop(sprintf("loss_resolution = \"tree\": %d tree(s) of census year %s have a dbh other than their cohort's at the interval start, e.g. %s has dbh %.3f but cohort %s has %.3f. The table and %s do not describe the same trees - were the initial cohorts binned (makeInitCohorts(dbh_binsize = ))?",
                         sum(d_bad), obs_years[k], where[which(d_bad)[1]], rows$dbh[which(d_bad)[1]],
                         cid[which(d_bad)[1]], L$dbh[idx][which(d_bad)[1]],
                         if (length(prev)) "the anchor" else "`init_cohort`"), call. = FALSE)
        }
        g = array(NaN, dim = dims); mort = array(NaN, dim = dims)
        g[idx] = as.numeric(rows$growth); mort[idx] = as.numeric(rows$died)
        g[is.na(g)] = NaN; mort[is.na(mort)] = NaN
        out[[as.character(step)]] = list(growth = g, mort = mort)
      }
      message(sprintf("[tree_data] %d trees keyed to their initial cohorts at %d census(es)", nrow(d), length(out)))
      out
    },

    # Observed trees -> padded [sites, obs_years, patches, max_trees] tensors.
    # Padding slots carry trees = 0 (no weight in competition) and NaN responses,
    # which the loss functions already mask out.
    #
    # The patch dimension is the number of distinct patchIDs in the table, not
    # the `patches` of fit()/predictTrees(). Light is computed per patch from the
    # trees in it, so patches holding no tree change no prediction; sized by
    # `patches` (25 in the BWI fits, where every tree has patchID 1) the tensors
    # carried 24 empty patches through competition, growth and mortality in
    # every batch.
    build_tree_obs = function(tree_data, site_ids, obs_years) {
      # copy(): as.data.table() returns the caller's own object when it already
      # is a data.table, and everything below adds columns by reference
      d = data.table::copy(data.table::as.data.table(tree_data))
      d[, .row := .I]          # so a prediction can be handed back on the caller's rows
      need = c("siteID", "year", "patchID", "species", "dbh")
      missing = setdiff(need, colnames(d))
      if (length(missing))
        stop(sprintf("`tree_data` needs columns %s; missing: %s.",
                     paste(need, collapse = ", "), paste(missing, collapse = ", ")), call. = FALSE)
      if (!"trees"  %in% colnames(d)) d[, trees  := 1]
      if (!"growth" %in% colnames(d)) d[, growth := NA_real_]
      if (!"died"   %in% colnames(d)) d[, died   := NA_real_]
      d[, s := match(siteID, site_ids)]
      d[, y := match(year, obs_years)]
      d[, p := match(patchID, sort(unique(patchID)))]
      d = d[!is.na(s) & !is.na(y) & !is.na(p)]
      if (!nrow(d)) stop("`tree_data` has no rows matching the sites and observation years of `data`.", call. = FALSE)
      data.table::setorder(d, s, y, p)
      d[, slot := seq_len(.N), by = .(s, y, p)]
      dims = c(length(site_ids), length(obs_years), max(d$p), max(d$slot))
      idx  = cbind(d$s, d$y, d$p, d$slot)
      fill_arr = function(col, fill) { a = array(fill, dim = dims); a[idx] = d[[col]]; a }
      message(sprintf("[tree_data] %d trees, up to %d per patch and observation year", nrow(d), dims[4]))
      # `idx` and `rows` say which tensor slot each input row went into, so a
      # prediction can be handed back on the caller's own rows (predict_trees).
      list(idx     = idx,
           rows    = d[, .(.row, siteID, year, patchID, s, y, p, slot)],
           dbh     = torch::torch_tensor(fill_arr("dbh", 0), dtype = torch::torch_float32()),
           trees   = torch::torch_tensor(fill_arr("trees", 0), dtype = torch::torch_float32()),
           species = torch::torch_tensor(fill_arr("species", 1L), dtype = torch::torch_int64()),
           growth  = torch::torch_tensor(fill_arr("growth", NaN), dtype = torch::torch_float32()),
           died    = torch::torch_tensor(fill_arr("died", NaN), dtype = torch::torch_float32()))
    },

    # Splits a (named) list of parameters into three groups by name prefix:
    #   - "loss":        par_loss_*_scale / par_loss_*_theta (loss-distribution
    #                    nuisance parameters; never touched by the forward
    #                    simulation itself)
    #   - "nn":          nn_growth / nn_mortality / nn_regeneration (the
    #                    env-effect network for each process - a small NN even
    #                    when the process is not a "hybrid" one)
    #   - "mechanistic": everything else (par_growth*, par_mortality*,
    #                    par_regeneration*, par_competition*, par_theta_recruits*,
    #                    ...) - the per-species/per-process mechanistic rates.
    # Name-based so this degrades gracefully for models that don't have all
    # three kinds of parameters (e.g. a purely mechanistic, non-hybrid model
    # ends up with everything in "mechanistic" and behaves like the old single
    # global clip).
    split_parameter_groups = function(params) {
      nm = names(params)
      is_loss = grepl("^par_loss_", nm)
      is_nn   = grepl("^nn_", nm)
      list(
        mechanistic = params[!is_loss & !is_nn],
        nn          = params[is_nn],
        loss        = params[is_loss]
      )
    },
    # Clips each parameter group (see split_parameter_groups()) to its own
    # gradient-norm budget instead of clipping the global norm across every
    # parameter at once - see the call site in fit() for why this matters.
    # `clip_norm` is either a single number (applied to every group) or a
    # named list/vector keyed by "mechanistic"/"nn"/"loss" (any group not
    # named falls back to 2.0).
    clip_grad_norm_grouped = function(params, clip_norm) {
      groups = private$split_parameter_groups(params)
      get_budget = function(group_name) {
        if (is.list(clip_norm) || !is.null(names(clip_norm))) {
          # coerce a named numeric vector to a list first: `[[` on an atomic
          # vector errors for a missing name, whereas on a list it returns
          # NULL, which is what the fallback below expects.
          cn = if (is.list(clip_norm)) clip_norm else as.list(clip_norm)
          val = cn[[group_name]]
          if (!is.null(val) && !is.na(val)) return(val)
          return(2.0)
        }
        return(clip_norm[1])
      }
      for (group_name in names(groups)) {
        g = groups[[group_name]]
        if (length(g) > 0) torch::nn_utils_clip_grad_norm_(g, get_budget(group_name))
      }
      invisible(NULL)
    },
    # Builds the lr_scheduler object requested via fit(lr_scheduler = ...).
    # `type` is either one of the built-in names below, or a
    # function(optimizer) returning a custom torch lr_scheduler (must
    # implement $step()). `params` overrides individual defaults, which are
    # otherwise scaled off `epochs`/`lr` so they're reasonable without tuning.
    build_lr_scheduler = function(optimizer, type, lr, epochs, params) {
      if (is.function(type)) return(type(optimizer))

      type = match.arg(type, c("step", "exponential", "cosine", "plateau"))

      if (type == "step") {
        defaults = list(step_size = max(1L, floor(epochs / 4)), gamma = 0.5)
        a = utils::modifyList(defaults, params)
        return(torch::lr_step(optimizer, step_size = a$step_size, gamma = a$gamma))

      } else if (type == "exponential") {
        # lr_multiplicative() multiplies the *current* lr by the factor on every
        # step(), so a constant factor < 1 gives the usual lr_t = lr_0 * gamma^t
        # exponential decay.
        defaults = list(gamma = 0.97)
        a = utils::modifyList(defaults, params)
        local_gamma = a$gamma
        return(torch::lr_multiplicative(optimizer, lr_lambda = function(epoch) local_gamma))

      } else if (type == "cosine") {
        defaults = list(T_max = epochs, eta_min = lr * 0.01)
        a = utils::modifyList(defaults, params)
        return(torch::lr_cosine_annealing(optimizer, T_max = a$T_max, eta_min = a$eta_min))

      } else if (type == "plateau") {
        defaults = list(factor = 0.5, patience = max(5L, floor(epochs / 20)), mode = "min")
        a = utils::modifyList(defaults, params)
        return(torch::lr_reduce_on_plateau(optimizer, mode = a$mode, factor = a$factor, patience = a$patience))
      }
    },
    # Loss of an INTERCEPT-ONLY model: what each response's loss would be if we
    # predicted the single best constant. This is the null deviance, and it is
    # what `weights = "auto"` divides by, so every term becomes a fraction of its
    # own baseline ("how much better than the mean").
    #
    # The constant used is the MLE intercept, which is the mean for every family
    # here (Gaussian/MSE, Poisson, negative binomial) and the pooled proportion
    # sum(n*y)/sum(n) for the binomial.
    #
    # The real loss functions are reused rather than reimplemented, so this
    # cannot drift from them, and any future family is handled for free. They are
    # called with 2D [sites, species] tensors during fitting (the reg/nbinom one
    # sizes theta off pred$shape[1]), so Y's [sites, years, species] slices are
    # flattened to match.
    # One line naming what each response actually contributes, because a hand-set
    # weight vector acts on losses whose raw magnitudes differ by orders of
    # magnitude and nothing used to say so (FINNverse/FINN#37).
    report_loss_shares = function(loss_spec, weights, base, auto) {
      if (is.null(base) || !any(is.finite(base))) return(invisible(NULL))
      contrib = weights * base
      ok = is.finite(contrib) & contrib > 0
      if (!any(ok)) return(invisible(NULL))
      share = contrib / sum(contrib[ok], na.rm = TRUE)
      cli::cli_alert_info(
        "share of the objective at the intercept-only baseline ({if (auto) 'scaled' else 'raw scale'}): {paste(sprintf('%s=%.2g%%', names(loss_spec)[ok], 100*share[ok]), collapse = ', ')}")
      invisible(share)
    },

    compute_baseline_losses = function(Y, loss_spec, class_obs = NULL, tree_obs = NULL, quant_obs = NULL) {
      nm  = names(loss_spec)
      out = stats::setNames(rep(NA_real_, 6), nm)
      failed = list()
      # Responses whose observation does not live in Y need their own
      # intercept-only baseline, otherwise weights = "auto" silently falls back
      # to 1 for them and unbalances every other term.
      alt = list()
      if (!is.null(class_obs) && identical(self$loss_agg$trees$type, "size_classes"))
        alt$trees = list(obs = class_obs, f = self$loss_trees_func)
      if (!is.null(quant_obs) && identical(self$loss_agg$dbh$type, "quantiles"))
        alt$dbh = list(obs = quant_obs, f = self$loss_dbh_func)
      if (!is.null(tree_obs) && identical(self$loss_agg$growth$type, "tree"))
        alt$growth = list(obs = tree_obs$growth, f = self$loss_growth_func)
      if (!is.null(tree_obs) && identical(self$loss_agg$mortality$type, "tree"))
        alt$mortality = list(obs = tree_obs$died, f = self$loss_mortality_func)
      torch::with_no_grad({
        for (i in 1:6) {
          # a response that is switched off has nothing to scale by, and saying so
          # for every such response is noise rather than a warning
          if (identical(as.character(loss_spec[[i]]), "none")) next
          out[i] = tryCatch({
            flat = function(t) t$reshape(c(-1, t$shape[length(t$shape)]))
            a    = alt[[nm[i]]]
            true = if (is.null(a)) flat(Y[,,,i]) else a$obs
            # slice 8 = n_at_risk (binomial), slice 9 = growth_n (weighted MSE).
            # The baseline must use the SAME weighting as the loss it scales.
            # The alternative observations are already unaggregated (one entry
            # per tree, or per size class), so they carry no n.
            n = if (!is.null(a)) NULL
                else if (identical(nm[i], "mortality")) flat(Y[,,,8])
                else if (identical(nm[i], "growth")) flat(Y[,,,9])
                else NULL
            mask = true$isnan()$bitwise_not()
            if (!as.logical(mask$max()$data())) stop("no observations")
            yv = true[mask]
            # A constant observation has a baseline of a few float ulps rather than
            # 0, so the `w <= 0` guard below never fires and 1/baseline explodes.
            # The tree branch already refuses a one-class target for this reason.
            if (as.numeric(yv$max() - yv$min()) == 0)
              stop("the observation is constant, so the intercept-only model is exact")
            mu = if (is.null(n)) {
              as.numeric(yv$mean())
            } else {
              wv = n[mask]$clamp(min = 0)
              as.numeric((yv * wv)$sum()) / max(as.numeric(wv$sum()), 1e-8)
            }
            if (!is.finite(mu)) stop("non-finite intercept")
            pred = torch::torch_full_like(true, mu)
            f = if (is.null(a)) self[[paste0("loss_", nm[i], "_func")]] else a$f
            v = if (is.null(n)) f(true, pred) else f(true, pred, n)
            as.numeric(v)
          }, error = function(e) {
            # A baseline we cannot compute must not silently become a weight of
            # Inf/NaN; fall back to 1 (i.e. leave that term unscaled) and say so --
            # once for all of them, because up to six separate warnings per fit is
            # noise rather than information.
            failed[[nm[i]]] <<- conditionMessage(e)
            NA_real_
          })
        }
      })
      if (length(failed))
        warning(sprintf("weights = 'auto': no baseline for %s; left unscaled.",
                        paste(sprintf("'%s' (%s)", names(failed), unlist(failed)), collapse = ", ")),
                call. = FALSE)
      out
    },

    create_loss_functions = function(loss, weights) {
      for(l in 1:6) {
          tmp_loss = loss[l]
          tmp_loss_name = names(loss)[l]
          # `func` is assigned inside the if/else chain below and only stored at
          # the end of the iteration. An unrecognised name would leave `func`
          # holding the PREVIOUS variable's closure and silently fit the wrong
          # likelihood, so reject it here instead.
          # "none" switches a response off entirely, which is how a response is
          # retired without changing the six-slot layout (e.g. dbh and ba once
          # diameter information enters through size-classed stem counts).
          supported = c("mse", "gaussian", "poisson", "binomial", "nbinom", "none")
          if(!tmp_loss %in% supported) {
            stop(sprintf("Unsupported loss '%s' for '%s'. Supported: %s.",
                         tmp_loss, tmp_loss_name, paste(supported, collapse = ", ")))
          }
          if(tmp_loss == "none") {
            self[[paste0("loss_", tmp_loss_name, "_func")]] = function(true, pred, n = NULL) 0.0
            next
          }
          # "nbinom" only has an implementation for these two responses.
          if(tmp_loss == "nbinom" && !tmp_loss_name %in% c("trees", "regeneration")) {
            stop(sprintf("loss 'nbinom' is only implemented for 'trees' and 'regeneration', not '%s'.",
                         tmp_loss_name))
          }
          func = NULL
          if(tmp_loss == "mse") {
            func =local({
              local_l = l
              local_weights = weights
              function(true, pred, n = NULL) {
                mask = true$isnan()$bitwise_not()
                if(!as.logical(mask$max()$data())) return(0.0)
                if(is.null(n)) return(torch::nnf_mse_loss(pred[mask], true[mask])$mean()*local_weights[local_l])
                # `true` is a MEAN of n trees, so its variance is sigma^2/n: an
                # observation from 32 trees is worth 32x one from a single tree.
                # Weighting the squared error by n is the Gaussian counterpart of
                # weighting the binomial by n_at_risk. Normalised by sum(w) so the
                # term stays O(1) and comparable to the unweighted form.
                w = n[mask]$clamp(min = 0)
                se = (pred[mask] - true[mask])$pow(2)
                return(((se*w)$sum()/w$sum()$clamp(min = 1e-8))*local_weights[local_l])
              }
            })
          } else if(tmp_loss == "binomial") {
            # Bernoulli/binomial negative log-likelihood for a proportion.
            # `pred` is a survival-probability complement in (0,1) (the processes
            # return it through a sigmoid), and `true` is the observed fraction
            # that died. nnf_binary_cross_entropy accepts fractional targets, so
            # this is the binomial log-likelihood up to the constant binomial
            # coefficient - the term that does not depend on the parameters.
            func = local({
              local_l = l
              local_weights = weights
              function(true, pred, n = NULL) {
                mask = true$isnan()$bitwise_not()
                if(as.logical(mask$max()$data())) {
                  # clamp p away from {0,1}: log(0) is -Inf and would poison the
                  # gradient of every parameter, not just this loss term.
                  p = pred[mask]$clamp(1e-6, 1 - 1e-6)
                  y = true[mask]$clamp(0, 1)
                  if(is.null(n)) return(torch::nnf_binary_cross_entropy(p, y)*local_weights[local_l])
                  # Weighting each observation by its trees-at-risk is what turns
                  # the Bernoulli term into the binomial one: it is identical to
                  # a cbind(died, survived) response, and to glm(prop ~ .,
                  # family = binomial, weights = n). Normalising by sum(w) rather
                  # than taking torch's mean(w * l) keeps the term O(1), so the
                  # `weights` argument stays on the same scale as before.
                  w = n[mask]$clamp(min = 0)
                  nll = torch::nnf_binary_cross_entropy(p, y, weight = w, reduction = "sum")
                  denom = w$sum()$clamp(min = 1e-8)
                  return((nll/denom)*local_weights[local_l])
                }
                else return(0.0)
              }
            })
          } else if(tmp_loss == "gaussian") {

            if(is.null(self[[paste0("par_loss_",tmp_loss_name, "_scale")]])) {
              self[[paste0("par_loss_",tmp_loss_name, "_scale")]] = torch::nn_parameter(torch_tensor(0.55))
            }

            func = local({
              local_tmp_loss_name = tmp_loss_name
              local_l = l
              local_weights = weights
              function(true, pred, n = NULL) {
                sigma_raw = self[[paste0("par_loss_",local_tmp_loss_name, "_scale")]]
                sigma = torch_log(1.0+torch_exp(sigma_raw) )
                mask = true$isnan()$bitwise_not()
                if(as.logical(mask$max()$data())) {
                  return(torch::distr_normal(pred[mask], scale = sigma)$log_prob(true[mask])$negative()$mean()*local_weights[local_l] + 0.001*(sigma**2))
                } else { return(0.0) }
              }
            })
          } else if(tmp_loss == "poisson") {
            func = local({
                local_l = l
                local_weights = weights
                # `trees_offset` keeps log_prob finite at a simulated zero; the
                # historical 1e-5 leaves the GRADIENT at y/1e-5, see FINN#53.
                local_offset = if (is.null(self$trees_offset)) 1e-5 else self$trees_offset
                function(true, pred, n = NULL) {
                  mask = true$isnan()$bitwise_not()
                  if(as.logical(mask$max()$data())) return(torch::distr_poisson(pred[mask] + local_offset)$log_prob(true[mask])$negative()$mean()*local_weights[local_l])
                  else return(0.0)
                }
              })
          } else if (tmp_loss == "nbinom") {
            if (tmp_loss_name == "trees") {

              if(is.null(self[[paste0("par_loss_",tmp_loss_name, "_theta")]] )) {
                self[[paste0("par_loss_",tmp_loss_name, "_theta")]] = torch::nn_parameter(torch_tensor(0.5413))
              }

              func =
                local({
                  local_tmp_loss_name = tmp_loss_name
                  local_l = l
                  local_weights = weights
                  local_offset = if (is.null(self$trees_offset)) 0.001 else self$trees_offset
                  function(true, pred, n = NULL) {
                    mask = true$isnan()$bitwise_not()
                    theta = 1.0/torch::nnf_softplus( self[[paste0("par_loss_",local_tmp_loss_name, "_theta")]] )
                    if(as.logical(mask$max()$data())) return(dnbinom_torch(pred[mask]+local_offset, true[mask], theta)$mean()*local_weights[local_l])
                    else return(0.0)
                }
              })
            } else if(tmp_loss_name == "regeneration") {
              func =
                local({
                  local_l = l
                  local_weights = weights
                  function(true, pred, n = NULL) {
                    mask = true$isnan()$bitwise_not()
                    theta = self$par_theta_recruits
                    theta = theta$squeeze(1)$`repeat`(c(pred$shape[1], 1))
                    ## Observation operator (inventory simulator): map both the
                    ## predicted recruitment DENSITY (r_mean_ha) and the observed
                    ## target into the sampling design's OBSERVATION space by the
                    ## inclusion weight `recruit_obs_weight` (w) before the nbinom
                    ## likelihood. For an angle-count design w = g(recruit_dbh)/BAF,
                    ## so density (0..~9000/ha) becomes a stable tallied-recruit
                    ## COUNT (0..~9). Default w = 1 reproduces the previous per-ha
                    ## behaviour exactly (fixed-area inventories, FINN-fia).
                    w = if (is.null(self$recruit_obs_weight)) 1.0 else self$recruit_obs_weight
                    if(as.logical(mask$max()$data())) return(dnbinom_torch(pred[mask]*w+0.001, true[mask]*w, theta[mask])$mean()*local_weights[local_l])
                    else return(0.0)
                  }
                })
            }
          }
          self[[paste0("loss_", tmp_loss_name, "_func")]] = func
      }
    },

    create_nn = function(obj, type, inputs) {
      hybrid = inherits(obj, "hybrid")
      if(is.null(self[[paste0("nn_", type)]])) {
        if(!hybrid) {
          nn =
            if(is.null(obj$NN)) build_NN(input_shape = inputs, output_shape = self$N_species, bias = TRUE, activation = "selu", hidden = obj$hidden, dropout = obj$dropout, last_activation = "linear")
            else obj$NN

          if(!is.null(obj$initEnv)) {
            for(i in 1:length(nn$parameters)) nn$parameters[[i]]$set_data( obj$initEnv[[i]] )
          }

        } else {
          if(type == "mortality") {
            if(obj$transformer) {
              nn = hybrid_transformer(num_species = self$N_species,
                                      num_env_vars = inputs,
                                      dgtl_embedder_dim = 4L,
                                      max_len = 500L,
                                      emb_dim=obj$emb_dim,
                                      num_heads=1L,
                                      num_layers=obj$encoder_layers,
                                      dropout=obj$dropout,
                                      dim_feedforward = obj$dim_feedforward)
            } else {
               nn = hybrid_DNN(num_species = self$N_species,
                               num_env_vars = inputs+2,
                               emb_dim=obj$emb_dim,
                               dropout=obj$dropout,
                               hidden = obj$hidden)
            }
          }

          if(type == "growth") {
            if(obj$transformer) {
            nn = hybrid_transformer(num_species = self$N_species,
                                    num_env_vars = inputs,
                                    dgtl_embedder_dim = 3L,
                                    max_len = 500L,
                                    emb_dim=obj$emb_dim,
                                    num_heads=1L,
                                    num_layers=obj$encoder_layers,
                                    dropout=obj$dropout,
                                    dim_feedforward = obj$dim_feedforward)
            } else {
              nn = hybrid_DNN(num_species = self$N_species,
                              num_env_vars = inputs,
                              emb_dim=obj$emb_dim,
                              dropout=obj$dropout,
                              hidden = obj$hidden)
            }

          }

          if(type == "regeneration") {
            # TODO: I guess for regeneration we should only recommend the non-transformer option, always the same number numbers of cohorts etc, maybe add a warning or message informing the user about it??

              nn = hybrid_DNN(num_species = self$N_species,
                              num_env_vars = inputs,
                              emb_dim=obj$emb_dim,
                              dropout=obj$dropout,
                              hidden = obj$hidden,
                              regeneration = TRUE)

          }

        }

        if(!obj$optimizeEnv) {
          for(p in nn$parameters) p$requires_grad_(FALSE)
        }
        self[[paste0("nn_", type)]] = nn
      }
    },

    extract_env_method = function(env) {
      # if the model was fitted with env_autoscale = TRUE, re-apply the SAME stored
      # standardization to the (raw) env before building the design matrices.
      if (!is.null(self$env_scaling)) env = apply_env_scaling(env, self$env_scaling)
      mortality_env = extract_env(self$mortality_formula, env)
      growth_env = extract_env(self$growth_formula, env)
      regeneration_env = extract_env(self$regeneration_formula, env)
      return(list(mortality_env=mortality_env, growth_env=growth_env, regeneration_env=regeneration_env))
    },

    add_process = function(obj, type = c("mortality", "growth", "regeneration", "competition")) {
      type = match.arg(type)

      hybrid = FALSE
      if(inherits(obj, "hybrid")) {
        hybrid = TRUE
        func = switch(type, growth= { growth_hybrid }, mortality = { mortality_hybrid }, regeneration = { regeneration_hybrid })
        obj$func = func

      }

      if(is.null(obj)) {
        func = switch(type, mortality = { mortality }, growth = { growth }, regeneration = {regeneration_capped}, competition = { competition })
        obj = createProcess(func = func)
      }

      private$setup_species_parameters(obj, type, hybrid)

      # Train env model or not
      self[[paste0("env_", type, "_optimized")]] = obj$optimizeEnv
      self[[paste0("process_", type)]] = obj

      # process specific parameters, e.g. theta for nbinom sampling of the recruits
      if(type == "regeneration") {
        self$sample_regeneration = obj$sample_regeneration
        self$regeneration_sampler = obj$regeneration_sampler
        self$par_theta_recruits = obj$dispersion_parameter
        # The Beverton-Holt cap is FINN's to create, not the user's to declare:
        # `cap` sets its resolution, and fit() replaces the placeholder start value
        # with the resolved prior (see resolve_cap_prior()).
        if(private$uses_capped_regeneration(obj$func)) {
          self$reg_cap = obj$cap
          self$reg_prior = obj$prior
          self$reg_prior_lambda = obj$prior_lambda
          # setup_species_parameters() ran first, so a cap declared the old way
          # through custom_parameters is already registered and stays authoritative
          # -- regeneration_adult() still needs that route for its other parameters.
          if(is.null(self$reg_logK)) {
            n_cap = if(identical(obj$cap, "species")) self$N_species else 1L
            self$register_parameter("reg_logK",
              torch::nn_parameter(torch::torch_full(n_cap, log(800), dtype = torch::torch_float32())))
          } else {
            self$reg_cap = if(self$reg_logK$numel() > 1L) "species" else "all"
          }
        }
      }

      if(type == "mortality") {
        self$mortality_draw = obj$mortality_draw
        self$mortality_ceiling = obj$mortality_ceiling
        self$mortality_floor = obj$mortality_floor
        # Set as the function's formal, not passed down forward(): the kernel reads
        # `self$record_raws` from its own environment, so a user wrapper around
        # FINN::mortality() cannot be rebound by set_environment() and is not a
        # reachable opt-out. NULL keeps the func's own default, so a custom
        # mortality func without the argument is untouched.
        if(!is.null(obj$growth_scale) && "growth_scale" %in% names(formals(self$mortality_func))) {
          fn = self$mortality_func
          formals(fn)$growth_scale = obj$growth_scale
          self$mortality_func = fn
        }
      }

      if(type == "competition") {
        self$n_quantiles = obj$n_quantiles
        self$continuous = obj$continuous
        self$light_mask = obj$light_mask
      }

    },
    # A user may wrap or rename the process function, so identify the capped
    # regeneration by the parameter it reads rather than by `identical()`.
    uses_capped_regeneration = function(func) {
      if(is.null(func)) return(FALSE)
      any(grepl("reg_logK", deparse(body(func)), fixed = TRUE))
    },

    # `prior = "max_obs_recruits"` is resolved here and not in finn(), because the
    # data only exist at fit() time. The resolved value is both the prior's centre
    # and the cap's starting value (FINNverse/FINN#58).
    resolve_cap_prior = function(obs_reg_per_ha_yr) {
      if(is.null(self$reg_logK)) return(invisible(NULL))
      # Resolve ONCE. fit() is called again to continue a calibration (cluster
      # segments, a second call with more epochs), and re-resolving would throw
      # the trained cap away and restart it at the prior.
      if(!is.null(self$reg_prior_value)) return(invisible(self$reg_prior_value))
      prior = self$reg_prior
      if(identical(prior, "max_obs_recruits")) {
        if(length(obs_reg_per_ha_yr) == 0L || all(!is.finite(obs_reg_per_ha_yr)))
          stop("regeneration_capped(prior = \"max_obs_recruits\") needs observed ",
               "recruitment to resolve the prior, and this fit scores none. Pass a ",
               "number instead, or score regeneration.", call. = FALSE)
        prior = max(obs_reg_per_ha_yr[is.finite(obs_reg_per_ha_yr)])
        if(!(prior > 0))
          stop("regeneration_capped(prior = \"max_obs_recruits\") resolved to ", prior,
               "; the observed recruitment is all zero, so the data cannot place the cap. ",
               "Pass a number instead.", call. = FALSE)
      }
      self$reg_prior_value = as.numeric(prior)
      self$reg_logK$set_data(torch::torch_full(self$reg_logK$shape, log(self$reg_prior_value),
                                               dtype = torch::torch_float32()))
      invisible(self$reg_prior_value)
    },

    # prior_lambda * (log K - log prior)^2 per cap, added AFTER the weights = "auto"
    # baseline scaling so its strength does not depend on the data's null deviance.
    # With cap = "species" the target is the differentiable mean of the current log
    # caps, so the prior penalises how far the species spread and leaves the level
    # to the data.
    cap_prior_penalty = function() {
      if(is.null(self$reg_logK) || is.null(self$reg_prior_value)) return(NULL)
      lambda = if(is.null(self$reg_prior_lambda)) 0 else self$reg_prior_lambda
      if(!(lambda > 0)) return(NULL)
      target = if(identical(self$reg_cap, "species") && self$reg_logK$numel() > 1L)
                 self$reg_logK$mean() else log(self$reg_prior_value)
      lambda * (self$reg_logK - target)$pow(2)$sum()
    },

    setup_species_parameters = function(obj, type, hybrid) {
      self[[paste0(type, "_func")]] = private$set_environment(obj$func)
      self[[paste0(type, "_formula")]] = obj$formula

      # ---- custom (extra) process parameters, declared via createProcess() -----
      # Registered here (during model construction) so they land in self$parameters
      # -> optimised by fit(), saved/restored by torch_save/torch_load, and reachable
      # from the (self-bound) custom func as self[[name]] -- no post-hoc attach or
      # reload-repair needed. See createProcess(custom_parameters = ...).
      if(!is.null(obj$custom_parameters)) {
        for(pname in names(obj$custom_parameters)) {
          spec = obj$custom_parameters[[pname]]
          if(!is.null(self[[pname]]) || is_reserved_model_field(pname))
            stop(sprintf("custom_parameters name '%s' clashes with a model field (one FINN sets itself, now or later in construction or fit()); choose another.", pname), call. = FALSE)
          init_t = torch::torch_tensor(as.numeric(spec$init), dtype = torch::torch_float32())
          if(isTRUE(spec$trainable)) self$register_parameter(pname, torch::nn_parameter(init_t))
          else                       self$register_buffer(pname, init_t)
        }
      }
      self$register_buffer(paste0("par_", type, "_upper"), torch::torch_tensor(get_par_boundary(obj, type, upper = TRUE)))
      self$register_buffer(paste0("par_", type, "_lower"), torch::torch_tensor(get_par_boundary(obj, type, upper = FALSE)))

      self[[paste0("par_", type, "_optimized")]] = obj$optimizeSpecies

      # to keep things simple, in case of hybrid, we will just ignore the species parameters for now....TODO!

      # if null, random initialisierung
      # else recalclate from init values the required values
      if(is.null(obj$initSpecies)) {
        # mortality's growth coefficient is initialised on the scale growth enters on;
        # a custom kernel without that formal keeps the historical range (scale 1)
        gs = if (!is.null(obj$growth_scale)) obj$growth_scale
             else if ("growth_scale" %in% names(formals(obj$func))) eval(formals(obj$func)$growth_scale)
             else 1
        if (!is.numeric(gs) || length(gs) != 1L || !is.finite(gs) || gs <= 0) gs = 1
        self[[paste0("par_", type)]] = init_species_parameters(type, self$N_species, growth_scale = gs)
      } else {
        self[[paste0("par_", type)]] = obj$initSpecies
      }


    },
    set_environment = function(fn) {
      environment(fn) = self$.__enclos_env__
      return(fn)
    }

  ),
  active = list(
    # TODO: reduce code redundancy!
    par_mortality = function(value) {
      if(missing(value)) {
        return(forward(self$par_mortality_unconstrained, self$par_mortality_upper, self$par_mortality_lower))
      } else {
        init = torch::torch_tensor( backward(value,  self$par_mortality_upper, self$par_mortality_lower) )
        if(self$par_mortality_optimized) self$par_mortality_unconstrained = torch::nn_parameter(init)
        else self$register_buffer("par_mortality_unconstrained", init) # if not trainable, the parameter should be created via register_buffer
      }
    },
    par_growth = function(value) {
      if(missing(value)) {
        return(forward(self$par_growth_unconstrained, self$par_growth_upper, self$par_growth_lower))
      } else {
        init = torch::torch_tensor( backward(value,  self$par_growth_upper, self$par_growth_lower) )
        if(self$par_growth_optimized) self$par_growth_unconstrained = torch::nn_parameter(init)
        else self$register_buffer("par_growth_unconstrained", init)
      }
    },
    par_regeneration = function(value) {
      if(missing(value)) {
        return(forward(self$par_regeneration_unconstrained, self$par_regeneration_upper, self$par_regeneration_lower))
      } else {
        init = torch::torch_tensor( backward(value,  self$par_regeneration_upper, self$par_regeneration_lower) )
        if(self$par_regeneration_optimized) self$par_regeneration_unconstrained = torch::nn_parameter(init)
        else self$register_buffer("par_regeneration_unconstrained", init)
      }
    },
    par_competition = function(value){
      if(missing(value)) {
        return(forward(self$par_competition_unconstrained, self$par_competition_upper, self$par_competition_lower))
      } else {
        init = torch::torch_tensor( backward(value,  self$par_competition_upper, self$par_competition_lower) )
        if(self$par_competition_optimized) self$par_competition_unconstrained = torch::nn_parameter(init)
        else self$register_buffer("par_competition_unconstrained", init)
      }
    },

    par_theta_recruits = function(value) {
      if(missing(value)) {
        return(1.0/(torch::nnf_softplus(self$par_theta_recruits_raw)+0.0001))
      } else {
        value = 1.0 / value - 0.0001
        if(length(value) == 1) value = rep(value, self$N_species)
        self$par_theta_recruits_raw = torch::nn_parameter( torch_tensor(log( exp(value) - 1.0 )) )
      }
    },

    par_mortality_r = function() {
      return( as.matrix(self$par_mortality))
    },
    par_growth_r = function() {
      return( as.matrix(self$par_growth))
    },
    par_regeneration_r = function() {
      return( as.matrix(self$par_regeneration))
    },
    par_competition_r = function() {
      return( as.matrix(self$par_competition))
    },
    parameters_r = function() {
      return(lapply(self$parameters, function(p) as.matrix(p)))
    }
  )
)


# gg = createProcess(func = growth, optimizeSpecies = TRUE)
#
# m = finn(N_species = 5, growth_process = gg, mortality_process = createProcess(func = mortality))
# library(data.table)
# sim_env_data <- data.table(expand.grid(siteID = 1:10, year = 1:5, temp = 0.9, precip = 0.7))
#
# m = finn(N_species = 5, growth_process = gg, mortality_process = createProcess(func = mortality))
# t = m$simulate(env = sim_env_data, batchsize = 10L)
#
# data = t$wide$site
# m = finn(N_species = 5, growth_process = gg, mortality_process = createProcess(func = mortality))
# m$fit(data = data, env = sim_env_data, lr = 0.001, epochs = 20L, folder = "checkpoints", checkpoints = 5L)




# period_length

# plot(t$wide$site$dbh)
#
# list(long = list(site = rbindlist( lapply(t, function(X) X$long$site ))),
#      wide = list(site = rbindlist( lapply(t, function(X) X$wide$site )))
#      )
# rbindlist( lapply(t, function(X) X$long$site ))
# rbindlist( lapply(t, function(X) X$wide$site ))
#
# nn2 = torch::torch_load( torch::torch_serialize(m) )
# nn2$simulate(sim_env_data)

#
# m$fit(X, Y, .....)
#
# m$fit(...)
#
# # von den Prozessen, brauchen wir:
# # function
# # parameters
# # data model
# # formula
# # limits/parameter ranges der parameter?!
#
#
#
#
#
#
#
# ff = FINN(Nspecies = 5)
# ff$mortality
#
# #
# A = 5
# ff = function() A
# library(torch)
# Test = nn_module("Test",
#                  initialize = function(ff = NULL) {
#                    self$ff = ff
#                    self$nn = nn_sequential(nn_linear(5, 5))
#                    self$par = nn_parameter(torch_tensor(5.0, requires_grad = TRUE))
#                    self$register_buffer("g", torch_tensor(5.0))
#                    self[["Test"]] = nn_parameter(torch_tensor(1.0))
#                    self$egal = NULL
#                    self$optimizer = torch::optim_ignite_adagrad(params = self$par)
#
#                  },
#                  forward = function(x) torch::torch_serialize(self$clone()),
#                  active = list(
#                    testPar = function(value) {
#                      if(missing(value)) {return(self$par+5.0)}
#                      else self$par = nn_parameter(torch_tensor(value))
#                    }
#                  )
# )
# #
# nn = Test(ff)
# torch_save(nn, path = "test.pt")
# nn$forward()
# class(nn)
#
# nn2 = torch_load("test.pt")
# nn2$parameters

# nn$parameters
# nn$test_NN
#
# nn$to(device = "mps")
#
# nn$optimizer$param_groups
# nn$test2 = nn_parameter(torch_tensor(3.0))
#
# nn$testPar = 10.0
# nn$testPar
#
# nn$parameters
# nn$g =torch_tensor(3.0, requires_grad = TRUE)
# nn$g
#
# nn$ff()
# nn$par
#
# nn2 = torch::torch_load( torch::torch_serialize(m) )
# nn2$parameters
# nn2$optimizer$param_groups
#
# torch::torch_save(nn, "nn.RDS")
# nn2 = torch_load("nn.RDS")
# nn2$ff()
# nn2$par
#
#
# binomial_coefficient = function(n, x){
#   return(torch::torch_exp(torch::torch_lgamma((n + 1.0)) -
#                             torch::torch_lgamma((x + 1.0)) -
#                             torch::torch_lgamma((n - x + 1.0))))
# }
#
# binomial_likelihood = function(x, n, p) {
#   binom_coeff = binomial_coefficient(n, x)
#   likelihood = binom_coeff * (p ** x) * ((1 - p) ** (n - x))
#   return(likelihood$log())
# }
# binomial_likelihood(7, 10, 0.3)
# dbinom(7, 10, 0.3, log = TRUE)

