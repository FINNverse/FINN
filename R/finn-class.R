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
finn = function(N_species,
                mortality_process = NULL,
                growth_process = NULL,
                regeneration_process = NULL,
                competition_process = NULL,
                recruits_dbh = 1.0) {
  finn_class(N_species = N_species,
             mortality_process = mortality_process,
             growth_process = growth_process,
             regeneration_process = regeneration_process,
             competition_process = competition_process,
             recruits_dbh = recruits_dbh)
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
#' a continuous rate. `mortality` is an observed *proportion* of trees that died
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
#' @param env (`data.table|data.frame`)\cr Data with environmental covariates must be passed as `data.table` or `data.frame`.
#' @param disturbance (`data.table|data.frame`)\cr Data with disturbance rates must be passed as `data.table` or `data.frame`.
#' @param patches (`integer(1)`)\cr Number of patches.
#' @param patch_size (`numeric(1)`)\cr Patch size.
#' @param init_cohort (`CohortMat`)\cr Initial cohort matrix of class `CohortMat`, created by [FINN::CohortMat]
#' @param epochs (`integer(1)`)\cr Number of iteration steps.
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
#' @param loss (`character(6)`)\cr Named vector of the different losses. Names should be `dbh`, `ba`, `trees`, `growth`, `mortality`, and `regeneration`. Supported losses are `mse`, `poisson`, `nbinom`, `gaussian`, and `binomial`. `binomial` is a Bernoulli/binomial negative log-likelihood intended for `mortality`; it expects both the prediction and the observation to be proportions in \[0, 1\].
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
#'   **Weights must account for the raw scale of each loss, not just its
#'   importance.** The six terms are summed, and their raw magnitudes differ by
#'   orders of magnitude: on the bundled FIA data `dbh` is an MSE in cm^2 (~430
#'   at the start of training) while `growth` is an MSE on a ratio (~0.03) — a
#'   factor of ~1e4. Weights chosen by importance alone are therefore dominated
#'   by whichever response happens to have the largest units, and the rest
#'   receive too little gradient to learn from. With equal weights, `dbh` takes
#'   ~87% of the FIA objective and `growth` ~1%; `growth` then never improves.
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
#'   Passing a `numeric(6)` uses those weights as-is and disables the scaling.
#'   On the FIA example, balancing the objective is worth ~+0.24 Spearman on
#'   held-out `growth` — far more than tuning `lr` or `epochs`. See the
#'   "Fitting FINN to forest inventory data" vignette.
#' @param optimizer (`torch_optimizer_generator`)\cr Optimizer from the `torch` package.
#' @param batchsize (`integer(1)`)\cr Batch size, model will be trained in random batch sizes of the data to preserve memory and improve convergence.
#' @param device (`character(1)`)\cr Should the model be fitted on the CPU or the GPU (Graphic card). Support is only for NVIDIA GPUs available.
#' @param update_step (`integer(1)`)\cr Number of steps for which the gradient should be calculated. Automatic differentation becomes slow for larger update steps and the risk of vanishing gradients increases.
#' @param start_time (`integer(1)`)\cr Starting from which year should the model be fitted. Can be used to use on burn-in.
#' @param plot_progress (`logical(1)`)\cr Plot fitting progress (losses) or not.
#' @param folder (`character(1)`)\cr Path to folder for saving checkpoint models. If `NULL`, no models will be saved during the training.
#' @param checkpoints (`integer(1)`)\cr Interval size in epochs for saving checkpoint models.
#' @param shuffle (`logical(1)`)\cr Shuffle data or not.
#' @param record_gradients (`logical(1)`)\cr Record the gradients of all parameters or not. Can get large for many epochs.
#' @param env_autoscale (`logical(1)`)\cr If `TRUE`, FINN z-standardizes the
#'   environmental predictors in `env` internally: the per-variable mean and
#'   standard deviation are learned from the training `env` and stored on the
#'   model, then re-applied automatically at every `predict()`/`simulate()` call.
#'   This lets you pass raw (untransformed) `env` for both calibration and
#'   prediction; FINN guarantees an identical transformation at both stages.
#'   Recommended (and the default) for numerical stability when predictors are on
#'   different scales. Set `FALSE` to use `env` exactly as supplied (e.g. if you
#'   have already standardized it yourself). The learned constants are available
#'   as `model$env_scaling`.
#' @param clip_norm (`numeric(1)|list()`)\cr Gradient-norm budget passed to
#'   `torch::nn_utils_clip_grad_norm_()`, applied separately to each of three
#'   parameter groups (mechanistic per-species rates, env-effect networks,
#'   loss-distribution nuisance parameters) rather than once globally across all
#'   parameters. A single number (default `2.0`) applies the same budget to
#'   every group; a named list/vector keyed by `"mechanistic"`/`"nn"`/`"loss"`
#'   overrides individual groups, e.g. `clip_norm = list(loss = 5, nn = 1)`.
#' @param ... Additional arguments passed to `optimizer`.
#'
#' @return The fitted `model`, invisibly. `fit()` trains the model in place, so
#'   the returned object is the one passed in, now carrying the training results
#'   (e.g. `$history`, `$loss_weights`, `$loss_baseline`).
#' @export
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
               loss = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"), #
               weights = "auto",
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
                      weights = weights,
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
                      ...))
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
#' @param device (`character(1)`)\cr Should the simulation run on the CPU or the GPU (Graphics card). Support is only available for NVIDIA GPUs.
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
#' @param device (`character(1)`)\cr Should the simulation run on the CPU or the GPU (Graphics card). Support is only available for NVIDIA GPUs.
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
    recruits_dbh = 1.0
  ) {
    self$N_species = N_species
    self$recruits_dbh = recruits_dbh
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
                     year_sequence = NULL){


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

    # create cohort ids
    cohort_ids = torch_tensor(array(
      1:(prod(species$shape)+1),
      dim = species$shape), dtype=torch_int32(), device = self$device
    )

    # init Result tensors
    Result = lapply(1:7,function(tmp) torch::torch_zeros(list(sites, time, self$N_species), device=self$device))
    names(Result) =  c("dbh","ba", "trees", "growth", "mort", "reg", "r_mean_ha")

    # running aggregate of the per-update_step loss, used so that the loss reported
    # back to fit() reflects all observation points hit during this forward() call,
    # not just whichever one happened to be evaluated last in the time loop.
    loss_total = torch::torch_zeros(7L, device = self$device)
    loss_count = 0L

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

      # empty rate objects/tensors
      light = torch_zeros(list(sites, time,  dbh$shape[3]), device=self$device)
      g = torch_zeros(list(sites, time, dbh$shape[3]), device=self$device)
      m = torch_zeros(list(sites, time, dbh$shape[3]), device=self$device)
      r = torch_zeros(list(sites, time, dbh$shape[3]), device=self$device)
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
          continuous = self$continuous
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
          continuous = self$continuous
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

        trees_dead = binomial_from_gamma(torch::torch_clamp(trees+trees$le(0.5)$float()+0.01, min = 1.0) , torch::torch_clamp(m, 0.01, 0.99))*trees$ge(0.5)$float()
        trees_dead = trees_dead + trees_dead$round()$detach() - trees_dead$detach()
        trees_before = trees
        # trees_dead = (trees*m)*trees$ge(0.5)$float()
        # trees_dead = trees_dead + trees_dead$round()$detach() - trees_dead$detach()
        # trees_before = trees

        #.unsqueeze(3) # TODO check!
        #trees$sub_(m)$clamp_(min = 0.0)
        trees = torch_clamp(trees - trees_dead, min = 0) #### TODO if trees = 0 then NA...prevent!
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
        continuous = self$continuous
      )

      # browser()
      ### Regeneration ####
      if(!inherits(self$process_regeneration, "hybrid"))pred = pred_reg
      else pred = env[["reg"]][,i,]

      r_mean_ha = self$regeneration_func(species = species,
                                          parReg = self$par_regeneration[,1],
                                          pred = pred,
                                          light = AL_reg)

      r_mean_patch = r_mean_ha*self$patch_size_ha

      if(self$sample_regeneration){
        r = rnbinom_torch(r_mean_patch, self$par_theta_recruits)
      }else if(!self$sample_regeneration){
        r = r_mean_patch
      }

      r = r + r$round()$detach() - r$detach()

      # New recruits
      dbh_new = ((r-1+0.1)/1e-3)$sigmoid() + (self$recruits_dbh - 1.0) # TODO: check!!! --> when r 0 dann dbh = 0, ansonsten dbh = 1 dbh[r==0] = 0
      trees_new = r
      species_new = torch_arange(1, sp, dtype=torch_int64(), device = self$device)$unsqueeze(1)$`repeat`(c(r$shape[1], r$shape[2], 1))

      # assign cohortIDs
      max_id = max(c(1,as_array(cohort_ids), na.rm = TRUE))
      new_cohort_id = torch_tensor(array(
        (max_id+1):(max_id+prod(r$shape)+1),
        dim = r$shape), dtype=torch_int32(), device = self$device
      ) #TODO check for performance

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
        # Masks to find alive cohorts
          mask = (trees > 0.5)$flatten(start_dim = 1, end_dim = 2)
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
          dbh = dbh*trees$gt(0.5)$float()
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

      ##### Rework #######
      # NOTE: truncation of the recurrent cohort state (below) is intentionally
      # decoupled from `i %in% year_sequence` - see comment further down.
      update_boundary = i > 0 && dbh$shape[3] != 0 && !is.null(y) && (i %% update_step == 0)
      if(update_boundary) {
        if(i %in% year_sequence) {
          tmp_index = which(year_sequence %in% i, arr.ind = TRUE)
          # browser()
          # #dbh
          loss[1] = self$loss_dbh_func(y[, tmp_index,,1], Result[[1]][,i,] )
          # ba
          loss[2] = self$loss_ba_func(y[, tmp_index,,2], Result[[2]][,i,] )
          # counts
          loss[3] = self$loss_trees_func(y[,tmp_index,,3], Result[[3]][,i,])


          # growth rates - check for NA in period_length, if not, then accumulate gradients?
          # currently we assume that they are constant over sites
          accumulate_gradients = y[,tmp_index,,7] |> as.matrix()
          period = unique(accumulate_gradients[,1])
          if(!is.na(period)) {
            # loss[1] = self$loss_dbh_func(y[, tmp_index,,1], Result[[1]][,(i-period+1):(i),]$mean(2) )
            # # ba
            # loss[2] = self$loss_ba_func(y[, tmp_index,,2], Result[[2]][,(i-period+1):(i),]$mean(2) )
            # # counts
            # loss[3] = self$loss_trees_func(y[,tmp_index,,3], Result[[3]][,(i-period+1):(i),]$mean(2))

            #browser()
            loss[4] = self$loss_growth_func(y[,tmp_index,,4], Result[[4]][,(i-period+1):(i),]$mean(2), y[,tmp_index,,9])
            # mort rates
            loss[5] = self$loss_mortality_func(y[,tmp_index,,5], Result[[5]][,(i-period+1):(i),]$mean(2), y[,tmp_index,,8])
            # reg rates ha
            loss[6] = self$loss_regeneration_func(y[,tmp_index,,6], Result[[7]][,(i-period+1):(i),]$sum(2) )#(Result[[7]][,(i-period+1),] - Result[[7]][,i,])$clamp(min = 0.0)  )
            self$obs_rec = y[,tmp_index,,6] |> as.matrix()
            self$pred_rec = Result[[7]][,(i-period+1):(i),]$sum(2) |> as.matrix()
            self$loss_raw = as.numeric(loss)
            # ---- Check loss before backward ----
            # (this branch used to backward() unconditionally, unlike the per-step
            # branch below; a non-finite period-averaged loss would otherwise inject
            # NaN gradients into every parameter with no diagnostic at all)
            if (!as.logical(loss$isfinite()$all()$item())) {
              cat("\n>>> Non-finite (period-averaged) loss detected at time step", i, "\n")
              for (k in 1:length(loss)) {
                cat("loss[", k, "] = ", as.numeric(loss[k]$item()),
                    " finite=", as.logical(loss[k]$isfinite()$item()), "\n")
              }
            } else {
              loss$sum()$backward()
              for(j in 1:7) Result[[j]] = Result[[j]]$detach()
              loss_total = loss_total + loss$detach()
              loss_count = loss_count + 1L
            }
          } else {

            # loss[1] = self$loss_dbh_func(y[, tmp_index,,1], Result[[1]][,i,] )
            # # ba
            # loss[2] = self$loss_ba_func(y[, tmp_index,,2], Result[[2]][,i,] )
            # # counts
            # loss[3] = self$loss_trees_func(y[,tmp_index,,3], Result[[3]][,i,])

            loss[4] = self$loss_growth_func(y[,tmp_index,,4], Result[[4]][,i,], y[,tmp_index,,9])
            # mort rates
            loss[5] = self$loss_mortality_func(y[,tmp_index,,5], Result[[5]][,i,], y[,tmp_index,,8])
            # reg rates ha
            loss[6] = self$loss_regeneration_func(y[,tmp_index,,6], Result[[7]][,i,])
            self$loss_raw = as.numeric(loss)
            # ---- Check loss before backward ----
            if (!as.logical(loss$isfinite()$all()$item())) {
              cat("\n>>> Non-finite loss detected at time step", i, "\n")

              # Inspect each component
              for (k in 1:length(loss)) {
                cat("loss[", k, "] = ", as.numeric(loss[k]$item()),
                    " finite=", as.logical(loss[k]$isfinite()$item()), "\n")
              }

              # Drop into debug mode

              # You can choose to skip backward() and continue:
              # next
            } else {
              loss$sum()$backward()
              for(j in 1:7) Result[[j]] = Result[[j]]$detach()
              loss_total = loss_total + loss$detach()
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
        dbh=dbh$detach()
        trees=trees$detach()
        species=species$detach()
        cohort_ids=cohort_ids$detach()
      }

      if(!is.null(y)) {
        #for(j in 1:3) Result[[j]] = Result[[j]]$detach()
        # if period_length = NA on all sites, it is assumed that we need only yearly gradients
        if(as.numeric(y[,,,7]$isnan()$bitwise_not()$sum()) < 0.5) for(j in 1:7) Result[[j]] = Result[[j]]$detach()

      }
      loss$detach_()

      if(verbose) cli::cli_progress_update()

    }

    # browser()

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

    # (parameters' requires_grad is restored via on.exit() above, even if this
    # function exits early due to an error)

    return(Result_out)
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
    device = match.arg(device)
    if(device == "gpu") device="cuda:0"
    self$device = device
    self$patch_size_ha = patch_size
    # browser()
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
      # batch's global site indices (`ind`) so that runs split across more than one
      # batch keep correct, non-colliding siteIDs. For the default single batch
      # (batchsize = all sites) this is the identity.
      ind_r = as.integer(torch::as_array(ind$cpu()))
      for (nm in intersect(names(long_b), c("site", "patch", "cohort"))) {
        long_b[[nm]][, siteID := ind_r[siteID]]
        wide_b[[nm]][, siteID := ind_r[siteID]]
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
                 loss = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"), #
                 weights = "auto",
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
    if(!any(loss %in% c("mse", "gaussian", "poisson", "nbinom"))) stop("Loss not supported")

    # Internal z-standardization of environmental predictors. Learn the
    # centre/scale from the (raw) training env now and store them on the model;
    # extract_env_method() re-applies them here and at every predict/simulate
    # call, so raw env can be supplied throughout. See compute_env_scaling().
    if (isTRUE(env_autoscale)) self$env_scaling = compute_env_scaling(env)

    if(is.null(data)) {
      print("No data. Switching into simulation modus...")
      return(self$predict(env = env, patches = patches, patch_size = patch_size, batchsize=batchsize, device=device))
    }

    self$train_env         = data.table::copy(data.table::as.data.table(env))
    self$train_init_cohort = if(is.null(init_cohort)) NULL else list(cohort = init_cohort)
    self$fit_id            = as.numeric(Sys.time())
    self$conditional_effects = NULL
    self$ale                 = NULL
    self$perm_importance     = NULL

    # setup loss functions
    # weights = "auto" (the default): epoch 1 runs unweighted so each loss term
    # can be MEASURED, then the weights are set to the inverse of those
    # magnitudes (see the rescale block in the epoch loop). A numeric vector is
    # used as-is and disables the rescaling.
    # `loss` (the argument) names the loss FAMILY per response. The batch loop
    # below reuses the name `loss` for a torch tensor of the current losses,
    # which shadows this argument - so snapshot the spec now, and use the
    # snapshot anywhere the families are needed later (e.g. the "auto" rescale).
    loss_spec = loss
    auto_weights = is.character(weights) && length(weights) == 1L && identical(weights, "auto")
    if (auto_weights) {
      weights = rep(1, 6)
    } else {
      if (!is.numeric(weights) || length(weights) != 6L) {
        stop("`weights` must be \"auto\" or a numeric vector of length 6 (dbh, ba, trees, growth, mortality, regeneration).")
      }
      self$loss_weights = weights
    }
    # NOTE: the loss closures are built further down, AFTER the response tensor Y
    # exists - "auto" needs the observed responses to compute the intercept-only
    # baseline it scales by.

    # setup data
    device = match.arg(device)
    if(device == "gpu") device="cuda:0"
    self$device = device
    self$patch_size_ha = patch_size

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
    # na.action restored on function exit via on.exit() set above.

    Y = torch::torch_cat(lapply(response, function(x) torch::torch_tensor(x, dtype=torch::torch_float32(), device="cpu")$unsqueeze(4)), 4)

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
    private$create_loss_functions(loss_spec, weights)   # unweighted for now
    if (auto_weights) {
      base = private$compute_baseline_losses(Y, loss_spec)
      w = 1/base
      w[!is.finite(w) | w <= 0] = 1        # no data / degenerate -> leave alone
      names(w) = names(loss_spec)
      weights = w
      self$loss_baseline = base
      self$loss_weights  = w
      private$create_loss_functions(loss_spec, weights)
      cli::cli_alert_info(
        "weights = 'auto': scaled by the intercept-only baseline ({paste(sprintf('%s=%.3g', names(w), w), collapse = ', ')})")
    }

    envs = private$extract_env_method(env)

    if(!is.null(disturbance)) {
      disturbance = extract_env(~0+intensity, disturbance)
    }

    year_sequence = which(levels(as.factor(env$year)) %in% levels(as.factor(data$year)), arr.ind = TRUE)

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
    self$gradients = list()

    # Rebuild the optimizer if this is the first fit() call, or if `lr` changed since
    # the optimizer was last (re)built. Previously the optimizer (and its lr) was only
    # ever constructed once per model object, so a second fit() call with a different
    # `lr` silently kept training at the old rate. Calling fit() again with the *same*
    # lr still reuses the existing optimizer instance (preserving e.g. Adam's momentum
    # state) so that resuming training for more epochs behaves as before.
    if(is.null(self$optimizer) || !isTRUE(all.equal(self$optimizer_lr, lr))) {
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

    for(epoch in 1:epochs){
      counter = 1
      batch_loss = matrix(NA, nrow = 10000, ncol = 7L)
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

        pred_tmp = self$forward(dbh = dbh,
                                trees = trees,
                                species = species,
                                env = list(x_mort, x_growth, x_reg),
                                disturbance = dist,
                                start_time = start_time,
                                y = y,
                                update_step = update_step,
                                verbose = FALSE,
                                year_sequence = year_sequence)

        # browser()
        pred = pred_tmp[[1]]
        loss = pred_tmp[[2]]
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
        if(epoch %% checkpoints == 0) torch::torch_save(self$clone(deep=FALSE), path = paste0(folder, "/", "epoch_", epoch, "model.pt"))
      }

      cli::cli_progress_update()

    }
    cli::cli_progress_done()


    if(torch::cuda_is_available()) torch::cuda_empty_cache()

    # ignore debugging method
    self$pred = list(long = pred2DF(list(Predictions = pred), "long"), wide = pred2DF(list(Predictions = pred), "wide"))
    # par is restored via on.exit() (set above, only when plot_progress = TRUE)

  },

  private = list(
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
    compute_baseline_losses = function(Y, loss_spec) {
      nm  = names(loss_spec)
      out = stats::setNames(rep(NA_real_, 6), nm)
      torch::with_no_grad({
        for (i in 1:6) {
          out[i] = tryCatch({
            flat = function(t) t$reshape(c(-1, t$shape[length(t$shape)]))
            true = flat(Y[,,,i])
            # slice 8 = n_at_risk (binomial), slice 9 = growth_n (weighted MSE).
            # The baseline must use the SAME weighting as the loss it scales.
            n = if (identical(nm[i], "mortality")) flat(Y[,,,8])
                else if (identical(nm[i], "growth")) flat(Y[,,,9])
                else NULL
            mask = true$isnan()$bitwise_not()
            if (!as.logical(mask$max()$data())) stop("no observations")
            yv = true[mask]
            mu = if (is.null(n)) {
              as.numeric(yv$mean())
            } else {
              wv = n[mask]$clamp(min = 0)
              as.numeric((yv * wv)$sum()) / max(as.numeric(wv$sum()), 1e-8)
            }
            if (!is.finite(mu)) stop("non-finite intercept")
            pred = torch::torch_full_like(true, mu)
            f = self[[paste0("loss_", nm[i], "_func")]]
            v = if (is.null(n)) f(true, pred) else f(true, pred, n)
            as.numeric(v)
          }, error = function(e) {
            # A baseline we cannot compute must not silently become a weight of
            # Inf/NaN; fall back to 1 (i.e. leave that term unscaled) and say so.
            warning(sprintf("weights = 'auto': could not compute a baseline for '%s' (%s); leaving it unscaled.",
                            nm[i], conditionMessage(e)), call. = FALSE)
            NA_real_
          })
        }
      })
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
          supported = c("mse", "gaussian", "poisson", "binomial", "nbinom")
          if(!tmp_loss %in% supported) {
            stop(sprintf("Unsupported loss '%s' for '%s'. Supported: %s.",
                         tmp_loss, tmp_loss_name, paste(supported, collapse = ", ")))
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
                function(true, pred, n = NULL) {
                  mask = true$isnan()$bitwise_not()
                  if(as.logical(mask$max()$data())) return(torch::distr_poisson(pred[mask] + 0.00001)$log_prob(true[mask])$negative()$mean()*local_weights[local_l])
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
                  function(true, pred, n = NULL) {
                    mask = true$isnan()$bitwise_not()
                    theta = 1.0/torch::nnf_softplus( self[[paste0("par_loss_",local_tmp_loss_name, "_theta")]] )
                    if(as.logical(mask$max()$data())) return(dnbinom_torch(pred[mask]+0.001, true[mask], theta)$mean()*local_weights[local_l])
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
                    if(as.logical(mask$max()$data())) return(dnbinom_torch(pred[mask]+0.001, true[mask], theta[mask])$mean()*local_weights[local_l])
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
            else nn

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
        func = switch(type, mortality = { mortality }, growth = { growth }, regeneration = {regeneration}, competition = { competition })
        obj = createProcess(func = func)
      }

      private$setup_species_parameters(obj, type, hybrid)

      # Train env model or not
      self[[paste0("env_", type, "_optimized")]] = obj$optimizeEnv
      self[[paste0("process_", type)]] = obj

      # process specific parameters, e.g. theta for nbinom sampling of the recruits
      if(type == "regeneration") {
        self$sample_regeneration = obj$sample_regeneration
        self$par_theta_recruits = obj$dispersion_parameter
      }

      if(type == "competition") {
        self$n_quantiles = obj$n_quantiles
        self$continuous = obj$continuous
      }

    },
    setup_species_parameters = function(obj, type, hybrid) {
      self[[paste0(type, "_func")]] = private$set_environment(obj$func)
      self[[paste0(type, "_formula")]] = obj$formula
      self$register_buffer(paste0("par_", type, "_upper"), torch::torch_tensor(get_par_boundary(obj, type, upper = TRUE)))
      self$register_buffer(paste0("par_", type, "_lower"), torch::torch_tensor(get_par_boundary(obj, type, upper = FALSE)))

      self[[paste0("par_", type, "_optimized")]] = obj$optimizeSpecies

      # to keep things simple, in case of hybrid, we will just ignore the species parameters for now....TODO!

      # if null, random initialisierung
      # else recalclate from init values the required values
      if(is.null(obj$initSpecies)) {
        self[[paste0("par_", type)]] = init_species_parameters(type, self$N_species)
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

