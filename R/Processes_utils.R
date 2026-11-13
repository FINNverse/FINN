#' Define a demographic process for FINN
#'
#' @description
#' Configures one demographic process (growth, mortality, regeneration or
#' competition) for use as `mortality_process`, `growth_process`,
#' `regeneration_process`, or `competition_process` in [FINN::finn()], either as a
#' fully mechanistic process with an explicit, interpretable functional form, or
#' -- if `hidden` is supplied -- as the first ("Level 1") level of hybridization
#' described by Pichler & Käber (2026): only the process' environmental-response
#' function is replaced by a small feed-forward neural network, while the rest of
#' the process equation remains mechanistic. To replace the entire process
#' equation with a neural network ("Level 2"), use [FINN::createHybrid()] instead.
#'
#' @details
#' Each demographic process in FINN is the product of (i) a process equation
#' `func` that operates on the cohort state (dbh, number of trees, available
#' light, species) and species-specific process parameters, and (ii) a
#' species- and process-specific environmental-response function that maps
#' site-level environmental predictors to a scalar effect on the process (see
#' [FINN::finn()] for the underlying equations). `createProcess()` configures
#' both parts:
#'
#' * `func` implements the mechanistic process equation itself. The package's
#'   default process functions ([FINN::growth], [FINN::mortality],
#'   [FINN::regeneration], [FINN::competition]) reproduce the equations described
#'   by Pichler & Käber (2026); a custom function with the same arguments can be
#'   passed instead to use a different functional form while keeping the process
#'   embedded in, and jointly calibrated with, the rest of the model.
#' * The environmental-response function is, by default (`hidden = NULL`), a
#'   linear/logistic niche function with one coefficient per environmental
#'   covariate (named in `formula`) and per species, comparable to a classic
#'   species distribution model. Setting `hidden` to a vector of hidden-layer
#'   sizes (e.g. `c(25L)`) instead replaces this function with a feed-forward
#'   neural network, while `func` and the species-specific process parameters
#'   remain mechanistic -- the "Level 1" hybridization described in the paper.
#'
#' `formula` selects which columns of the `env` data (passed to [FINN::fit()] or
#' [FINN::predict.finn_class()]) enter the environmental-response function;
#' `initSpecies`/`initEnv` allow supplying custom starting values for the process
#' parameters/environmental-response model instead of the package's random
#' initialization; `optimizeSpecies`/`optimizeEnv` control whether these
#' parameters are estimated during [FINN::fit()] or held fixed at their initial
#' values; and `upper`/`lower` set box constraints (on the natural
#' process-parameter scale) within which the species-specific process parameters
#' are constrained during optimization.
#'
#' @param formula (`formula`)\cr Environmental predictors for this process'
#'   environmental-response function, evaluated against the `env` data. Default
#'   `NULL` uses `~.` (all available environmental covariates).
#' @param func (`function`)\cr Mechanistic process equation, e.g. [FINN::growth],
#'   [FINN::mortality], [FINN::regeneration], or [FINN::competition], or a custom
#'   function with the same arguments. Required.
#' @param initSpecies (`matrix`)\cr Optional custom initial values for the
#'   species-specific process parameters. Default `NULL` draws random initial
#'   values.
#' @param initEnv (`list`)\cr Optional custom initial values for the parameters of
#'   the environmental-response network/regression. Default `NULL` draws random
#'   initial values.
#' @param hidden (`integer()`)\cr Hidden-layer sizes of the feed-forward neural
#'   network that replaces the environmental-response function (Level 1 hybrid).
#'   Default `NULL` keeps the environmental response mechanistic (linear/logistic).
#' @param optimizeSpecies (`logical(1)`)\cr Should the species-specific process
#'   parameters be estimated during [FINN::fit()]? Default `FALSE`.
#' @param optimizeEnv (`logical(1)`)\cr Should the parameters of the
#'   environmental-response function be estimated during [FINN::fit()]? Default
#'   `TRUE`.
#' @param inputNN (`integer(1)`)\cr Input dimension of the environmental-response
#'   network. Default `NULL` infers it from `formula`/`env`.
#' @param outputNN (`integer(1)`)\cr Output dimension of the environmental-response
#'   network. Default `NULL` uses the number of species.
#' @param dispersion_parameter (`numeric(1)`)\cr Initial dispersion parameter of
#'   the negative binomial distribution used to sample recruits (the size
#'   `theta`), or its fixed value with `optimizeTheta = FALSE`. Only used when
#'   this process is the regeneration process. See [recruitment].
#' @param optimizeTheta (`logical(1)`)\cr Should the recruit dispersion
#'   `par_theta_recruits` be estimated during [FINN::fit()]? Default `TRUE`.
#'   With `FALSE` it stays at `dispersion_parameter` for the whole fit (a buffer,
#'   not a parameter), in the recruit draw and in the regeneration likelihood
#'   alike; `theta_floor` still applies. Only used when this process is the
#'   regeneration process.
#'
#'   Use `FALSE` when the recruitment level is to be calibrated: a trained
#'   `theta` falls on sparse recruitment data and leaves the level biased (see
#'   `theta_floor`). `dispersion_parameter = 1` is FINN's start value and held
#'   in tests on the German national forest inventory; a higher value such as
#'   `50` brings the draw close to a Poisson draw and gave the lowest holdout
#'   error there. See [recruitment].
#' @param cap,prior,prior_lambda The Beverton-Holt cap of
#'   [FINN::regeneration_capped()]: `cap = "all"` (default, one cap shared by all
#'   species) or `"species"`; `prior = "max_obs_recruits"` (default, the largest
#'   observed recruitment rate in the fitted data, recruits per ha and year) or a
#'   positive number; `prior_lambda` (default `0.1`) the strength of the pull of
#'   `log K` towards `log prior`. See [FINN::regeneration_capped()]. Only used
#'   when this process is the regeneration process and its `func` uses the cap.
#' @param NN (`nn_module`)\cr Optional custom `torch` module overriding the default
#'   environmental-response network architecture.
#' @param upper (`numeric()`)\cr Upper boundaries (natural scale) for the
#'   species-specific process parameters. For the competition process the
#'   defaults are `c(log(2.3), 2.0)`: the height parameter is capped so that no
#'   species grows taller than 130 m (see [FINN::height()]).
#' @param lower (`numeric()`)\cr Lower boundaries (natural scale) for the
#'   species-specific process parameters.
#' @param dropout (`numeric(1)`)\cr Dropout rate of the environmental-response
#'   network. Ignored unless `hidden` is set.
#' @param sample_regeneration (`logical(1)`)\cr Should recruits actually be sampled
#'   from the negative binomial regeneration distribution (`TRUE`), or only its
#'   expectation used (`FALSE`)? Only used when this process is the regeneration
#'   process. With `FALSE` the expected number of recruits per patch is rounded,
#'   so wherever it is below 0.5 no recruit enters at all: a silent zero for
#'   sparse species on small patches. See [recruitment].
#' @param n_quantiles (`integer(1)`)\cr Number of height classes used to discretize
#'   cohorts when computing shading/light availability. Only used when this
#'   process is the competition process and `continuous = FALSE`.
#' @param light_mask (`logical(1)`)\cr Take the height quantiles of the binned
#'   competition over the LIVE cohorts only. The cohort dimension is padded to the
#'   widest site x patch in the batch and an empty slot has height 0.001, so with
#'   enough padding every live cohort is searchsorted into the top bin, where
#'   nothing shades anything: measured on a 251 m2/ha stand of ten live cohorts,
#'   100 padded slots put all ten at light 1.00. Only the quantile INPUT is masked,
#'   so every cohort keeps its real height as the bin query. `TRUE` is the default
#'   since 0.2.0, because the padded width is set by the widest site in the BATCH:
#'   the same stand of ten cohorts gets light 1.00 throughout beside a neighbour
#'   holding a hundred, and its correct gradient down to 0.159 beside a neighbour
#'   holding ten. Batch composition is not part of the model, so the unmasked
#'   quantile made the fitted objective depend on `batchsize`. `FALSE` reproduces
#'   results published before 0.2.0. Ignored when `continuous = TRUE`, which
#'   compares pairs and is correct at any padding.
#' @param continuous (`logical(1)`)\cr Compute shading continuously for every pair
#'   of cohorts (`TRUE`) instead of binning cohorts into `n_quantiles` height
#'   classes (`FALSE`, default). Only used when this process is the competition
#'   process.
#' @param mortality_draw (`character(1)`)\cr How the number of dead trees is drawn
#'   from the cohort's mortality probability. Only used when this process is the
#'   mortality process. `"binomial"` (default since 0.2.0) is the exact binomial:
#'   each tree in the cohort gets its own reparameterised Bernoulli trial, which is
#'   integer-valued, matches the binomial in both moments -- measured 12.011 / 8.387
#'   against an exact 12.000 / 8.400 at `n = 40`, `p = 0.3` -- and keeps a gradient
#'   path to `p` through a straight-through sigmoid. It expands each cell to the
#'   largest cohort in the batch, so its memory grows with the largest `trees` value
#'   in the data.
#'   `"poisson"` draws a Bernoulli where
#'   the cohort holds one tree or fewer and a Poisson capped at the cohort size
#'   otherwise: integer-valued, exact in expectation for small `n*p`, and no
#'   pathwise gradient -- which costs nothing now that `fit(rb_gradient = TRUE)`
#'   takes the gradient through the expectation anyway. The Bernoulli branch
#'   matters because capping a Poisson at a one-tree cohort delivers only 88 % of
#'   the calibrated probability at the `p = 0.26` of a suppressed stem.
#'   `"gamma"` (the behaviour of every published FINN result) matches a Gamma to the binomial's first two moments and `rsample()`s
#'   it, which gives a pathwise gradient -- the hazard a discrete draw does not
#'   cover -- but the forward pass rounds that draw straight-through, so the
#'   gradient flows through a quantity the simulation discards and the rounded
#'   count delivers a measured 0.68-0.72 of the calibrated deaths for one-tree
#'   cohorts at annual probabilities of 0.01-0.10.
#'   See [FINN::binomial_from_poisson()] for the measured effect on the simulated
#'   stand.
#' @param theta_floor (`numeric(1)`)\cr Lower bound on the negative-binomial
#'   recruit dispersion `par_theta_recruits`, applied wherever it is read: the
#'   recruit draw and the regeneration likelihood. Default `0.5` since 0.2.0;
#'   `0` or `NULL` leaves it unbounded, as in 0.1.0. Only used when this object
#'   is the regeneration process.
#'
#'   theta is trained, and it collapses, because lowering it widens the
#'   regeneration likelihood and so cheapens that term without improving the
#'   predicted mean. The same parameter sets the variance of the *simulated*
#'   recruit draw, `mu + mu^2/theta`, which then turns bursty: on the German NFI
#'   it reached 0.035 by epoch 50 and the simulated stand overshot observed
#'   recruitment 133-fold, with basal-area and stem losses flat until epoch 25
#'   and running away after. The Beverton-Holt cap does not prevent this -- it
#'   bounds the mean, and the variance around a capped mean is still
#'   `K + K^2/theta`.
#'
#'   The floor guarantees stability only. It does not calibrate the recruitment
#'   level: on the German national forest inventory, out of sample on a
#'   spatially blocked holdout, a floored `theta` left the predicted recruitment
#'   level about 12 times too high, while a fixed `theta` (`optimizeTheta =
#'   FALSE`) or a Poisson regeneration likelihood (`fit(loss_family = )`) gave 2
#'   to 2.8 times. The floor is a guard rail for users who set nothing, not the
#'   recommended setting for calibrating the recruitment level; use
#'   `optimizeTheta = FALSE` or `loss_family` for that. See [recruitment].
#' @param mortality_floor,mortality_ceiling (`numeric(1)`)\cr The interval the
#'   per-cohort mortality probability is clamped to **per simulation step**,
#'   before the death draw. Only used when this process is the mortality process;
#'   the per-tree mortality loss is not bounded by either.
#'
#'   The defaults are `0.001` and `0.999` (0.1 % and 99.9 % per step). Until
#'   0.2.0 the interval was hard-coded to `[0.01, 0.99]`, which is why a
#'   simulated low-mortality stand -- a mature unthinned stand, a slow-growing
#'   species, an equilibrium run -- lost stems the observations did not: at a
#'   1 % floor it loses 1 % of every cohort per step however low the fitted rate.
#'
#'   The floor exists because `binomial_from_gamma()`'s shape `n*p/(1-p)`
#'   degenerates as `p -> 0`. Measured with that draw, it stays finite with an
#'   unbiased mean gradient down to `1e-4`, but single-draw gradients grow
#'   heavy-tailed as `n*p` shrinks -- the largest per-tree gradient is about ten
#'   times the old default's at `1e-3` and about a hundred times at `1e-4` --
#'   so `1e-3` is the lowest tested value and is the new default.
#' @param reg_floor (`numeric`)\cr A small recruitment that every species gets
#'   in every patch and time step, whatever the light and the environment, in
#'   trees per ha and time step. Only used when this process is the
#'   regeneration process ([FINN::regeneration()], [FINN::regeneration_capped()]
#'   and [FINN::regeneration_adult()] add it to their expected recruitment:
#'   `mean = regP * exp(env) + reg_floor`). Default `1e-3`, which keeps the
#'   recruitment likelihood finite where the light response closes and adds
#'   next to nothing. One value for all species or one per species. Before
#'   0.2.0 the floor was hard-coded to `0.2`, which with many species adds a
#'   large recruitment that does not depend on anything (65 species x 0.2 = 13
#'   trees per ha and year); models saved before then keep that value when
#'   reloaded. `finn(reg_floor = )`, the argument of earlier versions, is
#'   still accepted with a warning and overrides this value.
#' @param growth_scale (`numeric(1)` or `NULL`)\cr Overrides the `growth_scale`
#'   default of the process `func`. Only used when this process is the mortality
#'   process and `func` declares that argument. `NULL` (the default) keeps the
#'   function's own default, which for [FINN::mortality()] is `0.01` -- growth
#'   enters the logit in %/yr. `growth_scale = 1` feeds it the raw annual relative
#'   rate, the behaviour of 0.1.0. The value is set as the function's formal rather
#'   than passed down the forward pass, so a custom mortality `func` without the
#'   argument keeps working unchanged.
#' @param regeneration_sampler (`character(1)`)\cr How recruits are drawn from the
#'   negative binomial. Only used when this process is the regeneration process and
#'   `sample_regeneration = TRUE`. `"refined"` (default since 0.2.0) samples
#'   exactly below `lambda = 2` and drops both guards below, neither of which
#'   applies to `torch_poisson`; see [FINN::rnbinom_torch_refined()].
#'   `"gaussian"` (the behaviour of every
#'   published FINN result) offsets the Gamma draw by `+0.001`, clamps the
#'   probability at `1e-4` and pushes the result through
#'   `lambda + eps*sqrt(lambda)`. The offset and the clamp guard the Gaussian
#'   branch against a negative draw and against a NaN gradient through `log(0)`;
#'   together they are a recruitment FLOOR of 0.0011 per cell-year, which at a FINN
#'   patch of 0.0168 ha is 0.066 recruits ha^-1 yr^-1 per species -- over a pool of
#'   a hundred-odd species, as large as a temperate inventory's whole ingrowth.
#' @param custom_parameters (`named list` or `NULL`)\cr Extra parameters that a
#'   custom `func` needs, declared here so they are registered on the model by
#'   [FINN::finn()], optimised jointly by [FINN::fit()] (when trainable), and
#'   round-tripped by `torch_save()`/`torch_load()`. Each element is either a
#'   numeric vector (a trainable parameter with that starting value) or a list with
#'   `init` (numeric, length 1 or `N_species`) and optional `trainable`
#'   (`logical`, default `TRUE`; `FALSE` makes it a fixed buffer). The custom
#'   function accesses each as `self[[name]]`. Parameterise in the *unconstrained*
#'   space and apply any link inside the function (e.g. store `log K` and use
#'   `exp(self$logK)`), mirroring FINN's own `par_*_unconstrained` convention.
#'
#' @return A list of class `"process"` containing the process definition and
#'   associated parameters, to be passed as `mortality_process`, `growth_process`,
#'   `regeneration_process`, or `competition_process` to [FINN::finn()].
#'
#' @references
#' Pichler, M., & Käber, Y. (2026). Inferring processes within dynamic forest
#' models using hybrid modelling. \emph{Methods in Ecology and Evolution}.
#' \doi{10.1111/2041-210x.70347}
#'
#' @seealso [FINN::createHybrid()], [FINN::finn()], [recruitment] for the
#'   recruitment options and how to combine them
#'
#' @examples
#' growth_process <- createProcess(formula = ~temperature + precipitation, func = growth)
#'
#' # the regeneration process with a smaller recruitment floor and a cap per species
#' regeneration_process <- createProcess(~temperature, func = regeneration_capped,
#'                                       reg_floor = 1e-4, cap = "species")
#'
#' @export
createProcess = function(formula = NULL, func, initSpecies = NULL, initEnv = NULL, hidden = NULL, optimizeSpecies = FALSE, optimizeEnv = TRUE, inputNN = NULL, outputNN = NULL, dispersion_parameter = 1.0, NN = NULL, upper = NULL, lower = NULL, dropout = 0.0, sample_regeneration = TRUE, n_quantiles = 10L, continuous = FALSE, light_mask = TRUE, custom_parameters = NULL, mortality_draw = c("binomial", "bernoulli", "poisson", "gamma"), regeneration_sampler = c("refined", "gaussian"), growth_scale = NULL,
                         cap = c("all", "species"), prior = "max_obs_recruits", prior_lambda = 0.1,
                         theta_floor = 0.5, optimizeTheta = TRUE,
                         mortality_floor = 0.001, mortality_ceiling = 0.999,
                         reg_floor = 1e-3) {
  out = list()
  out$reg_floor = check_reg_floor(reg_floor)
  out$optimizeTheta = validate_optimize_theta(optimizeTheta)
  out$mortality_draw = match.arg(mortality_draw)
  out$mortality_floor = check_mortality_limit(mortality_floor, "mortality_floor")
  out$mortality_ceiling = check_mortality_limit(mortality_ceiling, "mortality_ceiling")
  if(out$mortality_floor >= out$mortality_ceiling)
    stop("`mortality_floor` must be below `mortality_ceiling`.", call. = FALSE)
  out$regeneration_sampler = match.arg(regeneration_sampler)
  # `growth` enters the logit as a bare divisor, so a bad value is silent rather than loud:
  # logical NA coerces to TRUE and divides by 1, which restores the pre-0.2.0 objective, and
  # 0 saturates mortality to {0, 1} across the stand without a NaN reaching a live cohort.
  if(!is.null(growth_scale) && !(is.numeric(growth_scale) && length(growth_scale) == 1L &&
                                 is.finite(growth_scale) && growth_scale > 0))
    stop("`growth_scale` must be a single positive finite number, or NULL.", call. = FALSE)
  out$cap = match.arg(cap)
  out$prior = validate_cap_prior(prior)
  out$prior_lambda = validate_prior_lambda(prior_lambda)
  out$theta_floor = validate_theta_floor(theta_floor)
  # the competition kernel gates on isTRUE(), so `light_mask = 1` would be accepted
  # here and then silently leave the quantiles unmasked
  if(!is.logical(light_mask) || length(light_mask) != 1L || is.na(light_mask))
    stop("`light_mask` must be TRUE or FALSE.", call. = FALSE)
  out$growth_scale = growth_scale

  # ---- custom (extra) process parameters -------------------------------------
  # A named list declaring extra parameters a custom `func` needs. Each element is
  # either a numeric vector (-> trainable) or a list with `init` (numeric, length 1
  # or N_species) and optional `trainable` (logical, default TRUE). finn() registers
  # each on the model (see setup_species_parameters), so fit() optimises the
  # trainable ones and torch_save/torch_load round-trips them; the custom func reads
  # it as `self[[name]]`. Parameterise in the unconstrained space and apply any link
  # inside the func (e.g. store log K, use exp(self$logK)).
  if(!is.null(custom_parameters)) {
    if(is.null(names(custom_parameters)) || any(names(custom_parameters) == ""))
      stop("`custom_parameters` must be a *named* list (one entry per parameter).", call. = FALSE)
    custom_parameters = lapply(custom_parameters, function(spec) {
      if(is.numeric(spec)) spec = list(init = spec)
      if(is.null(spec$init) || !is.numeric(spec$init))
        stop("each custom parameter needs a numeric `init`.", call. = FALSE)
      spec$trainable = if(is.null(spec$trainable)) TRUE else isTRUE(spec$trainable)
      spec
    })
  }
  out$custom_parameters = custom_parameters

  if(!is.null(formula)){
    mf = match.call()
    m = match("formula", names(mf))
    if(inherits(mf[m]$formula, "name")) mf[m]$formula = eval(mf[m]$formula, envir = parent.env(environment()))
    formula = stats::as.formula(mf[m]$formula)
    #X = stats::model.matrix(formula, data)
  } else {
    formula = stats::as.formula("~.")
    #X = stats::model.matrix(formula, data)
  }

  out$formula = formula
  out$func = func
  out$custom = TRUE
  if(isNamespace(environment(func))) out$custom = getNamespaceName(environment(func)) != "FINN"
  out$initSpecies = initSpecies
  out$optimizeSpecies = optimizeSpecies
  out$optimizeEnv = optimizeEnv
  out$inputNN = inputNN
  out$outputNN = outputNN
  out$dispersion_parameter = dispersion_parameter
  out$sample_regeneration = sample_regeneration
  out$NN = NN
  out$upper = upper
  out$lower = lower
  out$dropout = dropout
  out$n_quantiles = n_quantiles
  out$continuous = continuous
  out$light_mask = light_mask
  if(!is.null(initEnv)) {
    if(!is.list(initEnv)) initEnv = list(initEnv) # must be a list!
  }
  out$initEnv = initEnv
  out$hidden = hidden

  if(out$custom) cli::cli_text("Custom function detected...")
  if(!is.null(out$hidden)) cli::cli_text("Neural Network detected...")
  class(out) = "process"
  return(out)
}


#' Define a hybrid (deep-neural-network) demographic process for FINN
#'
#' @description
#' Configures a process (growth, mortality, or regeneration) in which the entire
#' process equation -- not just its environmental-response function -- is
#' replaced by a deep neural network (DNN), for use as `growth_process`,
#' `mortality_process`, or `regeneration_process` in [FINN::finn()]. This is the
#' second ("Level 2") level of hybridization described by Pichler & Käber (2026);
#' to replace only the environmental-response function while keeping the rest of
#' the process mechanistic ("Level 1"), use [FINN::createProcess()] instead.
#' `competition_process` must always be created with [FINN::createProcess()], as
#' the competition (light availability) process does not support full
#' replacement by a DNN.
#'
#' @details
#' The network receives the same cohort- and site-level information that the
#' corresponding mechanistic process equation would use -- diameter at breast
#' height, number of trees, available light, species identity and the site's
#' environmental predictors (plus the growth rate, for the mortality process) --
#' and predicts the process output directly. Species identity is passed through a
#' learned embedding layer rather than treated as a categorical covariate with
#' per-species coefficients, so the network can learn low-dimensional "contrasts"
#' between species or plant functional types (PFTs). An inverse-link function
#' appropriate to each process is applied to the network's output: a sigmoid for
#' mortality (a per-tree death probability, as for [FINN::mortality]), and an
#' exponential for growth and regeneration (always-positive rates; for
#' regeneration this is the mean of the negative binomial distribution from which
#' recruits are drawn, as for [FINN::regeneration]).
#'
#' Two network architectures are available, selected with `transformer`. The
#' default (`transformer = TRUE`) is a small transformer encoder (`encoder_layers`
#' layers, embedding dimension `emb_dim`, feed-forward dimension
#' `dim_feedforward`) that embeds each cohort, its species and the site's
#' environment and attends across the cohorts of a patch. Setting
#' `transformer = FALSE` instead uses a feed-forward network with hidden layers
#' `hidden` (default two layers of 50 units each) and `dropout`, matching the
#' architecture used by Pichler & Käber (2026) for the Barro Colorado Island case
#' study (where dropout was set to 10%).
#'
#' As with [FINN::createProcess()], hybrid processes are calibrated jointly,
#' end-to-end, with the remaining mechanistic or hybrid processes via
#' [FINN::fit()], rather than pre-trained in isolation and plugged in afterwards.
#' `optimize` controls whether the network's weights are estimated during fitting
#' or kept fixed at their (random) initial values, and `dispersion_parameter`/
#' `sample_regeneration` have the same meaning as in [FINN::createProcess()] and
#' are only used when this object is the regeneration process.
#'
#' @param formula (`formula`)\cr Environmental predictors passed to the network,
#'   evaluated against the `env` data. Default `NULL` uses `~.` (all available
#'   environmental covariates).
#' @param optimize (`logical(1)`)\cr Should the network's weights be estimated
#'   during [FINN::fit()]? Default `TRUE`.
#' @param dispersion_parameter (`numeric(1)`)\cr Initial dispersion parameter of
#'   the negative binomial distribution used to sample recruits, or its fixed
#'   value with `optimizeTheta = FALSE`. Only used when this object is the
#'   regeneration process. See [recruitment].
#' @param optimizeTheta (`logical(1)`)\cr Should the recruit dispersion be
#'   estimated during [FINN::fit()]? Default `TRUE`; `FALSE` keeps it at
#'   `dispersion_parameter`, which is the choice when the recruitment level is
#'   to be calibrated; see [FINN::createProcess()] and [recruitment]. Only used
#'   when this object is the regeneration process.
#' @param cap,prior,prior_lambda Accepted for symmetry with
#'   [FINN::createProcess()], but not used in 0.2.0: the hybrid regeneration
#'   network has no Beverton-Holt cap. See [FINN::regeneration_capped()] for the
#'   cap of the mechanistic process.
#' @param NN (`nn_module`)\cr Optional custom `torch` module overriding the
#'   default network architecture entirely.
#' @param dropout (`numeric(1)`)\cr Dropout rate of the network.
#' @param encoder_layers (`integer(1)`)\cr Number of transformer encoder layers.
#'   Only used when `transformer = TRUE`.
#' @param hidden (`integer()`)\cr Hidden-layer sizes of the feed-forward network.
#'   Only used when `transformer = FALSE`.
#' @param sample_regeneration (`logical(1)`)\cr Should recruits actually be sampled
#'   from the negative binomial regeneration distribution (`TRUE`), or only its
#'   expectation used (`FALSE`)? Only used when this object is the regeneration
#'   process. With `FALSE` the expected number of recruits per patch is rounded,
#'   so wherever it is below 0.5 no recruit enters at all. See [recruitment].
#' @param mortality_draw (`character(1)`)\cr How the number of dead trees is drawn,
#'   `"binomial"` (default since 0.2.0), `"poisson"` or `"gamma"`; see
#'   [FINN::createProcess()]. Only used
#'   when this object is the mortality process.
#' @param mortality_floor,mortality_ceiling (`numeric(1)`)\cr The interval the
#'   mortality probability is clamped to per simulation step before the death
#'   draw, defaults `0.001` and `0.999`; see [FINN::createProcess()]. Only used
#'   when this object is the mortality process.
#' @param regeneration_sampler (`character(1)`)\cr How recruits are drawn,
#'   `"refined"` (default since 0.2.0) or `"gaussian"`; see [FINN::createProcess()]. Only used
#'   when this object is the regeneration process.
#' @param theta_floor (`numeric(1)`)\cr Lower bound on the recruit dispersion
#'   `par_theta_recruits`, default `0.5` since 0.2.0, `0` or `NULL` to leave it
#'   unbounded as in 0.1.0; see [FINN::createProcess()]. Only used when this
#'   object is the regeneration process. The floor guarantees stability only;
#'   it leaves the recruitment level biased upward and is not the recommended
#'   setting for calibrating it (use `optimizeTheta = FALSE` or
#'   `fit(loss_family = )`). See [recruitment].
#' @param transformer (`logical(1)`)\cr Use a transformer encoder architecture
#'   (`TRUE`, default) or a feed-forward network (`FALSE`).
#' @param emb_dim (`integer(1)`)\cr Embedding dimension for the species/cohort/
#'   environment features.
#' @param dim_feedforward (`integer(1)`)\cr Dimension of the feed-forward
#'   sublayer within each transformer encoder layer. Only used when
#'   `transformer = TRUE`.
#'
#' @return A list of class `"hybrid"` containing the process definition and
#'   associated parameters, to be passed as `growth_process`, `mortality_process`,
#'   or `regeneration_process` to [FINN::finn()].
#'
#' @references
#' Pichler, M., & Käber, Y. (2026). Inferring processes within dynamic forest
#' models using hybrid modelling. \emph{Methods in Ecology and Evolution}.
#' \doi{10.1111/2041-210x.70347}
#'
#' @seealso [FINN::createProcess()], [FINN::finn()], [recruitment]
#'
#' @examples
#' growth_process <- createHybrid(formula = ~temperature + precipitation)
#'
#' @export
createHybrid = function(formula = NULL, optimize = TRUE, dispersion_parameter = 1.0, NN = NULL, dropout = 0.3, encoder_layers = 1L, hidden = c(50L, 50L), sample_regeneration = TRUE, transformer = TRUE, emb_dim = 20L, dim_feedforward = 256L, mortality_draw = c("binomial", "bernoulli", "poisson", "gamma"), regeneration_sampler = c("refined", "gaussian"),
                        cap = c("all", "species"), prior = "max_obs_recruits", prior_lambda = 0.1,
                        theta_floor = 0.5, optimizeTheta = TRUE,
                        mortality_floor = 0.001, mortality_ceiling = 0.999) {
  out = list()
  out$optimizeTheta = validate_optimize_theta(optimizeTheta)
  out$mortality_draw = match.arg(mortality_draw)
  out$mortality_floor = check_mortality_limit(mortality_floor, "mortality_floor")
  out$mortality_ceiling = check_mortality_limit(mortality_ceiling, "mortality_ceiling")
  if(out$mortality_floor >= out$mortality_ceiling)
    stop("`mortality_floor` must be below `mortality_ceiling`.", call. = FALSE)
  out$regeneration_sampler = match.arg(regeneration_sampler)
  if(!is.null(formula)){
    mf = match.call()
    m = match("formula", names(mf))
    if(inherits(mf[m]$formula, "name")) mf[m]$formula = eval(mf[m]$formula, envir = parent.env(environment()))
    formula = stats::as.formula(mf[m]$formula)
    #X = stats::model.matrix(formula, data)
  } else {
    formula = stats::as.formula("~.")
    #X = stats::model.matrix(formula, data)
  }

  out$formula = formula
  out$custom = TRUE
  out$optimizeSpecies = FALSE
  out$optimizeEnv = optimize
  out$cap = match.arg(cap)
  out$prior = validate_cap_prior(prior)
  out$prior_lambda = validate_prior_lambda(prior_lambda)
  out$theta_floor = validate_theta_floor(theta_floor)
  out$dispersion_parameter = dispersion_parameter
  out$sample_regeneration = sample_regeneration
  out$NN = NN
  out$dropout = dropout
  out$encoder_layers = encoder_layers
  out$hidden = hidden
  out$transformer = transformer
  out$emb_dim = emb_dim
  out$dim_feedforward = dim_feedforward
  class(out) = "hybrid"
  return(out)
}

#' Extract Environmental Data for a Process
#'
#' This function extracts and transforms environmental data according to a specified process object. The environmental data is formatted into an array suitable for input into the simulation model.
#'
#' @param formula formula object `createProcess`.
#' @param env A data frame containing environmental data, which must include columns "siteID" and "year".
#'
#' @return An array of environmental data formatted according to the process object.
#'
#' @noRd
extract_env = function(formula, env) {
  process_env = stats::model.matrix(formula, env[,-c("siteID", "year")])
  env_names = colnames(process_env)
  env_array = climateDF2array(climate_dt =  cbind(env[,c("siteID", "year")],process_env), env_vars = env_names)
  return(env_array)
}

#' Compute per-predictor environmental scaling
#'
#' Computes the centering and scaling constants of every numeric environmental
#' predictor in \code{env} (all columns except the keys \code{siteID} and
#' \code{year}), so they can be re-applied unchanged to new data at prediction
#' time. Every supported mode is an affine map \code{(x - center) / scale}, so the
#' stored \code{(variable, center, scale)} table is applied by
#' \code{\link{apply_env_scaling}} and inverted by ALE unchanged, whatever the
#' mode. A predictor with (near-)zero standard deviation is given \code{scale = 1}
#' (centred only) to avoid division by zero.
#'
#' @param env A \code{data.frame}/\code{data.table} with \code{siteID}, \code{year}
#'   and the environmental predictor columns.
#' @param spec how to scale each predictor. One of:
#'   * `TRUE` -- z-standardise every predictor (`center = mean`, `scale = sd`);
#'     the default and the historical behaviour.
#'   * `FALSE` -- handled by the caller (no scaling table computed).
#'   * a length-1 string applied to all predictors, or a vector/list with one
#'     entry per predictor. Entries may be:
#'       - `"auto"`     z-standardise (mean / sd),
#'       - `"identity"` (alias `"none"`) leave unchanged (center 0, scale 1),
#'       - `"0to1"`     min-max to `[0, 1]` (center min, scale range),
#'       - a **function** `f(x)` returning `list(center=, scale=)` -- a custom
#'         affine scaler (e.g. robust: `function(x) list(center=median(x),
#'         scale=IQR(x))`). Non-affine transforms are intentionally not supported
#'         because ALE needs an invertible mapping.
#'   Unnamed vectors/lists are matched to predictors by position; named ones by
#'   variable name (predictors not named default to `"auto"`).
#' @return A \code{data.frame} with columns \code{variable}, \code{center},
#'   \code{scale}; or \code{NULL} if there are no numeric predictors. This is the
#'   object stored on a fitted model as \code{model$env_scaling} when it is fit
#'   with \code{env_autoscale = TRUE}.
#' @seealso \code{\link{apply_env_scaling}}
#' @export
compute_env_scaling = function(env, spec = TRUE) {
  cols = setdiff(colnames(env), c("siteID", "year"))
  cols = cols[vapply(cols, function(cc) is.numeric(env[[cc]]), logical(1))]
  if (!length(cols)) return(NULL)

  # a named entry that matches no predictor is a typo, and silently falling back to
  # "auto" for the predictor the user meant to name is the worst possible answer
  if (!is.logical(spec) && !is.null(names(spec))) {
    unknown = setdiff(names(spec)[nzchar(names(spec))], cols)
    if (length(unknown))
      stop(sprintf("`env_autoscale` names '%s', which is not an environmental predictor in `env` (%s).",
                   paste(unknown, collapse = "', '"), paste(cols, collapse = ", ")), call. = FALSE)
  }

  resolve_mode = function(col, i) {
    if (is.logical(spec)) return(if (isTRUE(spec)) "auto" else "identity")
    if (!is.null(names(spec))) {                      # named -> match by variable
      if (col %in% names(spec)) return(spec[[col]])
      return("auto")
    }
    if (length(spec) == 1L) return(spec[[1L]])        # one mode for all
    if (i <= length(spec)) return(spec[[i]])          # positional
    "auto"
  }
  affine = function(x, mode) {
    if (is.function(mode)) {
      out = mode(x)
      if (!is.list(out) || !all(c("center", "scale") %in% names(out)))
        stop("A scaling function must return list(center=, scale=).", call. = FALSE)
      return(list(center = out$center, scale = out$scale))
    }
    switch(as.character(mode),
      auto     = list(center = mean(x, na.rm = TRUE), scale = stats::sd(x, na.rm = TRUE)),
      identity = ,
      none     = list(center = 0, scale = 1),
      `0to1`   = list(center = min(x, na.rm = TRUE),
                      scale  = diff(range(x, na.rm = TRUE))),
      stop("Unknown env scaling mode '", mode,
           "'. Use 'auto', 'identity'/'none', '0to1', or a function.", call. = FALSE))
  }
  rows = lapply(seq_along(cols), function(i) {
    cs = affine(env[[cols[i]]], resolve_mode(cols[i], i))
    sc = cs$scale
    if (!is.finite(sc) || abs(sc) < .Machine$double.eps) sc = 1   # constant -> centre only
    data.frame(variable = cols[i], center = as.numeric(cs$center),
               scale = as.numeric(sc), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Apply stored z-standardization to environmental predictors
#'
#' Re-applies the constants learned by \code{compute_env_scaling()} to \code{env}.
#' The transformation uses the \emph{stored} mean/sd only and never recomputes them
#' from \code{env}, so calibration and prediction use an identical transformation
#' (the usual pitfall of \code{scale()} inside a model formula is avoided).
#'
#' @param env A \code{data.frame}/\code{data.table} of raw environmental data.
#' @param scaling The \code{data.frame} returned by \code{compute_env_scaling()},
#'   or \code{NULL} (in which case \code{env} is returned unchanged). For a model
#'   fit with \code{env_autoscale = TRUE} this is \code{model$env_scaling}.
#' @return \code{env} with the scaled predictor columns, as a \code{data.table}.
#' @seealso \code{\link{compute_env_scaling}}, \code{\link{fit}}
#' @examples
#' # reproduce the standardization a model fit with env_autoscale = TRUE used:
#' env <- data.frame(siteID = 1:3, year = 1L, temp = c(4, 6, 8), prec = c(700, 850, 1000))
#' scaling <- compute_env_scaling(env)
#' apply_env_scaling(env, scaling)
#' @export
apply_env_scaling = function(env, scaling) {
  if (is.null(scaling)) return(env)
  env = data.table::copy(data.table::as.data.table(env))
  missing_vars = setdiff(scaling$variable, colnames(env))
  if (length(missing_vars))
    stop("env is missing environmental variable(s) the model was scaled on: ",
         paste(missing_vars, collapse = ", "))
  for (i in seq_len(nrow(scaling)))
    set(env, j = scaling$variable[i],
        value = (env[[scaling$variable[i]]] - scaling$center[i]) / scaling$scale[i])
  env
}

backward = function(value, upper, lower) {
  ratio <- t((t(value) - as.numeric(lower)) / (as.numeric(upper) - as.numeric(lower)))
  # `qlogis` maps a value sitting exactly on a bound to +/-Inf, and one outside
  # the bounds to NaN. Either silently poisons the whole model (every downstream
  # tensor becomes NaN and the simulation collapses to zero). Clamp strictly
  # inside (0, 1) so initialisation stays finite, and warn if an initial value
  # was at or beyond its [lower, upper] range so the user can widen the bounds.
  if (any(!is.na(ratio) & (ratio <= 0 | ratio >= 1))) {
    warning("an initial parameter value is at or outside its [lower, upper] bound; ",
            "clamping it to keep the model finite. Consider passing `lower`/`upper` ",
            "to createProcess() to widen the range.", call. = FALSE)
  }
  eps <- 1e-6
  ratio <- pmin(pmax(ratio, eps), 1 - eps)
  stats::qlogis(ratio)
}
forward =  function(value, upper, lower) torch::torch_sigmoid(value) * (upper - lower) + lower

init_species_parameters = function(type, N_species, growth_scale = 1){
  initLight = stats::runif(N_species, 0.45, 0.55)
  if(type == "competition") {
    init =
      cbind(initLight,
            stats::runif(N_species, min = 0.19, 0.21)) # 0.2 corresponds to 0 light at 50m2/ha basal area)
  }
  if(type == "growth") {
    init =
      cbind(initLight,
            stats::runif(N_species, min = 0.045, 0.055))
  }
  if(type == "mortality") {
    # the third column multiplies growth / growth_scale, so the range has to move with
    # growth_scale or the same draw means something 1/growth_scale times stronger: at the
    # 0.01 default an untouched +-0.2 puts the logit up to 8 units into saturation on the
    # FIA fixture, killing 93 to 99.9 % of a 25 cm stand in one year (FINNverse/FINN#39)
    init =
      cbind(stats::runif(N_species, min = -0.2, 0.2),
            stats::runif(N_species, min = -0.2, 0.2),
            stats::runif(N_species, min = -0.2, 0.2) * growth_scale)
  }
  if(type == "regeneration") {
    init = matrix(initLight, ncol = 1L)
  }
  return(init)
}


get_par_boundary = function(obj, type, upper = TRUE) {
  if(upper) {
    if(is.null(obj$upper)) {
      upper = switch(type,
                     #mortality = { c(0.93, 4.00) },
                     mortality = { c(10.0, 10.0, 10.0) },
                     growth = {c(0.99, 4.0)},
                     regeneration = { 0.99 },
                     # height parameter: height() approaches (exp(p) - 1) * 100 m as dbh grows, so
                     # log(2.3) caps trees at 130 m, the physical limit of tree height (Koch et al. 2004)
                     competition = { c(log(2.3), 2.0)})
    } else {
      upper = obj$upper
    }
    return(upper)
  } else {
    if(is.null(obj$lower)) {
      lower = switch(type,
                     mortality = { c(-10.0, -10.0, -10.0) },
                     #mortality = { c(0.011, 0.00) },
                     growth =    {c(0.01, 0.01)},
                     regeneration = { 0.01 },
                     competition = { c(0.3, 0.0)})
    } else {
      lower = obj$lower
    }
    return(lower)
  }
}

# The recruitment cap's prior is specified identically on createProcess() and
# createHybrid(); validating it in one place keeps the two from drifting apart.
# The additive recruitment floor: non-negative, finite, one value or one per
# species. Checked in createProcess() (length 1 or any, N_species unknown) and
# again when the process is attached to a model.
check_reg_floor = function(reg_floor, N_species = NULL) {
  ok = is.numeric(reg_floor) && all(is.finite(reg_floor)) && all(reg_floor >= 0) &&
    (if (is.null(N_species)) length(reg_floor) >= 1L else length(reg_floor) %in% c(1L, N_species))
  if (!ok)
    stop("`reg_floor` must be non-negative finite numbers, of length 1 or `N_species`.", call. = FALSE)
  reg_floor
}

validate_cap_prior = function(prior) {
  if(identical(prior, "max_obs_recruits")) return(prior)
  if(is.numeric(prior) && length(prior) == 1L && is.finite(prior) && prior > 0) return(as.numeric(prior))
  stop('`prior` must be "max_obs_recruits" or a single positive finite number ',
       '(recruits ha^-1 yr^-1).', call. = FALSE)
}

# A lower bound on the recruit dispersion, not an upper one: theta only ever
# collapses, because lowering it cheapens the regeneration likelihood.
validate_optimize_theta = function(optimizeTheta) {
  if(!is.logical(optimizeTheta) || length(optimizeTheta) != 1L || is.na(optimizeTheta))
    stop("`optimizeTheta` must be TRUE or FALSE.", call. = FALSE)
  optimizeTheta
}

validate_theta_floor = function(theta_floor) {
  if(is.null(theta_floor)) return(NULL)
  if(is.numeric(theta_floor) && length(theta_floor) == 1L &&
     is.finite(theta_floor) && theta_floor >= 0) {
    return(if(theta_floor > 0) as.numeric(theta_floor) else NULL)
  }
  stop("`theta_floor` must be a single finite number >= 0, or NULL ",
       "(0 and NULL both leave the recruit dispersion unbounded, as in 0.1.0).",
       call. = FALSE)
}

validate_prior_lambda = function(prior_lambda) {
  if(is.numeric(prior_lambda) && length(prior_lambda) == 1L &&
     is.finite(prior_lambda) && prior_lambda >= 0) return(as.numeric(prior_lambda))
  stop("`prior_lambda` must be a single finite number >= 0 (0 switches the prior off).",
       call. = FALSE)
}
