#' @title Recruitment in FINN: the options and how to choose them
#'
#' @description
#' Recruitment in FINN is set in several places: the regeneration process
#' function and its arguments in [FINN::createProcess()] (or
#' [FINN::createHybrid()]), two arguments of [FINN::finn()], and the likelihood
#' and recruitment mode of [FINN::fit()], [FINN::predict.finn_class()] and
#' [FINN::simulateForest()]. This page says what each option does, how they act
#' together, and which combination to use for which purpose. The individual help
#' pages carry the full argument descriptions.
#'
#' One simulated year of recruitment runs in three steps. The regeneration
#' process returns an EXPECTED recruitment per hectare for every site, patch and
#' species. This rate is turned into a whole number of recruits per patch, which
#' enter the stand as new cohorts at `finn(recruits_dbh = )`. During [FINN::fit()]
#' the expected rate (not the drawn recruits), summed over an observation
#' interval, is compared with the observed `reg` by the regeneration likelihood.
#'
#' @section The recruitment process:
#' The process is the `func` of the regeneration process passed to
#' `finn(regeneration_process = )`. All mechanistic variants share the same core:
#' a light response gated by the species parameter `parReg`, times the
#' environmental response `exp(env)`, plus the additive floor
#' `createProcess(reg_floor = )` of the regeneration process (default `1e-3`
#' stems per ha and year).
#'
#' * [FINN::regeneration()]: the core alone. The expected recruitment is not
#'   bounded; it grows with the environmental response.
#' * [FINN::regeneration_capped()]: the default since 0.2.0, used when
#'   `regeneration_process` is `NULL`. The core is capped by a Beverton-Holt
#'   term, `K m / (K + m)`. `createProcess(cap = "all")` (default) gives one cap
#'   `K` shared by all species, `cap = "species"` one per species. The cap carries
#'   a prior, `prior_lambda * (log K - log prior)^2`, with
#'   `prior = "max_obs_recruits"` (the largest observed recruitment per ha and
#'   year) and `prior_lambda = 0.1` by default. The cap bounds the MEAN
#'   recruitment only. It does not bound the variance of the recruit draw, which
#'   around a capped mean is still `K + K^2 / theta` (see the dispersion section
#'   below).
#' * `regeneration_saturation()`: the same function as `regeneration_capped()`
#'   under its older name, kept as a deprecated alias for one release.
#' * [FINN::regeneration_adult()]: the capped process with a seed source. The
#'   recruitment of a species is multiplied by `f(A) = A / (A + A_half)`, where
#'   `A` is the density of its own adults (stems above an adult threshold `x`)
#'   pooled over the patches of the site. A species without adults on a site
#'   recruits only at `reg_floor`. `x` and `A_half` are declared through
#'   `createProcess(custom_parameters = )`; see its help page.
#'
#' A regeneration process built with [FINN::createHybrid()] replaces the whole
#' equation by a network. In 0.2.0 it has no Beverton-Holt cap: `createHybrid()`
#' accepts `cap`, `prior` and `prior_lambda`, but the network does not use them.
#'
#' @section From expected recruitment to recruits:
#' The expected rate per ha times `patch_size` (in ha) is the expected number
#' of recruits per patch. `createProcess(sample_regeneration = )` decides how it
#' becomes a whole number:
#'
#' * `sample_regeneration = TRUE` (default): recruits are drawn from a negative
#'   binomial with that mean and size `theta` (the dispersion, next section),
#'   using `createProcess(regeneration_sampler = )`. `"refined"` (default since
#'   0.2.0) draws the Gamma-mixed Poisson exactly where its rate is below 2 and
#'   uses a Gaussian approximation above. `"gaussian"` reproduces results
#'   published before 0.2.0; its guards act as a recruitment floor of 0.0011 per
#'   cell and year (see [FINN::rnbinom_torch_refined()]).
#' * `sample_regeneration = FALSE`: the expected number is used and rounded.
#'   Wherever the expected number per patch is below 0.5 the rounding delivers
#'   **no recruit at all**. At small patches and low rates this is a silent zero
#'   for every sparse species, however well its rate is calibrated.
#'
#' @section The recruit dispersion theta:
#' `par_theta_recruits` is the negative-binomial size `theta`, one value per
#' species. It sets the variance of the recruit draw, `mu + mu^2 / theta`, and
#' the width of the negative-binomial regeneration likelihood. Three arguments of
#' [FINN::createProcess()] (and [FINN::createHybrid()]) control it:
#'
#' * `dispersion_parameter`: the starting value, or the fixed value when
#'   `theta` is not trained. Default `1`.
#' * `optimizeTheta`: `TRUE` (default) trains `theta`. `FALSE` keeps it at
#'   `dispersion_parameter` for the whole fit, in the draw and in the
#'   likelihood alike.
#' * `theta_floor`: a lower bound applied wherever `theta` is read, trained or
#'   not. Default `0.5` since 0.2.0; `0` or `NULL` removes it. A fixed
#'   `dispersion_parameter` below the floor is raised to the floor.
#'
#' With the default `fit(rb_gradient = TRUE)`, `theta` is set only by the
#' regeneration likelihood. On sparse, zero-heavy recruitment data the
#' likelihood optimum of `theta` is very small: a small `theta` widens the
#' likelihood and so lowers the loss without improving the predicted mean, and
#' the same small `theta` makes the simulated draws bursty. A free `theta` can
#' therefore run away (FINNverse/FINN#63).
#'
#' What the floor does and does not do: it guarantees that the fit stays
#' stable, but it does not calibrate the recruitment level. On the German
#' national forest inventory, out of sample on a spatially blocked holdout, a
#' floored `theta` left the predicted recruitment level about 12 times too high,
#' with an error 2.8 times that of predicting the mean. A fixed `theta` (1 or 50,
#' `optimizeTheta = FALSE`) or a Poisson regeneration likelihood gave a level 2
#' to 2.8 times too high and an error at the level of the mean, and a free
#' `theta` ran away by four orders of magnitude. The floor is a guard rail for
#' users who set nothing, not the recommended setting for calibrating the
#' recruitment level.
#'
#' @section The likelihood and the observation operator:
#' `fit(loss_family = c(regeneration = ))` chooses the regeneration likelihood.
#' `"nbinom"` (default) uses `theta` as above. `"poisson"` has no dispersion
#' parameter and so cannot collapse; `theta` then only shapes the recruit draw.
#'
#' `finn(recruit_obs_weight = w)` is an observation operator: it maps the
#' predicted and the observed recruitment per ha alike into the count space of
#' the sampling design before they are compared, so that the likelihood sees
#' the counts that were actually tallied. For an angle-count inventory
#' `w = g(recruits_dbh) / BAF`. For a fixed-radius small-tree plot that counts
#' the stems in a diameter band below the recruitment threshold,
#'
#' `w = plot area (ha) x residence time in the band (years) / observation period (years)`,
#'
#' where the residence time is the band width divided by the diameter growth of
#' the stems in it. One tallied stem then stands for `1 / w` recruits per ha over
#' the observation period. The residence time is an assumption and should be
#' estimated from the data where possible: `w` scales every absolute recruitment
#' number the fit returns, though not the ranking between model configurations.
#'
#' In 0.2.0 `recruit_obs_weight` is applied **only** under `"nbinom"`. With
#' `loss_family["regeneration"] = "poisson"` (or `"mse"`, `"gaussian"`) the term
#' scores the per-ha values as they are and `recruit_obs_weight` has no effect,
#' without a warning; the Poisson variance then is that of a per-ha value, not of
#' the tallied count (FINNverse/FINN#65). The Poisson term adds
#' `fit(trees_offset = )` (default `0.01`) to the prediction.
#'
#' @section Feeding observed recruitment:
#' `fit(loss_recruit = "observed")` drives the simulated stand with the observed
#' recruitment instead of the regeneration process. The observed `reg` divided
#' by `period_length` is the expected recruitment per ha and year in every year
#' of that interval; times `patch_size` it is drawn per patch by the recruit
#' sampler, as above. Growth, mortality and competition run and are trained as
#' usual. The regeneration term is weighted 0 and every regeneration parameter
#' (species parameters, environment response, cap, `theta`, custom parameters)
#' is frozen for that call. `predict()` and `simulateForest()` take
#' `recruits = "observed"` and the table `obs_recruits` to simulate the same way. With
#' observed recruitment the model is not an autonomous simulator: every
#' simulated year must lie in an observed interval, and bare-ground, succession
#' or projection runs are refused.
#'
#' The intended use is a staged calibration when the inventory samples small
#' trees too sparsely to calibrate recruitment jointly with the rest. First fit
#' the demography (growth, mortality, competition) with
#' `fit(loss_recruit = "observed")`; then fit the regeneration process on its own, with
#' the other processes frozen (`optimizeSpecies`/`optimizeEnv = FALSE`, or
#' through `loss_weights`).
#'
#' One caveat: the fed per-site count is a noisy estimate, and it is re-sampled
#' as if it were the mean. The observation noise of a sparse small-tree count
#' therefore enters the simulated stand. On the German national forest
#' inventory one tallied small tree stands for about 1 600 stems per ha per
#' decade, and the fed counts carry a variance comparable to that of the whole
#' stem-count column, which costs accuracy in the simulated stem counts.
#'
#' @section Choosing:
#' * **Default settings** (`regeneration_capped`, trained `theta` with
#'   `theta_floor = 0.5`, `"nbinom"`): the fit is stable, but do not read the
#'   fitted recruitment level as calibrated.
#' * **For a calibrated recruitment level**: keep `theta` fixed
#'   (`createProcess(optimizeTheta = FALSE)`, with `dispersion_parameter` 1 or
#'   higher), or use `fit(loss_family = c(..., regeneration = "poisson"))`.
#'   Mind that `recruit_obs_weight` is ignored under `"poisson"` (above).
#' * **When small trees are sampled too sparsely**, and in particular when the
#'   spatial pattern of recruitment is not predictable from the predictors (on
#'   the German national forest inventory no configuration ranked sites better
#'   than a Spearman correlation of about 0.08 out of sample): use the staged
#'   route of the previous section.
#' * **Check before interpreting recruitment**: compare the predicted and the
#'   observed total recruitment, and the regeneration loss in `m$history`
#'   against its intercept-only baseline (with `weights = "auto"`, a value of 1
#'   is no better than predicting the mean).
#'
#' @seealso [FINN::createProcess()], [FINN::createHybrid()], [FINN::finn()],
#'   [FINN::fit()], [FINN::simulateForest()], [FINN::regeneration()],
#'   [FINN::regeneration_capped()], [FINN::regeneration_adult()]
#'
#' @name recruitment
#' @aliases FINN-recruitment
NULL
