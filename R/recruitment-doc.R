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
#' * [FINN::regeneration_adult_effect()]: the adult density as an effect on
#'   recruitment instead of a gate,
#'   `f(A) = f0 + (1 - f0) A / (A + A_half)` with `f0 = plogis(reg_adult_logitf0)`,
#'   the share of the potential recruitment a species gets without adults on
#'   the site (seed rain from outside, advance regeneration). Prefer it to
#'   `regeneration_adult()` whenever a species must be able to establish where
#'   it has no adults: colonisation from an empty stand, and any free run from
#'   bare ground, need `f0 > 0`. `f0 -> 1` removes the adult effect, `f0 -> 0`
#'   is `regeneration_adult()`.
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
#' What the floor does: it keeps the fit stable, and in a fit run to
#' convergence it costs nothing we could resolve. In one converged pair of fits
#' on the German national forest inventory (2000 epochs, about 126 000 optimizer
#' steps), the model with the floor predicted the recruitment level on a
#' spatially blocked holdout at 1.06 times the observed and the one with a fixed
#' `theta` of 1 (`optimizeTheta = FALSE`) at 1.04; their recruitment likelihoods
#' were within 0.003 of each other (0.897 and 0.894 relative to an
#' intercept-only null), and they agreed on basal area, stem numbers, growth and
#' the mortality AUC to three decimals. Early in a fit the two settings differ:
#' after 50 of those 2000 epochs (3 150 optimizer steps on that inventory) a
#' floored `theta` still put the level about 12 times too high, a fixed `theta`
#' 2 to 3 times, and a free `theta` had run away by four orders of magnitude.
#' Count the budget in optimizer steps, not epochs: on a small data set 50
#' epochs can be 50 steps. So read the recruitment level only once
#' the parameters have settled (`$param_history`), whichever setting you use; for
#' short fits, a fixed `theta` gets the level right sooner.
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
#'   `theta_floor = 0.5`, `"nbinom"`): the fit is stable, and in one converged
#'   pair of fits on the German national forest inventory the recruitment level
#'   was calibrated as well as with a fixed `theta`. Do not read the level from a
#'   fit stopped early; check that `$param_history` has settled.
#' * **For short fits**, or when the recruitment level should not depend on how
#'   long the fit runs: keep `theta` fixed (`createProcess(optimizeTheta =
#'   FALSE)`, with `dispersion_parameter` 1 or higher), or use
#'   `fit(loss_family = c(..., regeneration = "poisson"))`. Mind that
#'   `recruit_obs_weight` is ignored under `"poisson"` (above).
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
#'   [FINN::regeneration_capped()], [FINN::regeneration_adult()],
#'   [FINN::regeneration_adult_effect()]
#'
#' @name recruitment
#' @aliases FINN-recruitment
NULL
