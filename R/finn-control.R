#' Settings of the optimizer used by fit()
#'
#' @description
#' Collects the settings of the optimization behind [FINN::fit()] in one object,
#' passed as `fit(control = finnControl(...))`. The defaults are robust in most
#' cases: a fit needs no `control` argument at all. The learning rate `lr`, the
#' number of `epochs` and the `batchsize` are changed much more often, so they
#' stay arguments of [FINN::fit()] itself.
#'
#' @details
#' [FINN::fit()] estimates the parameters by gradient descent: in every epoch it
#' simulates the sites, compares the simulation with the observations (see
#' [finn_loss]), and moves every parameter a small step in the direction that
#' makes the comparison better. `finnControl()` sets how those steps are taken.
#'
#' * `optimizer` is the rule that turns a gradient into a step. The default,
#'   Adam ([torch::optim_ignite_adam]), adapts the step size per parameter and
#'   works well without tuning. Any `torch` optimizer generator can be passed;
#'   extra arguments for it, such as `weight_decay`, go into `...`.
#' * `lr_scheduler` lowers the learning rate during the fit. The default
#'   `"none"` keeps `lr` constant. `"step"` divides the learning rate by a factor
#'   every `step_size` epochs, `"exponential"` multiplies it by `gamma` every
#'   epoch, `"cosine"` lets it fall smoothly to `eta_min` over `T_max` epochs,
#'   and `"plateau"` lowers it by `factor` after `patience` epochs without an
#'   improvement of the total loss. A schedule helps when the losses keep
#'   oscillating late in a fit. Advanced: a `function(optimizer)` that returns a
#'   `torch` learning-rate scheduler (an object with a `$step()` method).
#' * `lr_scheduler_params` tunes the chosen schedule, e.g.
#'   `list(step_size = 50, gamma = 0.5)` for `"step"`, `list(T_max = 200)` for
#'   `"cosine"`, `list(factor = 0.5, patience = 20)` for `"plateau"`. Entries
#'   that are not given take defaults scaled to `epochs` and `lr`.
#' * `clip_norm` limits how large a single step can be. Gradients are clipped to
#'   this norm, separately for three groups of parameters (the species
#'   parameters of the mechanistic processes, the environmental-response
#'   networks, and the parameters of the likelihoods), so that a group with
#'   large gradients cannot swamp the others. A single number applies the same
#'   budget to every group; a named list, e.g. `list(mechanistic = 2, nn = 1,
#'   loss = 5)`, sets groups individually.
#'
#' Before 0.2.0 these were arguments of [FINN::fit()] itself. The old form is
#' still accepted with a warning, e.g. `fit(lr_scheduler = "step")` becomes
#' `fit(control = finnControl(lr_scheduler = "step"))`.
#'
#' @param optimizer (`function`)\cr A `torch` optimizer generator. Default
#'   [torch::optim_ignite_adam].
#' @param lr_scheduler (`character(1)` or `function`)\cr Learning-rate schedule:
#'   `"none"` (default), `"step"`, `"exponential"`, `"cosine"`, `"plateau"`, or a
#'   `function(optimizer)` returning a `torch` scheduler. See Details.
#' @param lr_scheduler_params (`list`)\cr Settings of the chosen schedule. See
#'   Details.
#' @param clip_norm (`numeric(1)` or named `list`)\cr Gradient-norm budget per
#'   parameter group. Default `2`.
#' @param ... Further arguments for `optimizer`, e.g. `weight_decay = 1e-4`.
#'
#' @return A list of class `finnControl`, to be passed to `fit(control = )`.
#'
#' @seealso [FINN::fit()], [finn_loss]
#'
#' @examplesIf torch::torch_is_installed()
#' # the defaults, which most fits can keep
#' finnControl()
#'
#' # halve the learning rate every 50 epochs and allow larger steps
#' ctrl <- finnControl(lr_scheduler = "step",
#'                     lr_scheduler_params = list(step_size = 50, gamma = 0.5),
#'                     clip_norm = 5)
#'
#' library(data.table)
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' m <- finn(N_species = Nsp, recruits_dbh = 12.9,
#'           growth_process       = createProcess(~ temp + prec, growth),
#'           mortality_process    = createProcess(~ temp + prec, mortality),
#'           regeneration_process = createProcess(~ temp + prec, regeneration_capped))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
#'     control = ctrl, plot_progress = FALSE)
#' @export
finnControl = function(optimizer = torch::optim_ignite_adam,
                       lr_scheduler = "none",
                       lr_scheduler_params = list(),
                       clip_norm = 2.0,
                       ...) {
  if (!is.function(optimizer))
    stop("`optimizer` must be a torch optimizer generator, e.g. torch::optim_ignite_adam.", call. = FALSE)
  if (!(is.function(lr_scheduler) ||
        (is.character(lr_scheduler) && length(lr_scheduler) == 1L &&
         lr_scheduler %in% c("none", "step", "exponential", "cosine", "plateau"))))
    stop('`lr_scheduler` must be one of "none", "step", "exponential", "cosine", "plateau", or a function(optimizer).', call. = FALSE)
  if (!is.list(lr_scheduler_params))
    stop("`lr_scheduler_params` must be a list.", call. = FALSE)
  if (!(is.numeric(clip_norm) || is.list(clip_norm)))
    stop("`clip_norm` must be a number or a named list.", call. = FALSE)
  optimizer_args = list(...)
  if (length(optimizer_args) && (is.null(names(optimizer_args)) || any(names(optimizer_args) == "")))
    stop("arguments for the optimizer passed through `...` must be named, e.g. weight_decay = 1e-4.", call. = FALSE)
  structure(list(optimizer = optimizer,
                 lr_scheduler = lr_scheduler,
                 lr_scheduler_params = lr_scheduler_params,
                 clip_norm = clip_norm,
                 optimizer_args = optimizer_args),
            class = "finnControl")
}

#' @export
print.finnControl = function(x, ...) {
  sched = if (is.function(x$lr_scheduler)) "<function>" else x$lr_scheduler
  cat("FINN optimizer settings (fit(control = ))\n")
  cat("  optimizer:           ", if (identical(x$optimizer, torch::optim_ignite_adam)) "Adam (torch::optim_ignite_adam)" else "<custom torch optimizer>", "\n")
  cat("  lr_scheduler:        ", sched, "\n")
  if (length(x$lr_scheduler_params))
    cat("  lr_scheduler_params: ", paste(names(x$lr_scheduler_params), unlist(x$lr_scheduler_params), sep = " = ", collapse = ", "), "\n")
  cat("  clip_norm:           ", if (is.list(x$clip_norm)) paste(names(x$clip_norm), unlist(x$clip_norm), sep = " = ", collapse = ", ") else x$clip_norm, "\n")
  if (length(x$optimizer_args))
    cat("  optimizer arguments: ", paste(names(x$optimizer_args), unlist(x$optimizer_args), sep = " = ", collapse = ", "), "\n")
  invisible(x)
}
