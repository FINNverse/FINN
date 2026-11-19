# Plot methods for a fitted model (finn_class) and for a simulation
# (finn_simulation, the list returned by simulateForest() / predict()).

# Names of the species-specific parameters of FINN's default process functions,
# in the column order of model$par_<process>; documented in ?createProcess.
default_parameter_names = list(
  growth       = list(func = "growth",       names = c("light requirement", "size decay")),
  mortality    = list(func = "mortality",    names = c("shade", "size", "growth")),
  regeneration = list(func = "regeneration", names = c("light requirement")),
  competition  = list(func = "competition",  names = c("height", "strength"))
)

# The labels of the columns of model$par_<type>: the documented names when the
# process uses FINN's own function (or regeneration_capped(), which shares the
# parameter of regeneration()), else "parameter 1", "parameter 2", ...
parameter_labels = function(model, type, k) {
  proc = model[[paste0("process_", type)]]
  def  = default_parameter_names[[type]]
  same_body = function(f, g) !is.null(f) && is.function(f) && identical(body(f), body(g))
  funcs = list(getExportedValue("FINN", def$func))
  if (type == "regeneration") funcs = c(funcs, list(FINN::regeneration_capped))
  known = !is.null(proc) && any(vapply(funcs, function(g) same_body(proc$func, g), logical(1)))
  if (known && length(def$names) == k) return(def$names)
  paste("parameter", seq_len(k))
}

# TRUE when the environmental response of a process is linear in its predictors:
# not a hybrid, no hidden layers and no custom network.
env_is_linear = function(proc) {
  !is.null(proc) && !inherits(proc, "hybrid") && length(proc$hidden) == 0L && is.null(proc$NN)
}

species_subset = function(species, N) {
  if (is.null(species)) return(seq_len(N))
  species = as.integer(species)
  if (anyNA(species) || any(species < 1L | species > N))
    stop(sprintf("`species` must be species numbers between 1 and %d.", N), call. = FALSE)
  sort(unique(species))
}

# One row per process parameter and species: the fitted value on the scale the
# parameter acts on, and the admissible range from createProcess(lower, upper).
process_parameter_table = function(model, species) {
  rows = list()
  for (type in c("growth", "mortality", "regeneration", "competition")) {
    mat = tryCatch(model[[paste0("par_", type, "_r")]], error = function(e) NULL)
    if (is.null(mat) || !length(mat)) next
    mat = as.matrix(mat)
    lower = as.numeric(torch::as_array(model[[paste0("par_", type, "_lower")]]$cpu()))
    upper = as.numeric(torch::as_array(model[[paste0("par_", type, "_upper")]]$cpu()))
    labs  = parameter_labels(model, type, ncol(mat))
    for (j in seq_len(ncol(mat))) {
      rows[[length(rows) + 1L]] = data.frame(
        process = type, parameter = labs[j], species = species,
        value = mat[species, j], lower = lower[j], upper = upper[j],
        stringsAsFactors = FALSE)
    }
  }
  if (!is.null(model$par_theta_recruits_raw)) {
    theta = as.numeric(torch::as_array(model$par_theta_recruits$cpu()))
    if (length(theta) == 1L) theta = rep(theta, model$N_species)
    floor_theta = if (is.null(model$theta_floor)) 0 else model$theta_floor
    rows[[length(rows) + 1L]] = data.frame(
      process = "regeneration", parameter = "recruit dispersion theta", species = species,
      value = theta[species], lower = floor_theta, upper = Inf,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  df = do.call(rbind, rows)
  span = df$upper - df$lower
  tol  = ifelse(is.finite(span), 1e-3 * span, 1e-3 * pmax(abs(df$lower), 1))
  df$at_bound = (df$value - df$lower) <= tol | (is.finite(df$upper) & (df$upper - df$value) <= tol)
  df
}

# The fitted environmental term of each linear process, per species, over the
# range each predictor covers in `env`; the other predictors at their mean.
niche_table = function(model, env, species, n_grid = 100L) {
  env = data.table::as.data.table(env)
  numeric_cols = setdiff(names(env), c("siteID", "year"))
  numeric_cols = numeric_cols[vapply(numeric_cols, function(cc) is.numeric(env[[cc]]), logical(1))]
  centre = lapply(env[, numeric_cols, with = FALSE], mean, na.rm = TRUE)
  out = list()
  for (type in c("growth", "mortality", "regeneration")) {
    proc = model[[paste0("process_", type)]]
    nn   = model[[paste0("nn_", type)]]
    form = model[[paste0(type, "_formula")]]
    if (is.null(proc) || is.null(nn) || is.null(form)) next
    if (!env_is_linear(proc)) {
      warning(sprintf(paste0("The environmental response of the %s process is a neural network, ",
                             "so there are no coefficients to show; its panel is skipped. ",
                             "Use ALE() to look at it."), type), call. = FALSE)
      next
    }
    vars = all.vars(form)
    if (!length(vars)) next
    missing_vars = setdiff(vars, numeric_cols)
    if (length(missing_vars)) {
      warning(sprintf("The %s formula uses %s, which is not a numeric column of `env`; its panel is skipped.",
                      type, paste(missing_vars, collapse = ", ")), call. = FALSE)
      next
    }
    p1 = nn$parameters[[1]]
    for (v in vars) {
      rng  = range(env[[v]], na.rm = TRUE)
      grid = data.table::as.data.table(centre)[rep(1L, n_grid)]
      data.table::set(grid, j = v, value = seq(rng[1], rng[2], length.out = n_grid))
      grid[, `:=`(siteID = 1L, year = 1L)]
      if (!is.null(model$env_scaling)) grid = apply_env_scaling(grid, model$env_scaling)
      X = stats::model.matrix(form, as.data.frame(grid[, numeric_cols, with = FALSE]))
      Xt = torch::torch_tensor(X, dtype = p1$dtype, device = p1$device)
      pred = torch::with_no_grad(as.matrix(torch::as_array(nn(Xt)$cpu())))
      x_raw = seq(rng[1], rng[2], length.out = n_grid)
      out[[length(out) + 1L]] = data.frame(
        process = type, predictor = v,
        x = rep(x_raw, times = length(species)),
        species = rep(species, each = n_grid),
        value = as.numeric(pred[, species, drop = FALSE]),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

#' Plot the fitted parameters of a FINN model
#'
#' Draws the fitted parameters themselves, per species, on the scale on which
#' they enter the model -- no simulation is involved. By default the process
#' parameters; `pars = "env"` switches to the environmental niches.
#'
#' * **Process parameters** (`pars = "process"`, the default): one panel per parameter of the
#'   growth, mortality, regeneration and competition processes, and the recruit
#'   dispersion `theta`; one point per species, on the constrained scale the
#'   process function uses. The admissible range set by
#'   `createProcess(lower, upper)` is the light band behind the points (for
#'   `theta` it starts at the floor `createProcess(theta_floor)`), so a species
#'   whose parameter ended the fit at a bound -- drawn as an open triangle -- is
#'   visible at a glance (see `?createProcess` on what such a parameter means).
#' * **Environmental niches** (`pars = "env"`): for each process whose
#'   environmental response is linear (no `hidden`, no custom `NN`, not
#'   [createHybrid()]), the fitted environmental term of every species across
#'   the range each predictor covers in `env`, with the other predictors held at
#'   their mean. The curve is the linear predictor the process function
#'   receives (`pred`): for FINN's default functions growth and regeneration use
#'   `exp(pred)` and mortality `plogis(pred + ...)`, so a curve higher up means
#'   faster growth, more recruits or more deaths. The x axis is on the original
#'   scale of the predictor; the stored standardisation of [fit()] is applied
#'   internally. A process with a network response is skipped with a warning:
#'   it has no coefficients to show, and [ALE()] is the way to look at it.
#'   Competition does not depend on the environment.
#'
#' @param x A model built by [finn()], usually fitted with [fit()].
#' @param pars Which parameters: `"process"` (default) the species-specific
#'   process parameters, `"env"` the environmental niches. The argument and its
#'   two values are those of the `plot()` method of FINN 0.1.0, so a call such
#'   as `plot(m, pars = "env")` keeps working.
#' @param species Species numbers to draw; `NULL` (default) all.
#' @param env Table whose predictor ranges the niches span; `NULL` (default)
#'   the environment the model was fitted on. Only used with `pars = "env"`.
#' @param ... Unused.
#' @return The `ggplot` object, invisibly (`NULL` when there is nothing to
#'   draw); it is printed.
#' @seealso [ALE()] for the effect of a network response, [summary.finn_class()]
#'   for importances, [plot.finn_simulation()] for a simulation.
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(N_species = Nsp, recruits_dbh = 12.9,
#'           growth_process       = createProcess(~ temp + prec, growth),
#'           mortality_process    = createProcess(~ temp + prec, mortality),
#'           regeneration_process = createProcess(~ temp + prec, regeneration_capped))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
#'     plot_progress = FALSE)
#' plot(m)
#' plot(m, pars = "env", species = 1:3)
#' @export
plot.finn_class = function(x, pars = c("process", "env"), species = NULL, env = NULL, ...) {
  pars = match.arg(pars)
  sp = species_subset(species, x$N_species)
  out = NULL

  if (pars == "process") {
    df = process_parameter_table(x, sp)
    if (is.null(df)) {
      message("This model has no species-specific process parameters to plot.")
    } else {
      df$panel   = factor(paste0(df$process, ": ", df$parameter),
                          levels = unique(paste0(df$process, ": ", df$parameter)))
      df$species = factor(df$species, levels = sp)
      bands = unique(df[, c("panel", "lower", "upper")])
      p = ggplot2::ggplot(df, ggplot2::aes(x = value, y = species)) +
        ggplot2::geom_rect(data = bands, inherit.aes = FALSE,
                           ggplot2::aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
                           fill = "grey85", alpha = 0.6) +
        ggplot2::geom_point(ggplot2::aes(shape = at_bound), size = 2) +
        ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2),
                                    labels = c(`FALSE` = "inside", `TRUE` = "at a bound"),
                                    name = NULL, drop = FALSE) +
        ggplot2::facet_wrap(~ panel, scales = "free_x") +
        ggplot2::labs(x = "fitted value (band: admissible range)", y = "species") +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
      print(p)
      out = p
    }
  }

  if (pars == "env") {
    env_tab = if (is.null(env)) x$train_env else env
    if (is.null(env_tab)) {
      warning("No environment to take the predictor ranges from: fit the model first or pass `env`.", call. = FALSE)
    } else {
      df = niche_table(x, env_tab, sp)
      if (is.null(df)) {
        message("No linear environmental response to plot.")
      } else {
        df$process = factor(df$process, levels = intersect(c("growth", "mortality", "regeneration"), unique(df$process)))
        df$species = factor(df$species, levels = sp)
        p = ggplot2::ggplot(df, ggplot2::aes(x = x, y = value, colour = species, group = species)) +
          ggplot2::geom_line() +
          ggplot2::facet_grid(process ~ predictor, scales = "free") +
          ggplot2::labs(x = "predictor (original scale)",
                        y = "environmental term (linear predictor of the process)",
                        colour = "species") +
          ggplot2::theme_minimal() +
          ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
        print(p)
        out = p
      }
    }
  }
  invisible(out)
}

#' Plot a FINN simulation
#'
#' Draws the trajectories returned by [simulateForest()] or [predict()]: one
#' panel per response, one line per species, averaged over the sites, with the
#' spread between sites (10 % to 90 % quantile) as a light band when there is
#' more than one site.
#'
#' @param x The object returned by [simulateForest()] or [predict()] (class
#'   `finn_simulation`).
#' @param response Responses to draw, any of `"dbh"`, `"ba"`, `"trees"`,
#'   `"growth"`, `"mort"`, `"reg"`, `"r_mean_ha"`; `NULL` (default) the stand
#'   variables `ba`, `dbh` and `trees`.
#' @param sites Site ids to include; `NULL` (default) all.
#' @param ... Unused.
#' @return The `ggplot` object, invisibly; it is printed.
#' @seealso [plot.finn_class()] for the fitted parameters.
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(N_species = Nsp, recruits_dbh = 12.9,
#'           growth_process       = createProcess(~ temp + prec, growth),
#'           mortality_process    = createProcess(~ temp + prec, mortality),
#'           regeneration_process = createProcess(~ temp + prec, regeneration_capped))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
#'     plot_progress = FALSE)
#' sim <- simulateForest(m, env = fia_env_dt[siteID %in% sites], init_cohort = init,
#'                       patches = 4, patch_size = 0.06)
#' plot(sim)
#' plot(sim, response = c("growth", "mort"), sites = 1)
#' @export
plot.finn_simulation = function(x, response = NULL, sites = NULL, ...) {
  d = x$long$site
  if (is.null(d) || !nrow(d))
    stop("This simulation holds no site table (`$long$site`) to plot.", call. = FALSE)
  levels_all = c("dbh", "ba", "trees", "growth", "mort", "reg", "r_mean_ha")
  if (is.null(response)) response = c("ba", "dbh", "trees")
  bad = setdiff(response, levels_all)
  if (length(bad))
    stop(sprintf("Unknown `response` %s; use any of %s.", paste(bad, collapse = ", "),
                 paste(levels_all, collapse = ", ")), call. = FALSE)
  d = data.table::as.data.table(d)[as.character(variable) %in% response]
  if (!is.null(sites)) {
    d = d[siteID %in% sites]
    if (!nrow(d)) stop("None of `sites` is in the simulation.", call. = FALSE)
  }
  agg = d[, list(mean = mean(value, na.rm = TRUE),
                 lo = stats::quantile(value, 0.1, na.rm = TRUE, names = FALSE),
                 hi = stats::quantile(value, 0.9, na.rm = TRUE, names = FALSE),
                 n_sites = data.table::uniqueN(siteID)),
          by = list(variable = as.character(variable), year, species)]
  agg[, variable := factor(variable, levels = intersect(levels_all, response))]
  agg[, species := factor(species)]
  p = ggplot2::ggplot(agg, ggplot2::aes(x = year, y = mean, colour = species, group = species))
  if (any(agg$n_sites > 1L))
    p = p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = species),
                                 alpha = 0.08, colour = NA, show.legend = FALSE)
  p = p + ggplot2::geom_line() +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::labs(x = "simulated year", y = "mean over sites", colour = "species") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  print(p)
  invisible(p)
}

#' @rdname plot.finn_simulation
#' @export
print.finn_simulation = function(x, ...) {
  print(unclass(x), ...)
  invisible(x)
}
