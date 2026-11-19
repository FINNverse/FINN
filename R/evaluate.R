# Comparing simulations with observations: matching and aggregation (as fit()
# scores them), skill measures, comparePlot(), evaluate(), and the bare-ground
# environment used by simulateForest(bare_ground = TRUE).

response_levels = c("dbh", "ba", "trees", "growth", "mort", "reg")

# The simulation as a long site table (siteID, year, species, variable, value),
# from a finn_simulation, its $long$site table, or a wide site table.
simulation_long = function(sim) {
  if (inherits(sim, "finn_simulation") || (is.list(sim) && !is.data.frame(sim) && !is.null(sim$long)))
    sim = sim$long$site
  if (is.null(sim) || !is.data.frame(sim))
    stop("`sim` must be the result of simulateForest() or predict(), or its site table.", call. = FALSE)
  d = data.table::as.data.table(sim)
  if (!all(c("siteID", "year", "species") %in% names(d)))
    stop("The simulation table needs the columns siteID, year and species.", call. = FALSE)
  if (!all(c("variable", "value") %in% names(d))) {
    vars = intersect(c(response_levels, "r_mean_ha"), names(d))
    if (!length(vars)) stop("The simulation table holds none of the simulated variables.", call. = FALSE)
    d = data.table::melt(d, id.vars = c("siteID", "year", "species"), measure.vars = vars,
                         variable.name = "variable", value.name = "value")
  }
  d[, variable := as.character(variable)]
  d[, c("siteID", "year", "species", "variable", "value"), with = FALSE]
}

# The observations as a wide table (one column per response), from a wide table
# as makeObsData() returns it or a long one with `variable` / `value`.
observations_wide = function(obs) {
  if (!is.data.frame(obs))
    stop("`obs` must be a table of observations (as makeObsData() returns it).", call. = FALSE)
  d = data.table::as.data.table(obs)
  if (!all(c("siteID", "year", "species") %in% names(d)))
    stop("The observation table needs the columns siteID, year and species.", call. = FALSE)
  if (all(c("variable", "value") %in% names(d))) {
    keys = intersect(c("siteID", "year", "species", "period_length", "n_at_risk", "n_died", "species_name"), names(d))
    extra = unique(d[, keys, with = FALSE])
    w = data.table::dcast(d, siteID + year + species ~ variable, value.var = "value", fun.aggregate = mean)
    if (length(setdiff(keys, c("siteID", "year", "species"))))
      w = merge(w, unique(extra), by = c("siteID", "year", "species"), all.x = TRUE)
    d = w
  }
  if (!"period_length" %in% names(d)) d[, period_length := 1L]
  d
}

# Pair every observation with the simulation aggregated as fit() scores it: the
# states (dbh, ba, trees) at the census year, growth and mortality compounded
# over the period_length years before it (geometric mean of 1 + g minus one, one
# minus the geometric-mean survival), and the expected recruits summed over them.
# Returns one row per site, census, species and response with `obs` and `pred`,
# and for mortality the counts n_at_risk / n_died and the interval probability.
match_simulation = function(sim, obs, quiet = FALSE) {
  s = simulation_long(sim)
  o = observations_wide(obs)

  note = function(...) if (!quiet) message(sprintf(...))
  miss_sites = setdiff(unique(o$siteID), unique(s$siteID))
  if (length(miss_sites)) {
    note("%d observed site(s) are not in the simulation and are left out.", length(miss_sites))
    o = o[!siteID %in% miss_sites]
  }
  if (!nrow(o)) stop("No observed site is in the simulation: simulate the sites of `obs` (same siteID).", call. = FALSE)
  max_year = s[, list(max_year = max(year)), by = siteID]
  o = merge(o, max_year, by = "siteID")
  late = o$year > o$max_year
  if (any(late)) {
    note("%d observation(s) lie after the last simulated year of their site and are left out.", sum(late))
    o = o[!late]
  }
  o[, max_year := NULL]
  miss_species = setdiff(unique(o$species), unique(s$species))
  if (length(miss_species)) {
    note("Species %s are not in the simulation and are left out.", paste(sort(miss_species), collapse = ", "))
    o = o[!species %in% miss_species]
  }
  if (!nrow(o)) stop("No observation matches the simulation in site, year and species.", call. = FALSE)

  win = unique(o[, list(siteID, census = year, period_length)])
  win = win[, list(year = seq.int(census - period_length + 1L, census)), by = list(siteID, census)]
  sw = merge(s, win, by = c("siteID", "year"), allow.cartesian = TRUE)
  pred = data.table::rbindlist(list(
    sw[variable %in% c("dbh", "ba", "trees") & year == census,
       list(siteID, year = census, species, variable, pred = value)],
    sw[variable == "growth",
       list(pred = exp(mean(log1p(value))) - 1), by = list(siteID, year = census, species, variable)],
    sw[variable == "mort",
       list(pred = 1 - exp(mean(log1p(-pmin(value, 1 - 1e-7))))), by = list(siteID, year = census, species, variable)],
    sw[variable == "r_mean_ha",
       list(variable = "reg", pred = sum(value)), by = list(siteID, year = census, species)]
  ), use.names = TRUE)

  vars = intersect(response_levels, names(o))
  keep = c("siteID", "year", "species", "period_length",
           intersect(c("n_at_risk", "n_died", "species_name"), names(o)))
  ol = data.table::melt(o[, c(keep, vars), with = FALSE], id.vars = keep, measure.vars = vars,
                        variable.name = "variable", value.name = "obs")
  ol[, variable := as.character(variable)]
  op = merge(ol, pred, by = c("siteID", "year", "species", "variable"))
  op = op[is.finite(obs) & is.finite(pred)]
  # the interval death probability the annual rate implies, for the AUC
  op[variable == "mort", p_interval := 1 - (1 - pred)^period_length]
  op[, variable := factor(variable, levels = response_levels)]
  op[]
}

# AUC of the interval death probability against the observed fates: every tree
# at risk is one case (n_died positives, n_at_risk - n_died negatives sharing the
# row's probability); ties count one half.
auc_counts = function(p, pos, neg) {
  g = data.table::data.table(p = p, pos = pos, neg = neg)[, list(pos = sum(pos), neg = sum(neg)), by = p][order(p)]
  if (sum(g$pos) == 0 || sum(g$neg) == 0) return(NA_real_)
  neg_below = cumsum(g$neg) - g$neg
  sum(g$pos * (neg_below + 0.5 * g$neg)) / (sum(g$pos) * sum(g$neg))
}

# Skill per response (and per species with by_species): R2 (1 - SSE/SST) first,
# then RMSE, Pearson r, bias, the ratio of the mean prediction to the mean
# observation and n; for mortality also the tree-level AUC and the ratio of
# predicted to observed deaths when the counts are there.
skill_table = function(op, by_species = FALSE) {
  by = c("variable", if (by_species) "species")
  m = op[, list(n = .N,
                r2 = if (.N > 1L && stats::var(obs) > 0) 1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2) else NA_real_,
                rmse = sqrt(mean((pred - obs)^2)),
                r = if (.N > 2L && stats::sd(obs) > 0 && stats::sd(pred) > 0) stats::cor(obs, pred) else NA_real_,
                bias = mean(pred - obs),
                ratio = if (mean(obs) != 0) mean(pred) / mean(obs) else NA_real_),
         by = by]
  if (all(c("n_at_risk", "n_died") %in% names(op))) {
    mo = op[variable == "mort" & is.finite(n_at_risk) & n_at_risk > 0,
            list(auc = auc_counts(p_interval, n_died, n_at_risk - n_died),
                 deaths_ratio = if (sum(n_died) > 0) sum(p_interval * n_at_risk) / sum(n_died) else NA_real_),
            by = by]
    m = merge(m, mo, by = by, all.x = TRUE)
  } else {
    m[, `:=`(auc = NA_real_, deaths_ratio = NA_real_)]
  }
  m[, variable := factor(variable, levels = response_levels)]
  data.table::setorderv(m, by)
  m[]
}

species_label = function(op) {
  if ("species_name" %in% names(op)) {
    lab = unique(op[, list(species, species_name)])
    stats::setNames(paste0(lab$species, " ", lab$species_name), lab$species)
  } else NULL
}

#' Compare a simulation with observations
#'
#' Pairs a simulation with the observed data and draws, per response, how well
#' they agree. The simulation is aggregated the way [fit()] scores it: the stand
#' states `dbh`, `ba` and `trees` at the census year, `growth` and `mort`
#' compounded over the `period_length` years before the census (the geometric
#' mean of `1 + g` minus one; one minus the geometric-mean survival), and the
#' expected recruits (`r_mean_ha`) summed over them for `reg`. Observations whose
#' site, year or species are not in the simulation are left out with a message.
#' Both tables may come in long (`variable`, `value`) or wide form.
#'
#' * `type = "scatter"` (default): observed against simulated for every
#'   response, with the 1:1 line and R2 (`1 - SSE/SST`) in each panel; for
#'   mortality the panel shows the interval death probability against the
#'   observed share of deaths and prints the tree-level AUC instead (each tree
#'   at risk is one case, `n_died` positives and `n_at_risk - n_died`
#'   negatives; this needs those two columns in `obs`).
#' * `type = "timeseries"`: the simulated trajectory over the years as a line
#'   and the observations as points at their census years, averaged over the
#'   selected sites.
#'
#' `plot(sim, obs)` is the same as `comparePlot(sim, obs)`.
#'
#' @param sim The result of [simulateForest()] or [predict()] for the sites of
#'   `obs`, started from the same initial state (or its site table).
#' @param obs Observations as [makeObsData()] returns them (wide), or long with
#'   `variable` and `value`; `period_length`, `n_at_risk` and `n_died` are used
#'   where present.
#' @param type `"scatter"` (default) or `"timeseries"`.
#' @param species `NULL` (default) pools all species in one panel per response;
#'   species numbers draw one row of panels per species.
#' @param response Responses to draw, any of `"dbh"`, `"ba"`, `"trees"`,
#'   `"growth"`, `"mort"`, `"reg"`; `NULL` (default) all that are observed.
#' @param sites Site ids to include; `NULL` (default) all.
#' @param ... Unused.
#' @return The `ggplot` object, invisibly; it is printed. The paired table is
#'   attached as attribute `"matched"` and the skill measures as `"skill"`.
#' @seealso [evaluate()], which stores the same comparison on the model, and
#'   [summary.finn_class()], which prints it.
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
#' sim <- simulateForest(m, env = fia_env_dt[siteID %in% sites], init_cohort = init)
#' comparePlot(sim, fia_obs_dt[siteID %in% sites])
#' plot(sim, fia_obs_dt[siteID %in% sites], type = "timeseries", species = 1:2)
#' @export
comparePlot = function(sim, obs, type = c("scatter", "timeseries"), species = NULL,
                       response = NULL, sites = NULL, ...) {
  type = match.arg(type)
  op = match_simulation(sim, obs)
  if (!is.null(sites)) {
    op = op[siteID %in% sites]
    if (!nrow(op)) stop("None of `sites` has a matching observation.", call. = FALSE)
  }
  if (!is.null(response)) {
    bad = setdiff(response, response_levels)
    if (length(bad)) stop(sprintf("Unknown `response` %s; use any of %s.", paste(bad, collapse = ", "),
                                  paste(response_levels, collapse = ", ")), call. = FALSE)
    op = op[as.character(variable) %in% response]
  }
  if (!is.null(species)) {
    sp_sel = species_subset_obs(species, op)   # outside [ ]: `species` is also a column
    op = op[species %in% sp_sel]
  }
  if (!nrow(op)) stop("Nothing left to compare after the selection.", call. = FALSE)
  op[, variable := droplevels(variable)]
  labs = species_label(op)
  op[, species_lab := if (is.null(labs)) paste("species", species) else labs[as.character(species)]]
  skill = skill_table(op, by_species = !is.null(species))

  facet = if (is.null(species)) ggplot2::facet_wrap(~ variable, scales = "free")
          else ggplot2::facet_grid(species_lab ~ variable, scales = "free")

  if (type == "scatter") {
    d = data.table::copy(op)
    # mortality: interval probability against the observed share of deaths
    if (all(c("n_at_risk", "n_died") %in% names(d))) {
      d[variable == "mort", `:=`(pred = p_interval,
                                 obs = ifelse(n_at_risk > 0, n_died / n_at_risk, NA_real_))]
    }
    d = d[is.finite(obs) & is.finite(pred)]
    sk = data.table::copy(skill)
    if (!is.null(species)) sk[, species_lab := if (is.null(labs)) paste("species", species) else labs[as.character(species)]]
    sk[, label := ifelse(variable == "mort" & is.finite(auc), sprintf("AUC = %.2f", auc),
                         ifelse(is.finite(r2), sprintf("R² = %.2f", r2), "R² = NA"))]
    p = ggplot2::ggplot(d, ggplot2::aes(x = pred, y = obs)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
      ggplot2::geom_point(alpha = 0.35, size = 1) +
      ggplot2::geom_text(data = sk, ggplot2::aes(label = label), x = -Inf, y = Inf,
                         hjust = -0.1, vjust = 1.4, size = 3.4, inherit.aes = FALSE) +
      facet +
      ggplot2::labs(x = "simulated", y = "observed",
                    caption = "mort: interval death probability vs observed share of deaths, with the tree-level AUC") +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  } else {
    s = simulation_long(sim)
    s[variable == "r_mean_ha", variable := "reg_annual"]
    s = s[siteID %in% unique(op$siteID) & variable %in% as.character(unique(op$variable))]
    if (!is.null(species)) s = s[species %in% unique(op$species)]
    by_sp = !is.null(species)
    sl = s[, list(value = mean(value, na.rm = TRUE)), by = c("year", "variable", if (by_sp) "species")]
    ob = op[, list(obs = mean(obs, na.rm = TRUE)), by = c("year", "variable", if (by_sp) "species")]
    if (!by_sp) {
      # all species pooled: stand totals for ba and trees, the mean over species otherwise
      sl = s[, list(value = if (variable[1] %in% c("ba", "trees")) sum(value, na.rm = TRUE) else mean(value, na.rm = TRUE)),
             by = list(siteID, year, variable)][, list(value = mean(value)), by = list(year, variable)]
      ob = op[, list(obs = if (as.character(variable[1]) %in% c("ba", "trees")) sum(obs, na.rm = TRUE) else mean(obs, na.rm = TRUE)),
              by = list(siteID, year, variable)][, list(obs = mean(obs)), by = list(year, variable)]
    } else {
      sl[, species_lab := if (is.null(labs)) paste("species", species) else labs[as.character(species)]]
      ob[, species_lab := if (is.null(labs)) paste("species", species) else labs[as.character(species)]]
    }
    sl[, variable := factor(variable, levels = levels(op$variable))]
    ob[, variable := factor(as.character(variable), levels = levels(op$variable))]
    p = ggplot2::ggplot() +
      ggplot2::geom_line(data = sl, ggplot2::aes(x = year, y = value)) +
      ggplot2::geom_point(data = ob, ggplot2::aes(x = year, y = obs), colour = "#C0392B", size = 2) +
      facet +
      ggplot2::labs(x = "simulated year", y = if (by_sp) "mean over sites" else "mean over sites (ba, trees: stand total)",
                    caption = "line: simulated (annual rates for growth and mort); points: observed at the census") +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  }
  attr(p, "matched") = op
  attr(p, "skill") = skill
  print(p)
  invisible(p)
}

species_subset_obs = function(species, op) {
  species = as.integer(species)
  bad = setdiff(species, unique(op$species))
  if (length(bad)) stop(sprintf("Species %s have no matching observation.", paste(bad, collapse = ", ")), call. = FALSE)
  species
}

# The environment of a bare-ground run: the mean of each predictor over the
# training environment (one site, siteID 1) or over each of `sites` (that
# site's mean, siteID kept), held constant for `years` years.
bare_ground_env = function(model, sites = NULL, years = 500L, env = NULL) {
  base = if (!is.null(env)) data.table::as.data.table(env) else model$train_env
  if (is.null(base))
    stop("No environment to take the mean climate from: fit the model first or pass `env`.", call. = FALSE)
  years = as.integer(years)
  if (length(years) != 1L || is.na(years) || years < 1L) stop("`years` must be one positive number.", call. = FALSE)
  base = data.table::as.data.table(base)
  preds = setdiff(names(base), c("siteID", "year"))
  preds = preds[vapply(preds, function(cc) is.numeric(base[[cc]]), logical(1))]
  if (is.null(sites)) {
    means = base[, lapply(.SD, mean, na.rm = TRUE), .SDcols = preds][, siteID := 1L]
  } else {
    bad = setdiff(sites, unique(base$siteID))
    if (length(bad)) stop(sprintf("`sites` %s are not in the environment the model was fitted on.", paste(bad, collapse = ", ")), call. = FALSE)
    means = base[siteID %in% sites, lapply(.SD, mean, na.rm = TRUE), by = siteID, .SDcols = preds]
  }
  out = means[rep(seq_len(nrow(means)), each = years)]
  out[, year := rep(seq_len(years), times = nrow(means))]
  data.table::setcolorder(out, c("siteID", "year", preds))
  out[]
}

#' Evaluate a fitted model
#'
#' An optional step after [fit()] that simulates the model and stores what is
#' needed to judge it on the model object, as `model$evaluation`; [summary()]
#' prints it and [comparePlot()] draws the same comparison. Calling
#' `fit(..., holdout = )` runs it right after the fit.
#'
#' * **Training data**: the model is simulated from the training initial state
#'   under the training environment and compared with the observations it was
#'   fitted to, aggregated as [fit()] scores them.
#' * **Holdout** (`holdout = list(data =, env =, init_cohort =)`): the same on
#'   data the model has not seen, e.g. a spatially separate set of sites.
#' * **Bare ground** (`bare_ground = TRUE`): a run from an empty stand for
#'   `years` years under the mean training climate held constant, or under the
#'   mean climate of each of `sites` (see [simulateForest()]). Optional, because
#'   a long run takes time.
#'
#' @param model A model fitted with [fit()].
#' @param holdout `NULL` (default) or a list with the holdout `data`
#'   (observations as for [fit()]), its `env` and its `init_cohort`, and
#'   optionally `disturbance`.
#' @param bare_ground Run the bare-ground simulation too. Default `FALSE`.
#' @param years Length of the bare-ground run in years. Default `500`.
#' @param sites `NULL` (default) one bare-ground run under the mean training
#'   climate; site ids for one run per site under that site's mean climate.
#' @param patches Patches per site; `NULL` (default) as in the fit.
#' @param train `NULL` (default): the training data the model keeps from
#'   [fit()]. A list like `holdout` for a model that keeps none (fitted before
#'   this version of FINN, e.g. read with `torch::torch_load()`).
#' @return The model, invisibly, with `model$evaluation`: a list with `train`
#'   and, when given, `holdout` (each with the paired table `matched` and the
#'   skill measures `skill` overall and `skill_species` per species), and
#'   `bare_ground` (the simulation).
#' @seealso [comparePlot()], [summary.finn_class()], [simulateForest()]
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' Nsp  <- max(fia_obs_dt$species)
#' train <- 1:3; hold <- 4:5
#' init_tr <- makeInitCohorts(fia_init_trees[siteID %in% train], Nspecies = Nsp)
#' init_ho <- makeInitCohorts(fia_init_trees[siteID %in% hold], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(N_species = Nsp, recruits_dbh = 12.9,
#'           growth_process       = createProcess(~ temp + prec, growth),
#'           mortality_process    = createProcess(~ temp + prec, mortality),
#'           regeneration_process = createProcess(~ temp + prec, regeneration_capped))
#' fit(m, data = fia_obs_dt[siteID %in% train], env = fia_env_dt[siteID %in% train],
#'     init_cohort = init_tr, epochs = 2L, patches = 4, patch_size = 0.06,
#'     plot_progress = FALSE)
#' evaluate(m, holdout = list(data = fia_obs_dt[siteID %in% hold],
#'                            env = fia_env_dt[siteID %in% hold], init_cohort = init_ho),
#'          bare_ground = TRUE, years = 50)
#' m$evaluation$holdout$skill
#' @export
evaluate = function(model, holdout = NULL, bare_ground = FALSE, years = 500L, sites = NULL,
                    patches = NULL, train = NULL) {
  if (!inherits(model, "finn_class")) stop("`model` must be a FINN model built by finn().", call. = FALSE)
  check_set = function(x, what) {
    if (!is.null(x) && (!is.list(x) || is.null(x$data) || is.null(x$env)))
      stop(sprintf("`%s` must be a list with `data`, `env` and `init_cohort` (and optionally `disturbance`).", what), call. = FALSE)
  }
  check_set(holdout, "holdout"); check_set(train, "train")
  if (is.null(train) && !is.null(model$train_data))
    train = list(data = model$train_data, env = model$train_env,
                 init_cohort = model$train_init_cohort$cohort, disturbance = model$train_disturbance)
  if (is.null(train) && is.null(holdout) && !isTRUE(bare_ground))
    stop("Nothing to evaluate: this model keeps no training data (it was fitted before this version of FINN or not at all); pass `train = list(data =, env =, init_cohort =)` or a `holdout`.", call. = FALSE)
  if (is.null(patches)) patches = if (is.null(model$patches)) 100L else model$patches

  one_set = function(data, env, init_cohort, disturbance = NULL) {
    sim = predict(model, env = env, disturbance = disturbance, init_cohort = init_cohort,
                  patches = patches)
    op = match_simulation(sim, data, quiet = TRUE)
    list(matched = op, skill = skill_table(op), skill_species = skill_table(op, by_species = TRUE))
  }

  ev = list()
  if (!is.null(train))
    ev$train = one_set(train$data, train$env, train$init_cohort, train$disturbance)
  if (!is.null(holdout)) {
    if (is.null(holdout$init_cohort))
      warning("`holdout` has no `init_cohort`: the holdout sites start from bare ground.", call. = FALSE)
    ev$holdout = one_set(holdout$data, holdout$env, holdout$init_cohort, holdout$disturbance)
  }
  if (isTRUE(bare_ground)) {
    ev$bare_ground = simulateForest(model, bare_ground = TRUE, years = years, sites = sites, patches = patches)
    ev$bare_ground_settings = list(years = years, sites = sites, patches = patches)
  }
  model$evaluation = ev
  invisible(model)
}
