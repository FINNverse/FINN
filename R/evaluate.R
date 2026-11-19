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
  o = data.table::copy(o)
  for (v in vars) data.table::set(o, j = v, value = as.numeric(o[[v]]))
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
# observation (for mortality weighted by the trees at risk) and n; for
# mortality also the tree-level AUC and the ratio of
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
                 deaths_ratio = if (sum(n_died) > 0) sum(p_interval * n_at_risk) / sum(n_died) else NA_real_,
                 # the ratio of the means weighted by the trees at risk: a plot with
                 # one tree must not count as much as one with fifty
                 ratio_w = if (sum(n_at_risk * obs) > 0) sum(n_at_risk * pred) / sum(n_at_risk * obs) else NA_real_),
            by = by]
    m = merge(m, mo, by = by, all.x = TRUE)
    m[variable == "mort" & is.finite(ratio_w), ratio := ratio_w]
    m[, ratio_w := NULL]
  } else {
    m[, `:=`(auc = NA_real_, deaths_ratio = NA_real_)]
  }
  m[, variable := factor(variable, levels = response_levels)]
  data.table::setorderv(m, by)
  m[]
}

# Species names, as a character vector named by the species number 1..N. Taken
# from, first found: `species_names` (a character vector in species order, a
# vector named by species number, or a table with `species` and
# `species_name`), the names a model keeps from fit(), a table `tab` with
# `species_name`; a species without a name is "species k".
species_names_lookup = function(N, species_names = NULL, model = NULL, tab = NULL) {
  out = stats::setNames(paste("species", seq_len(N)), seq_len(N))
  src = species_names
  if (is.null(src) && !is.null(model)) src = model$species_names
  if (is.null(src) && is.data.frame(tab) && all(c("species", "species_name") %in% names(tab)))
    src = unique(as.data.frame(tab)[, c("species", "species_name")])
  put = function(sp, nm) {
    sp = suppressWarnings(as.integer(sp)); ok = !is.na(sp) & sp >= 1L & sp <= N & !is.na(nm) & nzchar(nm)
    out[sp[ok]] <<- as.character(nm[ok])
  }
  if (is.data.frame(src)) put(src$species, src$species_name)
  else if (!is.null(names(src))) put(names(src), src)
  else if (length(src)) put(seq_along(src), src)
  out
}

#' Compare a simulation with observations
#'
#' Pairs a simulation with the observed data and draws, per response, how well
#' they agree. The simulation is aggregated the way [fit()] scores it: the stand
#' states `dbh`, `ba` and `trees` at the census year, `growth` and `mort`
#' compounded over the `period_length` years before the census (the geometric
#' mean of `1 + g` minus one; one minus the geometric-mean survival), and the
#' expected recruits (`r_mean_ha`) summed over them for `reg` (see [finn_loss]
#' and, for the observed columns, [finn_data]). Observations of a site or a
#' species that is not in the simulation, or of a year after the last simulated
#' one, are left out with a message.
#' Both tables may come in long (`variable`, `value`) or wide form.
#'
#' * `type = "scatter"` (default): observed against simulated for `dbh`, `ba`,
#'   `trees` and `growth`, with the 1:1 line and R2 (`1 - SSE/SST`) in each
#'   panel; for mortality the ROC curve of the interval death probability with
#'   its AUC (each tree at risk is one case, `n_died` positives and
#'   `n_at_risk - n_died` negatives; without those two columns in `obs`, a
#'   scatter of the annual rates); for recruitment, whose site-level values are
#'   too noisy for a correlation, the ratio of simulated to observed
#'   recruitment summed over all sites, per species (the `n_species` species
#'   with most observed recruitment).
#' * `type = "timeseries"`: the simulated trajectory over the years as a line
#'   and the observations as points at their census years, averaged over the
#'   selected sites. With all species pooled, `ba` and `trees` are the stand
#'   totals and the other responses the mean over the species.
#'
#' `plot(sim, obs)` is the same as `comparePlot(sim, obs)` (see
#' [plot.finn_simulation()]).
#'
#' @param sim The result of [simulateForest()] or
#'   [`predict()`][predict.finn_class] for the sites of
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
#' @param n_species Number of species in the recruitment panel. Default `12`.
#' @param species_names Species names for the labels: a character vector in
#'   species order, a vector named by species number, or a table with
#'   `species` and `species_name`. `NULL` (default): the names the simulation
#'   carries from the model, else the `species_name` column of `obs`.
#' @param x A `finn_comparison` object (the scatter view of `comparePlot()`).
#' @param ... Unused.
#' @return `type = "scatter"`: a `finn_comparison` object (a list of the
#'   `ggplot` panels `scatter`, `mort` and `reg`, which `print()` arranges on
#'   one page), invisibly; `type = "timeseries"`: the `ggplot` object,
#'   invisibly. It is printed. The paired table is
#'   attached as attribute `"matched"` and the skill measures as `"skill"`
#'   (per response, and per species when `species` is given; the measures are
#'   described in [summary.finn_class()]).
#' @seealso [evaluate()], which stores the same comparison on the model, and
#'   [summary.finn_class()], which prints it; [plot.finn_simulation()] for the
#'   simulation alone.
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
                       response = NULL, sites = NULL, n_species = 12L, species_names = NULL, ...) {
  type = match.arg(type)
  op = match_simulation(sim, obs)
  N = max(c(op$species, simulation_long(sim)$species))
  if (is.null(species_names)) species_names = attr(sim, "species_names")
  labs = species_names_lookup(N, species_names, tab = obs)
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
  op[, species_lab := labs[as.character(species)]]
  skill = skill_table(op, by_species = !is.null(species))

  facet = if (is.null(species)) ggplot2::facet_wrap(~ variable, scales = "free")
          else ggplot2::facet_grid(species_lab ~ variable, scales = "free")

  theme_finn = ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (type == "scatter") {
    panels = list()
    # (1) the states and growth: observed against simulated, R2 in each panel
    sc = op[variable %in% c("dbh", "ba", "trees", "growth")]
    if (nrow(sc)) {
      sc[, variable := droplevels(variable)]
      sk = skill[variable %in% c("dbh", "ba", "trees", "growth")]
      if (!is.null(species)) sk[, species_lab := labs[as.character(species)]]
      sk[, label := ifelse(is.finite(r2), sprintf("R\u00b2 = %.2f", r2), "R\u00b2 = NA")]
      fct = if (is.null(species)) ggplot2::facet_wrap(~ variable, scales = "free", nrow = 1)
            else ggplot2::facet_grid(species_lab ~ variable, scales = "free")
      panels$scatter = ggplot2::ggplot(sc, ggplot2::aes(x = pred, y = obs)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
        ggplot2::geom_point(alpha = 0.35, size = 1) +
        ggplot2::geom_text(data = sk, ggplot2::aes(label = label), x = -Inf, y = Inf,
                           hjust = -0.1, vjust = 1.4, size = 3.4, inherit.aes = FALSE) +
        fct + ggplot2::labs(x = "simulated", y = "observed") + theme_finn
    }
    # (2) mortality: the ROC curve of the interval death probability, with the AUC
    mo = op[variable == "mort"]
    if (nrow(mo)) {
      if (all(c("n_at_risk", "n_died") %in% names(mo)) && any(mo$n_at_risk > 0, na.rm = TRUE)) {
        grp = if (is.null(species)) NULL else "species_lab"
        roc = mo[is.finite(n_at_risk) & n_at_risk > 0,
                 roc_counts(p_interval, n_died, n_at_risk - n_died), by = grp]
        auc_lab = if (is.null(species)) sprintf("AUC = %.2f", skill[variable == "mort"]$auc)
                  else NULL
        panels$mort = ggplot2::ggplot(roc, ggplot2::aes(x = fpr, y = tpr)) +
          ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
          (if (is.null(species)) ggplot2::geom_line() else ggplot2::geom_line(ggplot2::aes(colour = species_lab))) +
          (if (!is.null(auc_lab)) ggplot2::annotate("text", x = 0.6, y = 0.1, label = auc_lab, size = 3.6)) +
          ggplot2::coord_equal() +
          ggplot2::labs(x = "false positive rate", y = "true positive rate", colour = NULL,
                        title = "mort: ROC of the interval death probability (trees at risk)",
                        subtitle = if (!is.null(species)) paste(sprintf("%s: AUC %.2f", labs[as.character(skill[variable == "mort"]$species)],
                                                                         skill[variable == "mort"]$auc), collapse = "; ")) +
          theme_finn + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
      } else {
        sk = skill[variable == "mort"]
        panels$mort = ggplot2::ggplot(mo, ggplot2::aes(x = pred, y = obs)) +
          ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
          ggplot2::geom_point(alpha = 0.35, size = 1) +
          ggplot2::labs(x = "simulated", y = "observed", title = sprintf("mort (no n_at_risk / n_died for an AUC): R\u00b2 = %.2f", sk$r2[1])) +
          theme_finn + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
      }
    }
    # (3) recruitment: simulated / observed total per species over all sites
    rg = op[variable == "reg"]
    if (nrow(rg)) {
      rs = rg[, list(obs = sum(obs), pred = sum(pred)), by = species]
      total_ratio = sum(rs$pred) / sum(rs$obs)
      rs = rs[obs > 0][order(-obs)]
      n_hidden = max(0L, nrow(rs) - n_species)
      rs = utils::head(rs, n_species)
      rs[, ratio := pred / obs]
      rs[, name := factor(labs[as.character(species)], levels = rev(labs[as.character(species)]))]
      rs[, share := obs / sum(rg$obs)]
      panels$reg = ggplot2::ggplot(rs, ggplot2::aes(x = ratio, y = name)) +
        ggplot2::geom_vline(xintercept = 1, colour = "grey55", linetype = "dashed") +
        ggplot2::geom_segment(ggplot2::aes(x = 1, xend = ratio, yend = name), colour = "grey70") +
        ggplot2::geom_point(ggplot2::aes(size = share)) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", ratio)), nudge_y = 0.35, size = 3) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_size_area(max_size = 4, guide = "none") +
        ggplot2::labs(x = "simulated / observed recruitment, all sites (log scale)", y = NULL,
                      title = sprintf("reg: total ratio %.2f%s", total_ratio,
                                      if (n_hidden) sprintf("; %d more species with less recruitment not shown", n_hidden) else ""),
                      subtitle = "point size: share of the observed recruitment") +
        theme_finn + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
    }
    p = structure(panels, class = "finn_comparison")
  } else {
    s = simulation_long(sim)
    s = s[variable != "reg"]                       # the drawn recruits; the expected ones are compared
    s[variable == "r_mean_ha", variable := "reg"]
    # the simulation over the same site-species pairs as the observations of
    # each response: a species absent from a site has dbh 0 in the simulation
    # and no observation, and must not pull the mean down
    pairs = unique(op[, list(siteID, species, variable = as.character(variable))])
    s = merge(s[variable %in% unique(pairs$variable)], pairs, by = c("siteID", "species", "variable"))
    ob = data.table::copy(op)
    # recruitment: observed per period, simulated per year
    ob[variable == "reg", obs := obs / period_length]
    by_sp = !is.null(species)
    if (!by_sp) {
      # all species pooled: stand totals for ba, trees and reg, the mean over species otherwise
      pool = function(v, x) if (v %in% c("ba", "trees", "reg")) sum(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
      sl = s[, list(value = pool(variable[1], value)), by = list(siteID, year, variable)][, list(value = mean(value)), by = list(year, variable)]
      ob = ob[, list(obs = pool(as.character(variable[1]), obs)), by = list(siteID, year, variable)][, list(obs = mean(obs)), by = list(year, variable)]
    } else {
      sl = s[, list(value = mean(value, na.rm = TRUE)), by = list(year, variable, species)]
      ob = ob[, list(obs = mean(obs, na.rm = TRUE)), by = list(year, variable, species)]
      sl[, species_lab := labs[as.character(species)]]
      ob[, species_lab := labs[as.character(species)]]
    }
    sl[, variable := factor(variable, levels = levels(op$variable))]
    ob[, variable := factor(as.character(variable), levels = levels(op$variable))]
    p = ggplot2::ggplot() +
      ggplot2::geom_line(data = sl, ggplot2::aes(x = year, y = value)) +
      ggplot2::geom_point(data = ob, ggplot2::aes(x = year, y = obs), colour = "#C0392B", size = 2) +
      facet +
      ggplot2::labs(x = "simulated year", y = if (by_sp) "mean over sites" else "mean over sites (ba, trees, reg: stand total)",
                    caption = "line: simulated; points: observed at the census (growth, mort, reg as annual rates)") +
      theme_finn
  }
  attr(p, "matched") = op
  attr(p, "skill") = skill
  print(p)
  invisible(p)
}

# ROC points of a probability against counts of positives and negatives that
# share it (every tree at risk one case).
roc_counts = function(p, pos, neg) {
  g = data.table::data.table(p = p, pos = pos, neg = neg)[, list(pos = sum(pos), neg = sum(neg)), by = p][order(-p)]
  if (sum(g$pos) == 0 || sum(g$neg) == 0) return(data.table::data.table(fpr = numeric(0), tpr = numeric(0)))
  data.table::data.table(fpr = c(0, cumsum(g$neg) / sum(g$neg)), tpr = c(0, cumsum(g$pos) / sum(g$pos)))
}

#' @rdname comparePlot
#' @export
print.finn_comparison = function(x, ...) {
  panels = unclass(x)[c("scatter", "mort", "reg")]
  panels = panels[!vapply(panels, is.null, logical(1))]
  if (!length(panels)) return(invisible(x))
  grid::grid.newpage()
  top = !is.null(panels$scatter)
  bottom = setdiff(names(panels), "scatter")
  rows = if (top && length(bottom)) c(1.1, 1) else 1
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(length(rows), max(1L, length(bottom)), heights = grid::unit(rows, "null"))))
  if (top) print(panels$scatter, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = seq_len(max(1L, length(bottom)))))
  for (i in seq_along(bottom))
    print(panels[[bottom[i]]], vp = grid::viewport(layout.pos.row = if (top) 2 else 1, layout.pos.col = i))
  grid::popViewport()
  invisible(x)
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
#' needed to judge it on the model object, as `model$evaluation`;
#' [`summary()`][summary.finn_class] prints it and [comparePlot()] draws the same
#' comparison. Calling [fit()] with its argument `holdout` runs it right after
#' the fit.
#'
#' * **Training data**: the model is simulated from the training initial state
#'   under the training environment and compared with the observations it was
#'   fitted to, aggregated the way [fit()] scores them (see [comparePlot()]).
#' * **Holdout** (`holdout = list(data =, env =, init_cohort =)`): the same on
#'   data the model has not seen, e.g. a spatially separate set of sites.
#' * **Bare ground** (`bare_ground = TRUE`): the run of
#'   `simulateForest(model, bare_ground = TRUE, years = , sites = )`, from an
#'   empty stand for `years` years under the mean climate of the fitted
#'   environment held constant, or under each site's mean climate with `sites`
#'   (see [simulateForest()]). Optional, because a long run takes time.
#'
#' @param model A model built by [finn()] and fitted with [fit()].
#' @param holdout `NULL` (default) or a list with the holdout `data`
#'   (observations as for [fit()]), its `env` and its `init_cohort`, and
#'   optionally `disturbance` (see [finn_data] for the four tables). Without
#'   `init_cohort` the holdout sites start from bare ground, with a warning.
#' @param bare_ground (`logical(1)`)\cr Run the bare-ground simulation too.
#'   Default `FALSE`.
#' @param years (`integer(1)`)\cr With `bare_ground = TRUE`: the length of the
#'   run in years, as in [simulateForest()]. Default `500`.
#' @param sites With `bare_ground = TRUE`, as in [simulateForest()]: `NULL`
#'   (default) one run under the mean climate of all sites, or site ids for
#'   one run per site under that site's mean climate.
#' @param patches (`integer(1)`)\cr Number of patches per site of a run that
#'   starts from bare ground; `NULL` (default) the patch count the model was
#'   fitted with (`model$patches`), or `100` for a model fitted before 0.2.0.
#'   Not used for a set with an `init_cohort`: its patch count comes from the
#'   initial cohorts.
#' @param train `NULL` (default): the training data the model keeps from
#'   [fit()]. A list like `holdout` for a model that keeps none (fitted before
#'   0.2.0, e.g. read with [torch::torch_load()]).
#' @return The model, invisibly; it is changed in place, as by [fit()], and
#'   carries `model$evaluation`: a list with `train` and, when given, `holdout`
#'   (each with the paired table `matched` and the skill measures `skill`
#'   overall and `skill_species` per species; the measures are described in
#'   [summary.finn_class()]), and with `bare_ground = TRUE` also `bare_ground`
#'   (the simulation, as [simulateForest()] returns it) and
#'   `bare_ground_settings` (`years`, `sites`, `patches`). A new call replaces
#'   the stored evaluation.
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
    # with initial cohorts the patch count is theirs; patches only for bare ground
    sim = if (is.null(init_cohort))
      predict(model, env = env, disturbance = disturbance, init_cohort = NULL, patches = patches)
    else predict(model, env = env, disturbance = disturbance, init_cohort = init_cohort)
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
