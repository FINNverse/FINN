
#' Resolve site, patch, and year indices for FINN inputs
#'
#' Joins tree, environment, and observation tables; assigns integer indices
#' for site, patch and year; standardizes species coding; and returns aligned
#' data plus optional initial cohorts.
#'
#' With \code{time_axis = "annual"} (the default) one simulation step is one
#' year: \code{year} becomes the number of years since the first inventory of
#' the site's patch (0 at the first inventory, e.g. 10 and 20 at the next two),
#' the first inventory is the initial state (the initial cohorts, not an
#' observation), the observations carry \code{period_length} in years, and
#' \code{env_dt} is returned with one row per site and year \code{1..Tmax},
#' \code{Tmax} being the longest time from a first to a last inventory. The
#' environment row of year \eqn{k} drives the simulated year that ends at
#' year \eqn{k}, so an inventory \eqn{k} years after the first is compared
#' with the stand after \eqn{k} simulated years. \code{env_dt} therefore needs
#' one row per site and calendar year from the year after the site's first
#' inventory to \code{Tmax} years after it; a climate that does not change is
#' repeated for every year.
#'
#' @param tree_dt data.table with tree-level data including siteName, patchName, year,
#'   species_name, dbh, status, living, and optional trees.
#' @param env_dt data.table with environment data including siteName and year
#'   (calendar year; see Details for the years needed).
#' @param obs_dt data.table with observations including siteName, patchName, year,
#'   species_name, and stand metrics (ba, dbh, trees, growth, mort, n_at_risk,
#'   n_died, reg, period_length), as returned by \code{makeObsData}.
#' @param createInitCohorts logical. If TRUE, build FINN initial cohorts from
#'   the living trees of the first inventory (\code{year == 0}; with
#'   \code{time_axis = "census"}, as in FINN 0.1.0, from \code{year == 1}).
#'   Default TRUE.
#' @param time_axis \code{"annual"} (default) or \code{"census"}.
#'   \code{"census"} is the axis of FINN 0.1.0: \code{year} is the index of
#'   the inventory (0, 1, 2, ...), one simulation step is one census interval
#'   (\code{period_length = 1}), \code{env_dt} keeps the inventory years, and
#'   the observations of the first inventory stay in \code{obs_dt}. Use it
#'   with \code{makeObsData(rates = "interval")}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{siteID_dt}: site/patch/year index map (\code{OrigYear} is the
#'     calendar year).
#'   \item \code{tree_dt}: tree table with indices and standardized species.
#'   \item \code{env_dt}: environment table with indices, the years of the
#'     simulation.
#'   \item \code{obs_dt}: site-aggregated observations by species, with
#'     \code{period_length}.
#'   \item \code{obs_dt_patches}: patch-level observations by species.
#'   \item \code{species_dt}: species lookup (\code{species}, \code{species_name}).
#'   \item \code{initCohorts} (optional): FINN \code{CohortMat} object.
#'   \item \code{initCohort_dt} (optional): trees used to build cohorts.
#' }
#'
#' @examples
#' library(data.table)
#' # eight Oregon FIA plots, three inventories each; see the data-preparation vignette
#' tree_dt <- fread(system.file("extdata", "example_tree_dt.csv", package = "FINN"))
#' env_dt  <- fread(system.file("extdata", "example_env_dt.csv",  package = "FINN"))
#' obs <- makeObsData(tree_dt[complete == TRUE], plotsize = 0.06, aggregate_by_site = FALSE,
#'                    fix_period_length = 10, NspeciesQuantile = 0.98)
#'
#' # Names become integer ids shared by every table, calendar years become years
#' # since the first inventory: 0 (the initial state), 10 and 20. Pass the tree
#' # list makeObsData() returned: it carries the same species lumping as obs_dt.
#' # createInitCohorts = FALSE leaves out the starting cohorts, which need torch.
#' inputs <- resolveSiteIDs(tree_dt = obs$tree_dt, obs_dt = obs$obs_dt, env_dt = env_dt,
#'                          createInitCohorts = FALSE)
#' inputs$species_dt
#' head(inputs$obs_dt)      # per site, year and species: what fit(data = ) takes
#' range(inputs$env_dt$year) # one row per site and simulated year: what fit(env = ) takes
#' head(inputs$siteID_dt)   # the map back to the inventory's own names
#' @examplesIf torch::torch_is_installed()
#' # with torch, the same call also builds the starting cohorts (a CohortMat)
#' inputs <- resolveSiteIDs(tree_dt = obs$tree_dt, obs_dt = obs$obs_dt, env_dt = env_dt)
#' inputs$initCohorts$dims  # sites, patches, cohorts
#'
#' @import data.table
#' @export
resolveSiteIDs <- function(tree_dt, env_dt, obs_dt, createInitCohorts = TRUE, time_axis = c("annual", "census")){
  time_axis <- match.arg(time_axis)
  annual <- time_axis == "annual"
  tree_dt <- copy(as.data.table(tree_dt))
  env_dt  <- copy(as.data.table(env_dt))
  obs_dt  <- copy(as.data.table(obs_dt))
  siteID_dt <- unique(env_dt[,.(siteName, year)])
  siteID_dt <- merge(siteID_dt, unique(tree_dt[,.(siteName, year)]), by = c("siteName", "year"), all = FALSE)
  siteID_dt <- merge(siteID_dt, unique(obs_dt[,.(siteName, patchName, year)]), by = c("siteName", "year"), all = FALSE)

  # Radix, not the locale's collation: `as.factor()` orders through LC_COLLATE, so the
  # same file gives one numbering under C and another under de_DE, and the numbering
  # indexes the fitted parameters.
  siteID_dt[, siteID := match(siteName, sort(unique(siteName), method = "radix"))]
  siteID_dt[, patchID := match(patchName, sort(unique(patchName), method = "radix")), by = siteID]
  siteID_dt[, period := as.integer(as.factor(year)), by = .(siteID, patchID)]
  # The simulated year of each inventory: years since the patch's first
  # inventory on the annual axis, the inventory index on the census axis.
  if (annual) siteID_dt[, step := as.integer(year - min(year)), by = .(siteID, patchID)]
  else        siteID_dt[, step := period - 1L]

  tree_dt <- merge(tree_dt, siteID_dt, by = c("siteName", "patchName", "year"), all = FALSE)

  species_levels <- sort(unique(tree_dt$species_name), method = "radix")
  species_levels <- c(species_levels[species_levels != "other"], "other")

  tree_dt[,species := as.integer(factor(species_name, levels = species_levels, labels = 1:length(species_levels))),]
  species_dt <- unique(tree_dt[,.(species,species_name)])

  tree_dt[, OrigYear := year,]
  tree_dt[, year := step,]
  tree_dt <- tree_dt[, -c("period", "step")]

  # siteID_dt has one row per site, patch and year; the environment is a site
  # property, so take the site x year rows once or every patch would duplicate
  # them and climateDF2array() refuses the table when fit() builds its arrays.
  if (annual) {
    # One row per site and simulated year 1..Tmax: the row of year k drives the
    # year that ends k years after the first inventory, so the inventory at
    # offset k is scored against env row k, and there is no row 0.
    site_first <- siteID_dt[, .(first = min(year)), by = .(siteName, siteID)]
    Tmax <- max(siteID_dt$step)
    if (Tmax < 1) stop("resolveSiteIDs(): no site has an inventory after its first one, so there is nothing to simulate.", call. = FALSE)
    need <- site_first[, .(step = seq_len(Tmax), year = first + seq_len(Tmax)), by = .(siteName, siteID)]
    env_dt[, .env_row := TRUE]
    env_dt <- merge(need, env_dt, by = c("siteName", "year"), all.x = TRUE)
    if (anyNA(env_dt$.env_row)) {
      miss <- env_dt[is.na(.env_row)]
      stop(sprintf(paste0("resolveSiteIDs(): `env_dt` has no row for %d site-year(s) of the annual time axis, e.g. site %s in %s. ",
                          "It needs one row per site and calendar year from the year after the site's first inventory to %d years after it ",
                          "(a climate that does not change is repeated for every year), or use time_axis = \"census\"."),
                   nrow(miss), miss$siteName[1], format(miss$year[1]), Tmax), call. = FALSE)
    }
    if (anyDuplicated(env_dt[, .(siteID, step)]))
      stop("resolveSiteIDs(): `env_dt` has more than one row for a site and year.", call. = FALSE)
    env_dt[, .env_row := NULL]
  } else {
    env_dt <- merge(env_dt, unique(siteID_dt[,.(siteName, siteID, year, step)]), by = c("siteName","year"), all = FALSE)
  }
  env_dt[, OrigYear := year,]
  env_dt[, year := step,]
  env_dt <- env_dt[, -"step"]
  data.table::setcolorder(env_dt, c("siteName", "year"))
  data.table::setorder(env_dt, siteID, year)

  # The unit of growth and mort in obs_dt (makeObsData(rates = )), read off
  # the data: where a cohort died partly over an interval of several years the
  # yearly and the interval rate differ. Elsewhere the two are the same number.
  if ("period_length" %in% names(obs_dt)) {
    chk <- obs_dt[is.finite(mort) & n_died > 0 & n_died < n_at_risk & is.finite(period_length) & period_length > 1]
    rates_in <- if (nrow(chk) && isTRUE(all(abs(chk$mort - obs_mort(chk$n_died, chk$n_at_risk, chk$period_length, "annual")) < 1e-8)))
      "annual" else "interval"
    setnames(obs_dt, "period_length", "period_length_in")
  } else {
    if (annual) warning("resolveSiteIDs(): `obs_dt` has no `period_length` column (it was not built by makeObsData() of FINN 0.2.0), so its growth and mortality are taken as rates over the whole interval, while the time axis is annual. Rebuild it with makeObsData().", call. = FALSE)
    rates_in <- "interval"
    obs_dt[, period_length_in := NA_real_]
  }
  if (!annual && rates_in == "annual")
    warning("resolveSiteIDs(): time_axis = \"census\" (one step per census interval) with yearly growth and mortality rates; build obs_dt with makeObsData(rates = \"interval\").", call. = FALSE)

  obs_dt <- merge(obs_dt, siteID_dt, by = c("siteName", "patchName", "year"), all = FALSE)
  obs_dt[, OrigYear := year,]
  obs_dt[, year := step,]
  if (annual) {
    # the interval each inventory closes, in years; the first inventory is the
    # initial state, not an observation
    steps <- unique(obs_dt[, .(siteID, patchID, year)])
    data.table::setorder(steps, siteID, patchID, year)
    steps[, period_length := year - data.table::shift(year, 1, type = "lag"), by = .(siteID, patchID)]
    obs_dt <- merge(obs_dt, steps, by = c("siteID", "patchID", "year"))
    obs_dt <- obs_dt[year > 0]
    obs_dt[, pl_rate := period_length]
  } else {
    obs_dt[, pl_rate := period_length_in]
    obs_dt[, period_length := data.table::fifelse(year > 0, 1, NA_real_)]
  }
  obs_dt <- obs_dt[, -c("period", "step", "period_length_in")]
  obs_dt <- merge(obs_dt, species_dt, by = "species_name", all = TRUE)
  out_cols_obs_dt <- c("siteID", "patchID", "year", "ba", "dbh", "trees", "growth", "growth_n", "mort", "n_at_risk", "n_died", "reg", "species", "species_name", "period_length")

  obs_dt_aggr <- obs_dt[,.(ba = mean(ba, na.rm = TRUE),
                        dbh = mean(dbh, na.rm = TRUE),
                        trees = mean(trees, na.rm = TRUE),
                        growth = stats::weighted.mean(growth, growth_n, na.rm = TRUE),
                        growth_n = sum(growth_n, na.rm = TRUE),
                        # mortality counts pool by summing; the rate is derived
                        # from the pooled counts, never averaged from rates.
                        n_at_risk = sum(n_at_risk, na.rm = TRUE),
                        n_died = sum(n_died, na.rm = TRUE),
                        reg = mean(reg, na.rm = TRUE),
                        period_length = max(period_length),
                        pl_rate = max(pl_rate)),
                    by = .(siteID, year, species, species_name)]
  # in the unit the observations came in (see rates_in above)
  obs_dt_aggr[, mort := obs_mort(n_died, n_at_risk, pl_rate, rates_in)]
  out_cols_obs_dt_aggr <- c("siteID", "year", "ba", "dbh", "trees", "growth", "growth_n", "mort", "n_at_risk", "n_died", "reg", "species", "species_name", "period_length")

  siteID_dt[, OrigYear := year, by = .(siteID, patchID)]
  if (annual) siteID_dt[, year := step]
  else        siteID_dt[, year := period, by = .(siteID, patchID)]
  out_cols_siteID_dt <- c("siteID", "patchID", "siteName", "patchName", "year", "OrigYear")

  out <- list(
    siteID_dt = as.data.table(siteID_dt[,..out_cols_siteID_dt], key = NULL),
    tree_dt = as.data.table(tree_dt, key = NULL),
    env_dt = as.data.table(env_dt, key = NULL),
    obs_dt = as.data.table(obs_dt_aggr[,..out_cols_obs_dt_aggr]),
    obs_dt_patches = as.data.table(obs_dt[,..out_cols_obs_dt]),
    species_dt = species_dt
    # initCohort_dt = as.data.table(initCohort_dt, key = NULL)
  )

  if(createInitCohorts){
    init_trees <- tree_dt[living == TRUE & year == (if (annual) 0 else 1)]
    if(!any(colnames(init_trees) %in% "trees")) init_trees[, trees := 1,]
    out[["initCohorts"]] <- makeInitCohorts(init_trees, Nspecies = uniqueN(obs_dt$species), treeID_table = FALSE)
    out[["initCohort_dt"]] <- init_trees
  }
  return(out)
}


# makeEnvData <- function(){
#
#   return(env_dt)
# }

#' Convert DBH to basal area
#'
#' Computes basal area in m^2 from diameter at breast height in cm:
#' \deqn{BA = \pi \left(\frac{\mathrm{DBH}}{200}\right)^2}
#'
#' @param dbh numeric vector of diameters (cm).
#' @return numeric vector of basal areas (m^2).
#' @export
dbh2ba <- function(dbh){
  return((dbh/200)^2 * pi)
}

# The mortality rate of a closed cohort from its counts. `rates = "interval"`:
# the fraction n_died / n_at_risk over the census interval. `rates = "annual"`:
# the constant yearly probability that compounds to that fraction over
# `period_length` years, 1 - (1 - n_died / n_at_risk)^(1 / period_length). NA
# where no tree was at risk; an interval without a length is left as it is.
obs_mort <- function(n_died, n_at_risk, period_length, rates = c("annual", "interval")) {
  rates <- match.arg(rates)
  m <- data.table::fifelse(n_at_risk > 0, n_died / n_at_risk, NA_real_)
  if (rates == "annual") {
    p <- as.numeric(period_length)
    ok <- !is.na(m) & is.finite(p) & p > 0
    m[ok] <- 1 - (1 - m[ok])^(1 / p[ok])
  }
  m
}

#' Create observation data from trees
#'
#' Derives stand metrics by site/patch/species/year from tree measurements,
#' filters sites, harmonizes years, optionally aggregates by site.
#'
#' @param tree_dt data.table of tree records with siteName, patchName, year,
#'   treeName, species_name, dbh, \code{status} (one of "alive", "new", "dead")
#'   and \code{living}. Mortality is derived from \code{status}; a \code{mort}
#'   column, if present, is ignored.
#' @param plotsize numeric plot area used to scale recruitment.
#' @param aggregate_by_site logical. Aggregate patches to site level. Default TRUE.
#' @param minNyears integer or NULL. Keep sites with at least this many years and
#'   equal counts across patches. Default 2.
#' @param fix_period_length integer or NULL. If set, drop sites whose inventory
#'   interval differs. Default NULL.
#' @param dbh_growth_thresh length-2 numeric or NULL. Drop sites where any tree
#'   \code{dbh} change falls outside this range. Default \code{c(-10, 50)}.
#' @param Npatches integer or NULL. If set, keep only sites with exactly this many patches.
#' @param Nspecies integer or NULL. Cap species to the top N (others merged to "other").
#' @param NspeciesQuantile numeric in (0,1] or NULL. Choose the smallest N covering
#'   this fraction of individuals; overrides \code{Nspecies} if supplied.
#' @param rates \code{"annual"} (default) or \code{"interval"}: the unit of
#'   \code{growth} and \code{mort}. \code{"annual"} returns them per year, for a
#'   model simulated in yearly steps: a tree's growth over an interval of
#'   \eqn{p} years is \eqn{(dbh/dbh_{before})^{1/p} - 1}, and the mortality of a
#'   cohort is \eqn{1 - (1 - n_{died}/n_{at\,risk})^{1/p}}, the constant yearly
#'   rates that compound to what was observed over the interval.
#'   \code{"interval"} returns them over the whole interval
#'   (\code{dbh/dbh_before - 1} and \code{n_died / n_at_risk}), the output of
#'   FINN 0.1.0, for a model in which one simulation step is one census
#'   interval. \eqn{p} is \code{year - year_before} in the units of
#'   \code{year}, so \code{year} must be the calendar year of the inventory.
#'
#' @details
#' FINN simulates the forest in steps of one year. With \code{rates = "annual"}
#' the observation table is in the units that \code{fit()} compares with such a
#' simulation: \code{fit()} averages the simulated yearly growth and mortality
#' over the \code{period_length} years before each census and sums the
#' simulated recruits over them (see \code{\link{finn_data}}). The counts
#' (\code{n_at_risk}, \code{n_died}, \code{growth_n}) and \code{reg} stay
#' quantities of the whole interval in both modes: the counts weight the
#' likelihood by the number of trees behind an observation, and \code{reg} is
#' compared with the recruits summed over the interval. To use the table,
#' express \code{year} as years since the first inventory, give \code{env} one
#' row per site and year from 1 to the last census, and build the initial
#' cohorts from the first inventory; \code{\link{resolveSiteIDs}} does all
#' three.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{obs_dt}: observations at site or patch level, one row per site
#'     (and patch), inventory \code{year} and species. \code{period_length} is
#'     the number of years since the site's previous inventory (\code{NA} at
#'     its first). \code{growth} is the mean relative diameter growth rate of
#'     the surviving trees (per year, or over the interval, see \code{rates}),
#'     and \code{growth_n} is how many trees that mean rests on. Both matter: means
#'     are aggregated across patches TREE-WEIGHTED (by \code{growth_n}), because
#'     that is the quantity the model predicts - it simulates
#'     \code{sum(g*trees)/sum(trees)} over the whole site, so an unweighted mean
#'     of patch means would be a different number. \code{growth_n} also travels
#'     with the response so \code{fit} can weight the likelihood by it: 34% of
#'     patch-level growth observations rest on a single tree while others average
#'     30+, and the variance of a mean is \eqn{\sigma^2/n}. Mortality comes
#'     back as a closed-cohort pair of counts, \code{n_at_risk} (trees alive at
#'     the start of the interval) and \code{n_died} (how many of them were dead
#'     at the end), plus the derived rate \code{mort} (\code{NA} where no
#'     cohort was at risk): the yearly rate
#'     \code{1 - (1 - n_died / n_at_risk)^(1 / period_length)}, or
#'     \code{n_died / n_at_risk} with \code{rates = "interval"}. \code{fit()}
#'     scores \code{mort} with the binomial loss weighted by \code{n_at_risk}.
#'     \code{reg} is the number of recruits per hectare over the interval.
#'   \item \code{tree_dt}: input trees with added growth fields and species
#'     recode: \code{year_before}, \code{dbh_before}, \code{period_length}
#'     (years since the tree's previous record), \code{rel_growth} (the
#'     increment over the interval, \code{dbh/dbh_before - 1}) and
#'     \code{growth} (the tree's growth rate in the unit of \code{rates}, the
#'     quantity behind \code{obs_dt$growth} and the \code{growth} column of
#'     \code{fit(tree_data = )}); an input column \code{growth} is replaced.
#' }
#'
#' @examples
#' library(data.table)
#' # A tree list of eight Oregon FIA plots with four patches and three inventories
#' # each, one row per tree and inventory: siteName, patchName, treeName, year,
#' # species_name, dbh (cm), status ("new", "alive", "dead") and living.
#' tree_dt <- fread(system.file("extdata", "example_tree_dt.csv", package = "FINN"))
#' head(tree_dt[, .(siteName, patchName, treeName, year, species_name, dbh, status)])
#'
#' obs <- makeObsData(
#'   tree_dt[complete == TRUE],    # fully remeasured plots only
#'   plotsize          = 0.06,     # ha per patch; recruits are scaled to per hectare
#'   aggregate_by_site = FALSE,    # keep the patches: resolveSiteIDs() needs them
#'   fix_period_length = 10,       # drop plots whose remeasurement interval is not 10 years
#'   NspeciesQuantile  = 0.98)     # species covering 98 % of the stems; the rest is "other"
#' head(obs$obs_dt)
#'
#' # mortality is a closed cohort: the trees alive at the start of the interval
#' # and how many of them had died by its end; `mort` is the yearly rate over the
#' # period_length years of the interval
#' obs$obs_dt[n_at_risk > 0, .(siteName, patchName, year, period_length, species_name,
#'                             n_at_risk, n_died, mort)][1:5]
#' # growth is a yearly rate too, and travels with the number of trees behind it
#' obs$obs_dt[growth_n > 0, .(siteName, year, species_name, growth, growth_n)][1:5]
#'
#' # the same table with the rates over the whole ten-year interval (FINN 0.1.0)
#' obs10 <- makeObsData(tree_dt[complete == TRUE], plotsize = 0.06, aggregate_by_site = FALSE,
#'                      fix_period_length = 10, NspeciesQuantile = 0.98, rates = "interval")
#' obs10$obs_dt[growth_n > 0, .(siteName, year, species_name, growth, growth_n)][1:5]
#'
#' @import data.table
#' @export
makeObsData <- function(tree_dt, plotsize, aggregate_by_site = TRUE, minNyears = 2, fix_period_length = NULL, dbh_growth_thresh = c(-10,50), Npatches = NULL, Nspecies = NULL, NspeciesQuantile = NULL, rates = c("annual", "interval")){
  rates <- match.arg(rates)
  tree_dt <- copy(as.data.table(tree_dt))
  data.table::setorder(tree_dt, treeName, year)
  tree_dt[, ":="(
    year_before = data.table::shift(year,1,type = "lag"),
    dbh_before = data.table::shift(dbh,1,type = "lag"),
    # state at the START of the interval ending at this record - the closed
    # cohort below is defined entirely in terms of these.
    status_before  = data.table::shift(status,1,type = "lag"),
    living_before  = data.table::shift(living,1,type = "lag"),
    species_before = data.table::shift(species_name,1,type = "lag")
    ), by = treeName]
  tree_dt[, period_length := year - year_before,]
  # exclude all plots without consistent period length
  tree_dt[, dbh_growth := dbh - dbh_before,]
  tree_dt[status == "alive", rel_growth := dbh/dbh_before - 1,]
  # `growth` is the per-tree rate behind obs_dt$growth: annual by default, so a
  # census interval of p years gives (dbh/dbh_before)^(1/p) - 1, the constant
  # yearly rate that compounds to the observed increment. That is what fit()
  # compares with the mean of the simulated yearly rates over the interval.
  # `rel_growth` stays the increment over the whole interval.
  tree_dt[, growth := NA_real_]
  if (rates == "annual") {
    tree_dt[status == "alive" & is.finite(period_length) & period_length > 0,
            growth := (dbh/dbh_before)^(1/period_length) - 1]
  } else {
    tree_dt[status == "alive", growth := rel_growth]
  }
  tree_dt[living == TRUE,trees := 1,]

  if(!is.null(fix_period_length)){
    excluded_sites <- unique(tree_dt[period_length != fix_period_length]$siteName)
    tree_dt <- tree_dt[!(siteName %in% excluded_sites)]
    message(paste0("excluded ", length(excluded_sites)," sites due to period length not beeing ", fix_period_length))
  }
  if(!is.null(dbh_growth_thresh)){
    excluded_sites <- unique(tree_dt[dbh_growth < dbh_growth_thresh[1] | dbh_growth > dbh_growth_thresh[2]]$siteName)
    tree_dt <- tree_dt[!(siteName %in% excluded_sites)]
    message(paste0("excluded ", length(excluded_sites)," sites due to dbh_growth being not within the range of ", dbh_growth_thresh[1], " and ", dbh_growth_thresh[2]))
  }

  if(!is.null(Nspecies)){
    alive_species <- names(sort(table(tree_dt[status %chin% c("alive","new")]$species_name), decreasing = TRUE))
    selected_species <- alive_species[seq_len(min(Nspecies, length(alive_species)))]
    message(paste0("Nspecies set to ", Nspecies))
  }else if(!is.null(NspeciesQuantile)){
    species_table = sort(table(tree_dt[status %chin% c("alive","new")]$species_name), decreasing = TRUE)
    cumsum_species = cumsum(species_table)/sum(species_table)
    Nspecies = which(cumsum_species >= NspeciesQuantile)[1]
    selected_species <- names(cumsum_species)[1:Nspecies]
    message(paste0("Nspecies set to ", Nspecies, " to cover ", round(NspeciesQuantile*100,1), "% of individuals"))
  }else{
    Nspecies = uniqueN(tree_dt[status %chin% c("alive","new")]$species_name)
    selected_species = unique(tree_dt[status %chin% c("alive","new")]$species_name)
    message(paste0("Nspecies set to all species (", Nspecies, ")"))
  }
  message(paste0("Selected species: ", paste(selected_species, collapse = ", ")))

  other_species <- unique(tree_dt[!(species_name %in% selected_species)]$species_name)
  tree_dt <- tree_dt[!(species_name %in% selected_species), species_name := "other",]
  Nspecies = Nspecies + 1
  message(paste0("Combined ", uniqueN(other_species), " unselected species into 'other'\nactuall number of species including 'others' is ", Nspecies))

  obs_dt <- tree_dt[,.(
    ba = sum(dbh2ba(dbh)*(living == TRUE), na.rm = TRUE),
    growth = mean(growth[living == TRUE], na.rm = TRUE),
    # `growth_n` is how many trees the growth mean rests on. 34% of FIA growth
    # observations come from a SINGLE tree while others average 30+, and the
    # variance of a mean is sigma^2/n - so the count has to travel with the
    # response for the loss to weight it (as n_at_risk does for mortality).
    growth_n = sum(is.finite(growth) & living == TRUE, na.rm = TRUE),
    dbh = mean(dbh[living == TRUE], na.rm = TRUE),
    trees = sum(living == TRUE, na.rm = TRUE),
    reg = sum(status == "new", na.rm = TRUE)/plotsize
  ),by= .(siteName,patchName, year, species_name)]

  ## --- mortality: a CLOSED COHORT of the trees alive at the interval start ----
  ## Mortality is reported as two counts, `n_at_risk` and `n_died`, rather than
  ## as a bare rate. That is the `cbind(died, survived)` response of a binomial
  ## GLM, and it is what the "binomial" loss consumes: the counts carry the
  ## sample size, so an observation is weighted by how many trees it actually
  ## summarises, and they aggregate by SUMMING (below, and across species when
  ## lumping) instead of by an unweighted mean of per-patch proportions.
  ##
  ## Both columns are pinned to the tree's state at the START of the interval
  ## (`living_before` / `species_before`). This matters: FIA re-identifies a
  ## tree's species between visits, so a rate built from a patch x species count
  ## of survivors could book a death against the species the tree ENDED as while
  ## the denominator counted the species it STARTED as - which is how the old
  ## `n_mort / trees_before` produced mort > 1. Pinning both to the start makes
  ## `n_died <= n_at_risk` true by construction.
  ##
  ## NOTE - this is a CLOSED cohort: trees that recruit during the interval
  ## appear in neither column, so recruit mortality is not scored. The
  ## alternative is an OPEN cohort (at_risk = alive at start + recruits during
  ## the interval; died = every death incl. those recruits). That does score
  ## recruit mortality, at the cost of biasing the rate downwards, because a
  ## recruit is not exposed for the whole interval. To switch, admit the "new"
  ## trees into both columns here and carry an exposure offset into the loss.
  mort_cohort <- tree_dt[living_before == TRUE, .(
    n_at_risk = .N,
    n_died    = sum(status == "dead", na.rm = TRUE)
  ), by = .(siteName, patchName, year, species_name = species_before)]

  obs_dt <- merge(obs_dt, mort_cohort,
                  by = c("siteName", "patchName", "year", "species_name"), all.x = TRUE)

  # 1) Build, per site, the full cartesian grid of patch × year × species
  grid <- obs_dt[, CJ(
    patchName    = unique(patchName),
    year         = unique(year),
    species_name = unique(obs_dt$species_name),
    unique = TRUE
  ), by = siteName]

  # 2) Key both tables on the join columns
  setkey(obs_dt,  siteName, patchName, year, species_name)
  setkey(grid,    siteName, patchName, year, species_name)

  # 3) Join existing data onto the full grid (fills non-matches with NA)
  obs_dt <- obs_dt[grid]

  obs_dt[is.na(ba), ba := 0]
  obs_dt[is.na(trees), trees := 0]
  obs_dt[is.na(reg), reg := 0]
  obs_dt[is.na(growth), growth := NA_real_]
  obs_dt[is.na(growth_n), growth_n := 0L]
  obs_dt[is.na(dbh), dbh := NA_real_]
  # No cohort at risk -> the rate is undefined, not zero. Keep the counts at 0 so
  # they still sum correctly on aggregation, and let `mort` carry the NA so the
  # loss masks the observation out.
  obs_dt[is.na(n_at_risk), n_at_risk := 0L]
  obs_dt[is.na(n_died),    n_died    := 0L]
  # The interval each row closes, in years: the time since the site's previous
  # inventory (NA at its first). It is the window fit() aggregates the simulated
  # yearly rates and recruits over, and the exponent of the annualisation.
  site_periods <- unique(obs_dt[, .(siteName, year)])
  data.table::setorder(site_periods, siteName, year)
  site_periods[, period_length := year - data.table::shift(year, 1, type = "lag"), by = siteName]
  obs_dt <- merge(obs_dt, site_periods, by = c("siteName", "year"), all.x = TRUE)
  obs_dt[, mort := obs_mort(n_died, n_at_risk, period_length, rates)]
  obs_dt = obs_dt[order(siteName, patchName, species_name, year)]

  if(!is.null(minNyears)){
    obs_dt[, NyearsPerPatch := uniqueN(year), by = .(siteName, patchName)]
    obs_dt[, sameYearsAllPatches := all(NyearsPerPatch == max(NyearsPerPatch)), by = siteName]
    obs_dt <- obs_dt[sameYearsAllPatches == TRUE & minNyears <= NyearsPerPatch]
    message(paste0("Filtered to sites with at least ", minNyears, " years and all patches having the same number of years:\n", paste0("Sites with ",names(table(obs_dt$NyearsPerPatch)), " Inventories: ",table(obs_dt$NyearsPerPatch), collapse = "\n")))
  }

  if(!is.null(Npatches)){
    obs_dt[, NpatchesPerSite := uniqueN(patchName), by = siteName]
    obs_dt <- obs_dt[NpatchesPerSite == Npatches]
    message(paste0("Filtered to sites with exactly ", Npatches, " patches: ", table(obs_dt$NpatchesPerSite)))
  }else{
    obs_dt[, NpatchesPerSite := uniqueN(patchName), by = siteName]
    message(paste0("Sites have between ", min(obs_dt$NpatchesPerSite), " and ", max(obs_dt$NpatchesPerSite), " patches"))
  }

  if(aggregate_by_site){
    obs_dt <- obs_dt[,.(
      ba = mean(ba, na.rm = TRUE),
      trees = mean(trees, na.rm = TRUE),
      dbh = mean(dbh, na.rm = TRUE),
      # Tree-weighted, NOT a plain mean of patch means: the model predicts
      # sum(g*trees)/sum(trees) over the whole site, so an unweighted mean of
      # patch means is comparing against a different quantity (and lets a
      # 1-tree patch count as much as a 30-tree one).
      growth = stats::weighted.mean(growth, growth_n, na.rm = TRUE),
      growth_n = sum(growth_n, na.rm = TRUE),
      # counts SUM across patches; the site rate is then the pooled rate. The
      # old `mean(mort)` weighted a 1-tree patch the same as a 50-tree one.
      n_at_risk = sum(n_at_risk, na.rm = TRUE),
      n_died = sum(n_died, na.rm = TRUE),
      reg = mean(reg, na.rm = TRUE),
      Npatches = uniqueN(patchName)
    ), by = .(siteName, year, species_name, period_length)]
    obs_dt[, mort := obs_mort(n_died, n_at_risk, period_length, rates)]
    data.table::setcolorder(obs_dt, setdiff(names(obs_dt), "period_length"))   # last, as in the patch table
    message("Aggregated obs_dt by siteName. For unaggregated data, set aggregate_by_site = FALSE")
  }else{
    obs_dt <- obs_dt[,.(siteName, patchName, year, species_name, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, period_length)]
    message("Kept obs_dt unaggregated by siteName. For aggregated data, set aggregate_by_site = TRUE")
  }

  return(list(obs_dt = obs_dt, tree_dt = tree_dt))
}

#' Make initial cohorts for FINN
#'
#' Aggregates initial trees into cohorts by DBH bins or exact DBH and
#' constructs a \code{FINN::CohortMat} object.
#'
#' @description
#' This function prepares \code{obs_df} in the exact format expected by
#' \code{FINN::CohortMat$new()} and calls it. The required long-table schema is:
#' \itemize{
#'   \item \strong{siteID} integer index of sites \eqn{1..S}.
#'   \item \strong{patchID} integer index of patches within site \eqn{1..P_s}.
#'   \item \strong{species} integer species code \eqn{1..sp}.
#'   \item \strong{dbh} numeric DBH in cm (either exact or binned midpoints).
#'   \item \strong{trees} integer count of trees in the cohort.
#' }
#'
#' @details
#' Internally calls:
#' \preformatted{
#' FINN::CohortMat$new(
#'   obs_df = <data.frame with columns siteID, patchID, species, dbh, trees>,
#'   dbh    = NULL,
#'   trees  = NULL,
#'   species= NULL,
#'   dims   = c(S, P, K),   # inferred from obs_df
#'   sp     = Nspecies,     # passed from argument
#'   device = "cpu"
#' )
#' }
#' Key fields of the resulting R6 object:
#' \itemize{
#'   \item \code{dbh}, \code{trees}, \code{species}: tensors per cohort.
#'   \item \code{dims}: integer vector \code{c(S, P, K)} for sites, patches, cohorts.
#'   \item \code{sp}: integer number of species.
#'   \item \code{device}: "cpu" or "cuda".
#'   \item \code{dbh_r}, \code{trees_r}, \code{species_r}: R arrays.
#'   \item \code{obsDF2arrays}: method converting \code{obs_df} into arrays.
#' }
#'
#' @param init_trees data.table of initial trees with columns
#'   \code{siteID}, \code{patchID}, \code{species}, \code{dbh}, \code{treeName};
#'   optional \code{trees} for pre-counted individuals.
#' @param dbh_binsize numeric or NULL. Bin width in cm. If NULL, keep exact DBH.
#' @param min_dbh numeric or NULL. Lower bound for binning. If NULL, uses min DBH.
#' @param Nspecies integer. Number of species levels passed to \code{sp}.
#' @param treeID_table logical. If TRUE, also return the cohort table used.
#' @param singleCohortTreeNames character vector or NULL. Tree names excluded from
#'   binning and kept as single-tree cohorts.
#'
#' @return If \code{treeID_table = FALSE}, a \code{FINN::CohortMat} object.
#'   If \code{TRUE}, a list with:
#'   \itemize{
#'     \item \code{initCohort}: the \code{CohortMat} object.
#'     \item \code{init_trees}: the \code{obs_df} passed to \code{CohortMat$new()}.
#'   }
#'
#' @seealso \code{\link[FINN]{CohortMat}}
#'
#' @import data.table
#' @export
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # one cohort per observed tree, which is what `loss_resolution = "tree"` needs
#' init <- makeInitCohorts(fia_init_trees[siteID %in% 1:3],
#'                         Nspecies = max(fia_obs_dt$species))
#' dim(init$dbh)   # sites x patches x cohorts
#'
#' # 5 cm diameter bins instead: fewer cohorts, so a much cheaper simulation
#' binned <- makeInitCohorts(fia_init_trees[siteID %in% 1:3], dbh_binsize = 5,
#'                           Nspecies = max(fia_obs_dt$species))
#' dim(binned$dbh)
makeInitCohorts <- function(init_trees, dbh_binsize = NULL, min_dbh = NULL, Nspecies, treeID_table = FALSE, singleCohortTreeNames = NULL){
  # browser()
  if(!is.null(dbh_binsize)){
    if(is.null(min_dbh)) min_dbh = min(init_trees$dbh, na.rm = TRUE)
    dbh_intervals = seq(min_dbh, max(init_trees$dbh, na.rm = TRUE) + dbh_binsize, by = dbh_binsize)
    init_trees <- init_trees[!(treeName %in% singleCohortTreeNames),.(
      trees = sum(trees)),
      by = .(siteID, patchID, species,
             dbh = as.numeric(as.character(
               cut(dbh, breaks = dbh_intervals, labels = dbh_intervals[-1]-dbh_binsize/2, include.lowest = TRUE)

               ))
             )
      ]
  }

  init_trees[,cohortID := 1:.N, by = .(siteID, patchID)]

  initCohort <- FINN::CohortMat$new(obs_df = init_trees, sp = Nspecies)

  if(treeID_table) out = list(initCohort = initCohort, init_trees = init_trees) else if(!treeID_table) out = initCohort
  return(out)
}

#
# makeInitCohorts2 <- function(
#     init_trees,
#     dbh_bins,                   # numeric vector (same length as tree_groups)
#     tree_groups,                # list of character vectors (tree names)
#     Nspecies,
#     min_dbh = NULL,             # optional: scalar lower bound for all binning
#     return_tree_table = FALSE   # formerly treeID_table
# ){
#   stopifnot(is.list(tree_groups), length(dbh_bins) == length(tree_groups))
#
#   # ensure data.table
#   init_trees <- data.table::as.data.table(init_trees)
#
#   required_cols <- c("siteID","patchID","species","dbh","treeName")
#   missing_cols  <- setdiff(required_cols, names(init_trees))
#   if (length(missing_cols))
#     stop("init_trees is missing required columns: ", paste(missing_cols, collapse=", "))
#
#   # count expression: support either explicit 'trees' counts or one-row-per-tree
#   has_trees_col <- "trees" %in% names(init_trees)
#   count_expr <- if (has_trees_col) quote(sum(trees, na.rm = TRUE)) else quote(.N)
#
#   # helper: cohort a subset of trees with a given bin size
#   cohort_subset <- function(sub, binsize){
#     sub <- data.table::copy(sub)
#
#     if (binsize == 0) {
#       # keep actual DBH; one row per tree with trees = 1 (or existing trees if provided)
#       if (!has_trees_col) sub[, trees := 1L]
#       out <- sub[, .(trees = eval(count_expr)),
#                  by = .(siteID, patchID, species, dbh)]
#       return(out[])
#     }
#
#     # Positive bin size: cut DBH to mid-bin values and aggregate
#     if (is.null(min_dbh)) {
#       min_cut <- min(sub$dbh, na.rm = TRUE)
#     } else {
#       min_cut <- min_dbh
#     }
#
#     # Build breaks and labels (midpoints)
#     max_dbh <- max(sub$dbh, na.rm = TRUE)
#     brks    <- seq(min_cut, max_dbh + binsize, by = binsize)
#     mids    <- brks[-length(brks)] + binsize/2
#
#     sub[, dbh_binned := as.numeric(as.character(
#       cut(dbh, breaks = brks, labels = mids, include.lowest = TRUE)
#     ))]
#
#     out <- sub[, .(trees = eval(count_expr)),
#                by = .(siteID, patchID, species, dbh = dbh_binned)]
#     out[]
#   }
#
#   # Build per-group cohorts, then combine
#   cohorts <- lapply(seq_along(dbh_bins), function(i){
#     grp_names <- tree_groups[[i]]
#     binsize   <- dbh_bins[i]
#     sub       <- init_trees[treeName %in% grp_names]
#     if (nrow(sub) == 0L) return(data.table::data.table(
#       siteID=character(), patchID=character(), species=character(),
#       dbh=numeric(), trees=integer()
#     ))
#     cohort_subset(sub, binsize)
#   })
#
#   init_tbl <- data.table::rbindlist(cohorts, use.names = TRUE)
#
#   # deterministic cohort IDs
#   init_tbl[, cohortID := .I]
#
#   # FINN object
#   initCohort <- FINN::CohortMat$new(obs_df = init_tbl, sp = Nspecies)
#
#   if (return_tree_table) {
#     list(initCohort = initCohort, init_trees = init_tbl[])
#   } else {
#     initCohort
#   }
# }


add_site_completeness_flags <- function(env_dt = NULL,
                                        tree_dt = NULL,
                                        obs_dt = NULL) {
  stopifnot(!is.null(env_dt) || !is.null(tree_dt) || !is.null(obs_dt))

  # Collect all provided inputs
  inputs <- list(env_dt = env_dt, tree_dt = tree_dt, obs_dt = obs_dt)
  inputs <- inputs[!vapply(inputs, is.null, logical(1))]

  # Extract (siteName, year) keys for each
  key_list <- lapply(inputs, function(DT) unique(DT[, .(siteName, year)]))

  # Find intersection across provided tables
  common_keys <- Reduce(
    function(x, y) merge(x, y, by = c("siteName", "year"), all = FALSE),
    key_list
  )
  common_years <- common_keys[, .(common = list(sort(unique(year)))), by = siteName]
  setkey(common_years, siteName)

  # Helper: add flag to one DT by reference
  add_flag <- function(DT) {
    sy <- unique(DT[, .(siteName, year)])[
      , .(years = list(sort(unique(year)))), by = siteName
    ]
    setkey(sy, siteName)

    comp <- sy[common_years][
      , .(siteName, is_complete = unlist(Map(identical, years, common)))
    ]

    DT[, complete := FALSE]
    DT[comp, on = "siteName",
       complete := fifelse(is.na(i.is_complete), FALSE, i.is_complete)]
    invisible(NULL)
  }

  # Apply to each provided table
  if (!is.null(env_dt))  add_flag(env_dt)
  if (!is.null(tree_dt)) add_flag(tree_dt)
  if (!is.null(obs_dt))  add_flag(obs_dt)

  invisible(NULL)
}
