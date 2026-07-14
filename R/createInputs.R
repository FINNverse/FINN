
#' Resolve site, patch, and year indices for FINN inputs
#'
#' Joins tree, environment, and observation tables; assigns integer indices
#' for site, patch, period; standardizes species coding; and returns aligned
#' data plus optional initial cohorts.
#'
#' @param tree_dt data.table with tree-level data including siteName, patchName, year,
#'   species_name, dbh, status, living, and optional trees.
#' @param env_dt data.table with environment data including siteName and year.
#' @param obs_dt data.table with observations including siteName, patchName, year,
#'   species_name, and stand metrics (ba, dbh, trees, growth, mort, n_at_risk,
#'   n_died, reg), as returned by \code{makeObsData}.
#' @param createInitCohorts logical. If TRUE, build FINN initial cohorts from
#'   trees with \code{living == TRUE} and \code{year == 1}. Default TRUE.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{siteID_dt}: site/patch/year index map.
#'   \item \code{tree_dt}: tree table with indices and standardized species.
#'   \item \code{env_dt}: environment table with indices.
#'   \item \code{obs_dt}: site-aggregated observations by species.
#'   \item \code{obs_dt_patches}: patch-level observations by species.
#'   \item \code{species_dt}: species lookup (\code{species}, \code{species_name}).
#'   \item \code{initCohorts} (optional): FINN \code{CohortMat} object.
#'   \item \code{initCohort_dt} (optional): trees used to build cohorts.
#' }
#'
#' @import data.table
#' @export
resolveSiteIDs <- function(tree_dt, env_dt, obs_dt, createInitCohorts = T){
  # browser()
  siteID_dt <- unique(env_dt[,.(siteName, year)])
  siteID_dt <- merge(siteID_dt, unique(tree_dt[,.(siteName, year)]), by = c("siteName", "year"), all = F)
  siteID_dt <- merge(siteID_dt, unique(obs_dt[,.(siteName, patchName, year)]), by = c("siteName", "year"), all = F)

  siteID_dt[, siteID := as.integer(as.factor(siteName)),]
  siteID_dt[, patchID := as.integer(as.factor(patchName)), by = siteID]
  siteID_dt[, period := as.integer(as.factor(year)), by = .(siteID, patchID)]

  tree_dt <- merge(tree_dt, siteID_dt, by = c("siteName", "patchName", "year"), all = F)

  species_levels <- sort(unique(tree_dt$species_name))
  species_levels <- c(species_levels[species_levels != "other"], "other")

  tree_dt[,species := as.integer(factor(species_name, levels = species_levels, labels = 1:length(species_levels))),]
  species_dt <- unique(tree_dt[,.(species,species_name)])

  tree_dt[, OrigYear := year,]
  tree_dt[, year := period-1,]
  tree_dt <- tree_dt[, -"period"]

  env_dt <- merge(env_dt, siteID_dt[,.(siteName, siteID, year, period)], by = c("siteName","year"), all = F)
  env_dt[, OrigYear := year,]
  env_dt[, year := period-1,]
  env_dt <- env_dt[, -"period"]

  obs_dt <- merge(obs_dt, siteID_dt, by = c("siteName", "patchName", "year"), all = F)
  obs_dt[, OrigYear := year,]
  obs_dt[, year := period-1,]
  obs_dt <- obs_dt[, -"period"]
  obs_dt <- merge(obs_dt, species_dt, by = "species_name", all = T)
  out_cols_obs_dt <- c("siteID", "patchID", "year", "ba", "dbh", "trees", "growth", "mort", "n_at_risk", "n_died", "reg", "species", "species_name")

  obs_dt_aggr <- obs_dt[,.(ba = mean(ba, na.rm = T),
                        dbh = mean(dbh, na.rm = T),
                        trees = mean(trees, na.rm = T),
                        growth = mean(growth, na.rm = T),
                        # mortality counts pool by summing; the rate is derived
                        # from the pooled counts, never averaged from rates.
                        n_at_risk = sum(n_at_risk, na.rm = T),
                        n_died = sum(n_died, na.rm = T),
                        reg = mean(reg, na.rm = T)),
                    by = .(siteID, year, species, species_name)]
  obs_dt_aggr[, mort := data.table::fifelse(n_at_risk > 0, n_died / n_at_risk, NA_real_)]
  out_cols_obs_dt_aggr <- c("siteID", "year", "ba", "dbh", "trees", "growth", "mort", "n_at_risk", "n_died", "reg", "species", "species_name")

  siteID_dt[, OrigYear := year, by = .(siteID, patchID)]
  siteID_dt[, year := period, by = .(siteID, patchID)]
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
    init_trees <- tree_dt[living == T & year == 1]
    if(!any(colnames(init_trees) %in% "trees")) init_trees[, trees := 1,]
    out[["initCohorts"]] <- makeInitCohorts(init_trees, Nspecies = uniqueN(obs_dt$species), treeID_table = F)
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
#'
#' @return A list with:
#' \itemize{
#'   \item \code{obs_dt}: observations at site or patch level. Mortality comes
#'     back as a closed-cohort pair of counts, \code{n_at_risk} (trees alive at
#'     the start of the interval) and \code{n_died} (how many of them were dead
#'     at the end), plus the derived rate \code{mort = n_died / n_at_risk}
#'     (\code{NA} where no cohort was at risk). The counts are the binomial
#'     response; pass them to \code{fit} with \code{mortality = "binomial"}.
#'   \item \code{tree_dt}: input trees with added growth fields and species recode.
#' }
#'
#' @import data.table
#' @export
makeObsData <- function(tree_dt, plotsize, aggregate_by_site = T, minNyears = 2, fix_period_length = NULL, dbh_growth_thresh = c(-10,50), Npatches = NULL, Nspecies = NULL, NspeciesQuantile = NULL){
  # browser()
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
  # tree_dt[rel_growth < 0, rel_growth := 0,]
  tree_dt[living == T,trees := 1,]

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
    alive_species <- names(sort(table(td[status %chin% c("alive","new")]$species_name), decreasing = TRUE))
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
    ba = sum(dbh2ba(dbh)*(living == T), na.rm = T),
    growth = mean(rel_growth[living == T], na.rm = T),
    dbh = mean(dbh[living == T], na.rm = T),
    trees = sum(living == T, na.rm = T),
    reg = sum(status == "new", na.rm = T)/plotsize
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
  obs_dt[is.na(dbh), dbh := NA_real_]
  # No cohort at risk -> the rate is undefined, not zero. Keep the counts at 0 so
  # they still sum correctly on aggregation, and let `mort` carry the NA so the
  # loss masks the observation out.
  obs_dt[is.na(n_at_risk), n_at_risk := 0L]
  obs_dt[is.na(n_died),    n_died    := 0L]
  obs_dt[, mort := data.table::fifelse(n_at_risk > 0, n_died / n_at_risk, NA_real_)]
  obs_dt = obs_dt[order(siteName, patchName, species_name, year)]

  if(!is.null(minNyears)){
    obs_dt[, NyearsPerPatch := uniqueN(year), by = .(siteName, patchName)]
    obs_dt[, sameYearsAllPatches := all(NyearsPerPatch == max(NyearsPerPatch)), by = siteName]
    obs_dt <- obs_dt[sameYearsAllPatches == T & minNyears <= NyearsPerPatch]
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
      ba = mean(ba, na.rm = T),
      trees = mean(trees, na.rm = T),
      dbh = mean(dbh, na.rm = T),
      growth = mean(growth, na.rm = T),
      # counts SUM across patches; the site rate is then the pooled rate. The
      # old `mean(mort)` weighted a 1-tree patch the same as a 50-tree one.
      n_at_risk = sum(n_at_risk, na.rm = T),
      n_died = sum(n_died, na.rm = T),
      reg = mean(reg, na.rm = T),
      Npatches = uniqueN(patchName)
    ), by = .(siteName, year, species_name)]
    obs_dt[, mort := data.table::fifelse(n_at_risk > 0, n_died / n_at_risk, NA_real_)]
    message("Aggregated obs_dt by siteName. For unaggregated data, set aggregate_by_site = FALSE")
  }else{
    obs_dt <- obs_dt[,.(siteName, patchName, year, species_name, ba, dbh, trees, growth, mort, n_at_risk, n_died, reg)]
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
#' @import data.table
#' @export
makeInitCohorts <- function(init_trees, dbh_binsize = NULL, min_dbh = NULL, Nspecies, treeID_table = F, singleCohortTreeNames = NULL){
  # browser()
  if(!is.null(dbh_binsize)){
    if(is.null(min_dbh)) min_dbh = min(init_trees$dbh, na.rm = T)
    dbh_intervals = seq(min_dbh, max(init_trees$dbh, na.rm = T) + dbh_binsize, by = dbh_binsize)
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
