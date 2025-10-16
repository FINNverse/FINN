
#' Resolve siteIDs into array dimensions needed for the FINN simulation
#'
#' This function processes raw site and tree data to create a mapping of siteIDs and patchIDs into integer indices suitable for array dimensions in the FINN model.
#'
#' @param site_raw A data.table containing raw site-level data with columns "siteID" and "year".
#' @param tree_raw A data.table containing raw tree-level data with columns "siteID", "patchID", and "year".
#' @param selection_priority A string indicating the priority for selecting sites. Options are "fixed patches" or "many years". Default is "fixed patches".
#' @param Npatches_fixed An integer specifying the number of patches to fix when selection_priority is "fixed patches". Default is 4.
#' @param all A logical indicating whether to include all sites and patches, even if they do not match. Default is FALSE.
#' @return A data.table with resolved siteID and patchID indices along with their original IDs and years.
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
  out_cols_obs_dt <- c("siteID", "patchID", "year", "ba", "dbh", "trees", "growth", "mort", "reg", "species", "species_name")

  obs_dt_aggr <- obs_dt[,.(ba = sum(ba, na.rm = T),
                        dbh = mean(dbh, na.rm = T),
                        trees = sum(trees, na.rm = T),
                        growth = mean(growth, na.rm = T),
                        mort = mean(mort, na.rm = T),
                        reg = mean(reg, na.rm = T)),
                    by = .(siteID, patchID, year, species, species_name)]
  out_cols_obs_dt_aggr <- c("siteID", "year", "ba", "dbh", "trees", "growth", "mort", "reg", "species", "species_name")

  siteID_dt[, OrigYear := year, by = .(siteID, patchID)]
  siteID_dt[, year := period, by = .(siteID, patchID)]
  out_cols_siteID_dt <- c("siteID", "patchID", "siteName", "patchName", "year", "OrigYear")

  out <- list(
    siteID_dt = as.data.table(siteID_dt[,..out_cols_siteID_dt], key = NULL),
    tree_dt = as.data.table(tree_dt, key = NULL),
    env_dt = as.data.table(env_dt, key = NULL),
    obs_dt = unique(as.data.table(obs_dt_aggr[,..out_cols_obs_dt_aggr], key = NULL)),
    obs_dt_patches = unique(as.data.table(obs_dt[,..out_cols_obs_dt], key = NULL)),
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


makeEnvData <- function(){

  return(env_dt)
}

dbh2ba <- function(dbh){
  return((dbh/200)^2 * pi)
}

makeObsData <- function(tree_dt, plotsize, aggregate_by_site, minNyears = 2, Npatches = NULL, Nspecies = NULL, NspeciesQuantile = NULL){
  # browser()
  tree_dt <- copy(as.data.table(tree_dt))
  tree_dt[, dbh_before := data.table::shift(dbh,1,type = "lag"), by = treeName]
  tree_dt[, dbh_growth := dbh - dbh_before,]
  tree_dt[status == "alive", rel_growth := dbh/dbh_before - 1,]
  tree_dt[rel_growth < 0, rel_growth := 0,]
  tree_dt[living == T,trees := 1,]

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
    n_mort = sum(mort == T, na.rm = T),
    reg = sum(status == "new", na.rm = T)/plotsize
  ),by= .(siteName,patchName, year, species_name)]

  obs_dt[, trees_before := data.table::shift(trees, 1, type = "lag"), by = .(siteName, patchName, species_name)]
  obs_dt[, mort := 1-((trees_before-n_mort)/trees_before)^(1/1),]

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
  obs_dt[is.na(mort), mort := NA_real_]
  obs_dt[is.na(dbh), dbh := NA_real_]
  obs_dt = obs_dt[order(siteName, patchName, species_name, year)]

  obs_dt <- obs_dt[,.(
    ba = mean(ba, na.rm = T),
    growth = mean(growth, na.rm = T),
    dbh = mean(dbh, na.rm = T),
    trees = mean(trees, na.rm = T),
    mort = mean(mort, na.rm = T),
    reg = mean(reg, na.rm = T)
  ),by= .(siteName, year, species_name)]

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
      mort = mean(mort, na.rm = T),
      reg = mean(reg, na.rm = T),
      Npatches = uniqueN(patchName)
    ), by = .(siteName, year, species_name)]
    message("Aggregated obs_dt by siteName. For unaggregated data, set aggregate_by_site = FALSE")
  }else{
    obs_dt <- obs_dt[,.(siteName, patchName, year, species_name, ba, dbh, trees, growth, mort, reg)]
    message("Kept obs_dt unaggregated by siteName. For aggregated data, set aggregate_by_site = TRUE")
  }

  return(list(obs_dt = obs_dt, tree_dt = tree_dt))
}


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

makeInitCohorts2 <- function(
    init_trees,
    dbh_bins,                   # numeric vector (same length as tree_groups)
    tree_groups,                # list of character vectors (tree names)
    Nspecies,
    min_dbh = NULL,             # optional: scalar lower bound for all binning
    return_tree_table = FALSE   # formerly treeID_table
){
  stopifnot(is.list(tree_groups), length(dbh_bins) == length(tree_groups))

  # ensure data.table
  init_trees <- data.table::as.data.table(init_trees)

  required_cols <- c("siteID","patchID","species","dbh","treeName")
  missing_cols  <- setdiff(required_cols, names(init_trees))
  if (length(missing_cols))
    stop("init_trees is missing required columns: ", paste(missing_cols, collapse=", "))

  # count expression: support either explicit 'trees' counts or one-row-per-tree
  has_trees_col <- "trees" %in% names(init_trees)
  count_expr <- if (has_trees_col) quote(sum(trees, na.rm = TRUE)) else quote(.N)

  # helper: cohort a subset of trees with a given bin size
  cohort_subset <- function(sub, binsize){
    sub <- data.table::copy(sub)

    if (binsize == 0) {
      # keep actual DBH; one row per tree with trees = 1 (or existing trees if provided)
      if (!has_trees_col) sub[, trees := 1L]
      out <- sub[, .(trees = eval(count_expr)),
                 by = .(siteID, patchID, species, dbh)]
      return(out[])
    }

    # Positive bin size: cut DBH to mid-bin values and aggregate
    if (is.null(min_dbh)) {
      min_cut <- min(sub$dbh, na.rm = TRUE)
    } else {
      min_cut <- min_dbh
    }

    # Build breaks and labels (midpoints)
    max_dbh <- max(sub$dbh, na.rm = TRUE)
    brks    <- seq(min_cut, max_dbh + binsize, by = binsize)
    mids    <- brks[-length(brks)] + binsize/2

    sub[, dbh_binned := as.numeric(as.character(
      cut(dbh, breaks = brks, labels = mids, include.lowest = TRUE)
    ))]

    out <- sub[, .(trees = eval(count_expr)),
               by = .(siteID, patchID, species, dbh = dbh_binned)]
    out[]
  }

  # Build per-group cohorts, then combine
  cohorts <- lapply(seq_along(dbh_bins), function(i){
    grp_names <- tree_groups[[i]]
    binsize   <- dbh_bins[i]
    sub       <- init_trees[treeName %in% grp_names]
    if (nrow(sub) == 0L) return(data.table::data.table(
      siteID=character(), patchID=character(), species=character(),
      dbh=numeric(), trees=integer()
    ))
    cohort_subset(sub, binsize)
  })

  init_tbl <- data.table::rbindlist(cohorts, use.names = TRUE)

  # deterministic cohort IDs
  init_tbl[, cohortID := .I]

  # FINN object
  initCohort <- FINN::CohortMat$new(obs_df = init_tbl, sp = Nspecies)

  if (return_tree_table) {
    list(initCohort = initCohort, init_trees = init_tbl[])
  } else {
    initCohort
  }
}


#' Create inputs for the forest model
#' This function prepares input data for FINN by processing site and tree data. It resolves siteIDs into array dimensions needed for the FINN simulation.
#' @param site_dt A data.table containing site-level data.
#' @param tree_dt A data.table containing tree-level data.
#' @param Nspecies An integer specifying the maximum number of species to consider.
#' @param patchsize A numeric value specifying the patch size for the simulations. Default is 0.06.
#' @param dbh_binsize A numeric value specifying the bin size for diameter at breast height (DBH) when creating initial cohorts. Default is 0.5.
#' @return A list containing the following elements:
#' \item{env_dt}{A data.table with environmental data.}
#' \item{obs_dt}{A data.table with observation data.}
#' \item{initCohort}{An object representing the initial cohorts of trees.}
#' \item{patchsize}{The patch size used for the simulations.}
#' \item{Nspecies}{The maximum number of species considered.}
#' @import data.table
#' @export
createInputs <- function(site_dt, tree_dt, Nspecies, patchsize = patchsize, dbh_binsize = 0.1){

  siteID_dt <- resolveSiteIDs(env_dt, obs_dt, initCohort)

  tree_dt <- merge(tree_dt, siteID_dt, by = "siteID", all.x = TRUE)
  site_dt <- merge(site_dt, siteID_dt, by = "siteID", all.x = TRUE)
  env_dt <- merge(env_dt, siteID_dt, by = "siteID", all.x = TRUE)
  obs_dt <- merge(obs_dt, siteID_dt, by = "siteID", all.x = TRUE)

  env_dt <- makeEnvData(site_dt)

  obs_dt <- makeObsData(tree_dt, Nspecies = Nspecies)

  initCohort <- makeInitCohorts(tree_dt, dbh_binsize, Nspecies = Nspecies, site_dt)

  if(checkInputs(env_dt, obs_dt, initCohort, patchsize) == FALSE) stop("Input check failed")

  return(list(env_dt = env_dt, obs_dt = obs_dt, initCohort = initCohort, patchsize = patchsize, Nspecies = Nspecies))
}



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
