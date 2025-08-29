
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
resolveSiteIDs <- function(site_dt, tree_dt, env_dt = NULL, obs_dt = NULL,initCohort_dt = NULL, years2periods = F,
                           selection_priority = "fixed patches",
                           Npatches_fixed = 4,
                           Nyears_fixed = NULL, species2IDs = T,
                           all = F, createInitCohorts = F){

  only_sites <- unique(site_dt[,.(siteID2 = siteID, year)])
  only_patches <- unique(tree_dt[,.(siteID2 = siteID, patchID2 = patchID, year)])
  siteID_dt <- merge(only_patches, only_sites, by.x = c("siteID2", "year"), by.y = c("siteID2", "year"), all = all)

  siteID_dt[, NpatchesPerSite := uniqueN(patchID2), by = siteID2]
  siteID_dt[, NpatchesPerYear := uniqueN(patchID2), by = .(siteID2, year)]
  siteID_dt[, Nyears := uniqueN(year), by = siteID2]
  siteID_dt[, NyearsPerPatch := uniqueN(year), by = .(siteID2, patchID2)]
  siteID_dt[, NperPatch := .N, by = .(siteID2, patchID2)]
  siteID_dt[, NperYear := .N, by = year]
  siteID_dt[, maxYears := max(Nyears), by = .(siteID2, patchID2)]

  if(selection_priority == "fixed patches"){
    siteID_dt <- siteID_dt[NpatchesPerYear == Npatches_fixed]
  }else if(selection_priority == "many years"){
    siteID_dt <- siteID_dt[NyearsPerPatch == max(NyearsPerPatch)]
  }
  if(!is.null(Nyears_fixed)){
    siteID_dt <- siteID_dt[Nyears >= Nyears_fixed]
  }

  if(!is.null(env_dt)) {
    siteID_dt <- merge(siteID_dt, unique(env_dt[,.(siteID2 = siteID, year)]), by = c("siteID2", "year"), all = all)
  }

  if(!is.null(obs_dt)) {
    siteID_dt <- merge(siteID_dt, unique(obs_dt[,.(siteID2 = siteID, patchID2 = patchID, year)]), by = c("siteID2", "patchID2", "year"), all = all)
  }

  if(!is.null(initCohort_dt)){
    siteID_dt <- merge(siteID_dt, unique(initCohort_dt[,.(siteID2 = siteID, patchID2 = patchID)]), by = c("siteID2", "patchID2"), all = all)
  }

  siteID_dt[, siteID := as.integer(as.factor(siteID2)),]
  siteID_dt[, patchID := as.integer(as.factor(patchID2)), by = siteID]
  siteID_dt[, period := as.integer(as.factor(year)), by = .(siteID, patchID)]

  names(tree_dt)[names(tree_dt) == "siteID"] <- "siteID2"
  names(tree_dt)[names(tree_dt) == "patchID"] <- "patchID2"
  tree_dt <- merge(tree_dt, siteID_dt, by = c("siteID2", "patchID2", "year"), all = all)
  if(species2IDs){
    tree_dt[,species := as.integer(as.factor(species_name)),]
    species_dt <- unique(tree_dt[,.(species,species_name)])
  }

  if(years2periods) {
    tree_dt[, year2 := year,]
    tree_dt[, year := period,]
    tree_dt <- tree_dt[, -"period"]
  }

  names(site_dt)[names(site_dt) == "siteID"] <- "siteID2"
  site_dt <- merge(site_dt, siteID_dt[,.(siteID2, siteID, year, period)], by = c("siteID2", "year"), all = all)
  if(years2periods) {
    site_dt[, year2 := year,]
    site_dt[, year := period,]
    site_dt <- site_dt[, -"period"]
  }

  if(!is.null(env_dt)){
    names(env_dt)[names(env_dt) == "siteID"] <- "siteID2"
    env_dt <- merge(env_dt, siteID_dt[,.(siteID2, siteID, year, period)], by = c("siteID2","year"), all = all)
    if(years2periods) {
      env_dt[, year2 := year,]
      env_dt[, year := period,]
      env_dt <- env_dt[, -"period"]
    }
  }

  if(!is.null(obs_dt)){
    names(obs_dt)[names(obs_dt) == "siteID"] <- "siteID2"
    names(obs_dt)[names(obs_dt) == "patchID"] <- "patchID2"
    obs_dt <- merge(obs_dt, siteID_dt, by = c("siteID2", "patchID2", "year"), all = all)
    if(years2periods) {
      obs_dt[, year2 := year,]
      obs_dt[, year := period,]
      obs_dt <- obs_dt[, -"period"]
    }
    if(species2IDs){
      obs_dt <- merge(obs_dt, species_dt, by = "species_name", all = T)
      out_cols_obs_dt <- c("siteID", "patchID", "year", "ba", "dbh", "trees", "growth", "mort", "reg", "species", "species_name")
    }else{
      out_cols_obs_dt <- c("siteID", "patchID", "year", "ba", "dbh", "trees", "growth", "mort", "reg", "species_name")
    }
  }

  if(!is.null(initCohort_dt)){
    names(initCohort_dt)[names(initCohort_dt) == "siteID"] <- "siteID2"
    names(initCohort_dt)[names(initCohort_dt) == "patchID"] <- "patchID2"
    initCohort_dt <- merge(initCohort_dt, siteID_dt, by = c("siteID2", "patchID2"), all = all)
  }

  if(years2periods) {
    siteID_dt[, year2 := year, by = .(siteID, patchID)]
    siteID_dt[, year := period, by = .(siteID, patchID)]
    out_cols_siteID_dt <- c("siteID", "patchID", "siteID2", "patchID2", "year", "year2")
  } else {
    out_cols_siteID_dt <- c("siteID", "patchID", "siteID2", "patchID2", "year", "period")
  }

  out <- list(
    siteID_dt = as.data.table(siteID_dt[,..out_cols_siteID_dt], key = NULL),
    tree_dt = as.data.table(tree_dt, key = NULL),
    site_dt = as.data.table(site_dt, key = NULL),
    env_dt = as.data.table(env_dt, key = NULL),
    obs_dt = as.data.table(obs_dt[,..out_cols_obs_dt], key = NULL),
    initCohort_dt = as.data.table(initCohort_dt, key = NULL)
  )

  if(createInitCohorts){
    init_trees <- tree_dt[living == T & year == 1 & dbh > 0]
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

makeObsData <- function(tree_dt, plotsize){
  tree_dt <- copy(as.data.table(tree_dt))
  tree_dt[, dbh_before := data.table::shift(dbh,1,type = "lag"), by = treeID]
  tree_dt[, dbh_growth := dbh - dbh_before,]
  tree_dt[status == "alive", rel_growth := dbh/dbh_before - 1,]
  tree_dt[rel_growth < 0, rel_growth := 0,]
  tree_dt[living == T,trees := 1,]
  tree_dt[, cohortID := as.integer(as.factor(treeID)), by = .(siteID,patchID)]

  obs_dt <- tree_dt[,.(
    ba = sum(dbh2ba(dbh)*(living == T)),
    growth = mean(rel_growth[living == T], na.rm = T),
    dbh = mean(dbh[living == T], na.rm = T),
    trees = sum(living == T, na.rm = T),
    n_mort = sum(mort == T, na.rm = T),
    reg = sum(status == "new", na.rm = T)/plotsize
  ),by= .(siteID,patchID, year, species_name)]


  obs_dt[, trees_before := data.table::shift(trees, 1, type = "lag"), by = .(siteID, patchID, species_name)]
  obs_dt[, mort := 1-((trees_before-n_mort)/trees_before)^(1/1),]

  return(obs_dt[,.(siteID, patchID, year, species_name, ba, dbh, trees, growth, mort, reg)])
}


makeInitCohorts <- function(init_trees, dbh_binsize = NULL, min_dbh = NULL, Nspecies, treeID_table = F){
  if(!is.null(dbh_binsize)){
    if(is.null(min_dbh)) min_dbh = min(init_trees$dbh, na.rm = T)
    dbh_intervals = seq(min_dbh, max(init_trees$dbh, na.rm = T) + dbh_binsize, by = dbh_binsize)
    init_trees <- init_trees[,.(
      trees = sum(trees)),
      by = .(siteID, patchID, species, treeID,
             dbh = as.numeric(as.character(
               cut(dbh, breaks = dbh_intervals, labels = dbh_intervals[-1]-dbh_binsize/2, include.lowest = TRUE)
               ))
             )
      ]
  }

  init_trees$cohortID = 1:nrow(init_trees)

  initCohort <- FINN::CohortMat$new(obs_df = init_trees, sp = Nspecies)

  if(treeID_table) out = list(initCohort = initCohort, init_trees = init_trees) else if(!treeID_table) out = initCohort
  return(out)
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
