library(R6)

#' @useDynLib FINN, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Convert observation data frame to arrays
#'
#' @param obs_dt data.frame The observation data table containing siteID, patchID, cohortID, species, dbh, and trees columns.
#' @param additional_cols character vector Optional. Additional columns to be included as arrays.
#' @return A list of arrays for species, dbh, trees, and additional columns.
#' @examples
#' obs_dt <- data.frame(siteID = c(1, 1, 2), patchID = c(1, 2, 1), cohortID = c(1, 1, 2), species = c("A", "B", "A"), dbh = c(10, 20, 30), trees = c(100, 200, 150), height = c(5, 10, 15))
#' result <- obsDF2arrays(obs_dt, additional_cols = c("height"))
#' @export
obsDF2arrays <- function(obs_dt, additional_cols = character(0)) {
  if(!identical(sort(as.integer(obs_dt$cohortID)), sort(as.integer(as.factor(obs_dt$cohortID))))) {
    obs_dt$cohortID = as.integer(as.factor(obs_dt$cohortID))
    warning("cohortID must be a complete sequence of integers starting from 1. Cohort IDs were reassigned.")
  }
  if(!identical(sort(as.integer(obs_dt$siteID)), sort(as.integer(as.factor(obs_dt$siteID))))) {
    obs_dt$siteID = as.integer(as.factor(obs_dt$siteID))
    warning("siteID must be a complete sequence of integers starting from 1. Site IDs were reassigned.")
  }
  if(!identical(sort(as.integer(obs_dt$patchID)), sort(as.integer(as.factor(obs_dt$patchID))))) {
    obs_dt$patchID = as.integer(as.factor(obs_dt$patchID))
    warning("patchID must be a complete sequence of integers starting from 1. Patch IDs were reassigned.")
  }
  result <- obsDF2arraysCpp(obs_dt, additional_cols)

  # Reshape vectors to arrays
  dim <- result$dim
  species_array <- array(result$species, dim = dim)
  dbh_array <- array(result$dbh, dim = dim)
  trees_array <- array(result$trees, dim = dim)

  # Handle additional arrays
  additional_arrays <- lapply(additional_cols, function(col) {
    array(result[[col]], dim = dim)
  })
  names(additional_arrays) <- additional_cols

  return(c(list(
    species = species_array,
    dbh = dbh_array,
    trees = trees_array
  ), additional_arrays))
}

#' Generate Cohorts Using Weibull Distribution
#'
#' This function generates cohort data for tree populations using the Weibull distribution based on specified tree counts, diameter at breast height (DBH) shape, and scale parameters.
#'
#' @param trees Integer vector specifying the number of trees for each cohort. If a single integer is provided, it will be replicated for each draw.
#' @param dbh_shape Numeric vector specifying the shape parameters of the Weibull distribution for DBH for each cohort. If a single numeric value is provided, it will be replicated for each draw.
#' @param dbh_scale Numeric vector specifying the scale parameters of the Weibull distribution for DBH for each cohort. If a single numeric value is provided, it will be replicated for each draw.
#' @param dbh_class_range Numeric value specifying the range of DBH classes. Default is 1.
#' @param siteID Integer vector specifying the site ID for each cohort. If a single integer is provided, it will be replicated for each draw. Default is 1.
#' @param patchID Integer vector specifying the patch ID for each cohort. If a single integer is provided, it will be replicated for each draw. Default is 1.
#' @param species Integer vector specifying the species ID for each cohort. If a single integer is provided, it will be replicated for each draw. Default is 1.
#' @return A data frame containing the cohort data with columns for site ID, patch ID, cohort ID, species, number of trees, and DBH.
#' @details The function generates cohort data by drawing samples from the Weibull distribution for each cohort based on the specified shape and scale parameters. The resulting DBH values are binned into classes, and the cohort data is generated accordingly.
#' @examples
#' \dontrun{
#' obs_df <- rweibull_cohorts(
#'   trees = c(300, 10),
#'   dbh_shape = c(3, 3),
#'   dbh_scale = c(1, 50),
#'   dbh_class_range = 0.1,
#'   siteID = c(1, 1),
#'   patchID = c(1, 1),
#'   species = c(3, 4)
#' )
#' CohortMat$new(obs_df = obs_df)
#' }
#' @importFrom stats rweibull
#' @export
rweibull_cohorts = function(
  trees = NULL, # integer
  dbh_shape = NULL, # numeric
  dbh_scale = NULL, # numeric
  dbh_class_range = 1, # numeric
  siteID = 1, # integer
  patchID = 1, # integer
  species = 1 # integer
  ){
  n_draws = max(unlist(lapply(list(trees, dbh_shape, dbh_scale), length)))
  if(length(trees) == 1) trees = rep(trees, n_draws)
  if(length(dbh_shape) == 1) dbh_shape = rep(dbh_shape, n_draws)
  if(length(dbh_scale) == 1) dbh_scale = rep(dbh_scale, n_draws)
  if(length(siteID) == 1) siteID = rep(siteID, n_draws)
  if(length(patchID) == 1) patchID = rep(patchID, n_draws)
  if(length(species) == 1) species = rep(species, n_draws)

  if(!all(unlist(lapply(list(trees, dbh_shape, dbh_scale), length)) == length(trees))) {
    stop("inputs must all be of equal length ")
    }
  i=1
  cohortDF <-data.frame()
  cohortIDs = 0
  for(i in seq_along(trees)){
    drawn_pop = rweibull(trees[i], dbh_shape[i], dbh_scale[i])
    class_ranges = seq(0, ceiling(max(c(dbh_class_range,drawn_pop))), dbh_class_range)
    dbh_classes <- cut(drawn_pop,breaks = class_ranges)
    trees_out = table(dbh_classes)
    cohortDF <-
      rbind(cohortDF,
            data.frame(
              siteID = siteID[i],
              patchID = patchID[i],
              species = species[i],
              trees = as.integer(trees_out[trees_out > 0]),
              dbh = class_ranges[-1][trees_out > 0]+dbh_class_range
              )
        )
  }
  # return(list(cohortTensor = CohortMat$new(obs_df = cohortDF), cohortDF = cohortDF))
  cohortDF$cohortID <- 1:nrow(cohortDF)
  cohortDF <- cohortDF[,c("siteID","patchID", "cohortID", "species", "dbh", "trees")]
  return(cohortDF)
}

#' Transform Arrays to Observation Data Table
#'
#' This function transforms arrays of species, dbh, and trees back into an observation data table.
#'
#' @param obs_array A list containing three arrays: species, dbh, and trees.
#'
#' @return A data.frame with columns siteID, patchID, cohortID, species, dbh, and trees.
#'
#' @examples
#' obs_array <- list(species = array(c("A", "B"), dim = c(2, 2, 2)), dbh = array(c(10, 20, 30, 40), dim = c(2, 2, 2)), trees = array(c(100, 200, 150, 250), dim = c(2, 2, 2)))
#' result <- array2obsDF(obs_array)
#'
#' @export
array2obsDF <- function(obs_array) {
  # Retrieve dimensions

  species_array = obs_array$species
  dbh_array = obs_array$dbh
  trees_array = obs_array$trees

  Nsites <- dim(species_array)[1]
  Npatches <- dim(species_array)[2]
  maxNcohorts <- dim(species_array)[3]

  # Initialize an empty data.frame
  obs_dt <- data.frame()

  # Populate data.frame
  for (site in 1:Nsites) {
    for (patch in 1:Npatches) {
      for (cohort in 1:maxNcohorts) {
        species <- species_array[site, patch, cohort]
        dbh <- dbh_array[site, patch, cohort]
        trees <- trees_array[site, patch, cohort]

        # Only add rows where dbh is not NA (assuming NA means no data for that cohort)
        if (!is.na(dbh)) {
          obs_dt <- rbind(
            obs_dt,
            data.frame(
              siteID = site,
              patchID = patch,
              cohortID = cohort,
              species = species,
              dbh = dbh,
              trees = trees
            )
          )
        }
      }
    }
  }

  return(obs_dt)
}


#' Cohort Matrix Class
#'
#' An R6 class for managing cohorts of trees in forest models. This class allows for the initialization, transformation, and manipulation of cohorts represented by arrays of dbh, trees, and species.
#'
#' @field dbh A tensor representing the diameter at breast height for each cohort.
#' @field trees A tensor representing the number of trees in each cohort.
#' @field species A tensor representing the species of each cohort.
#' @field dims A vector representing the dimensions of the arrays (sites, patches, cohorts).
#' @field sp An integer representing the number of species.
#' @field device A character string specifying the device to use ('cpu' or 'cuda').
#' @field dbh_r A numeric array representing the diameter at breast height in R array format.
#' @field trees_r A numeric array representing the number of trees in R array format.
#' @field species_r An integer array representing the species in R array format.
#' @field device_r A character string specifying the device in R format ('cpu' or 'cuda').
#'
#' @param obs_df A data frame containing columns "siteID", "patchID", "species", "dbh", and "trees". If provided, it will be used to initialize the tensors.
#' @param dbh A tensor or array representing the diameter at breast height. Defaults to `self$dbh`.
#' @param trees A tensor or array representing the number of trees. Defaults to `self$trees`.
#' @param species A tensor or array representing the species. Defaults to `self$species`.
#' @param dims A numeric vector representing the dimensions of the arrays (sites, patches, cohorts). Defaults to `self$dims`.
#' @param sp An integer representing the number of species. Defaults to `self$sp`.
#' @param device A character string specifying the device to use ('cpu' or 'cuda'). Defaults to `self$device`.
#'
#' @export
CohortMat = R6::R6Class("CohortMat", public = list(
  dbh = NULL,
  trees = NULL,
  species = NULL,
  dbh_r = NULL,
  trees_r = NULL,
  species_r = NULL,
  dims = c(50, 30, 10),
  sp = 10,
  device = 'cpu',
  device_r = "cpu",

  #' @description
  #' Initialize the CohortMat class
  #' @param obs_df A data frame containing columns "siteID", "patchID", "species", "dbh", and "trees". If provided, it will be used to initialize the tensors.
  #' @param dbh A tensor or array representing the diameter at breast height. Defaults to `self$dbh`.
  #' @param trees A tensor or array representing the number of trees. Defaults to `self$trees`.
  #' @param species A tensor or array representing the species. Defaults to `self$species`.
  #' @param dims A numeric vector representing the dimensions of the arrays (sites, patches, cohorts). Defaults to `self$dims`.
  #' @param sp An integer representing the number of species. Defaults to `self$sp`.
  #' @param device A character string specifying the device to use ('cpu' or 'cuda'). Defaults to `self$device`.
  initialize = function(obs_df = NULL, dbh = self$dbh, trees = self$trees, species = self$species, dims = self$dims, sp = self$sp, device = self$device) {
    if (!is.null(obs_df)) {
      if (!all(c("siteID", "patchID", "species", "dbh", "trees") %in% colnames(obs_df)))
        stop('c("siteID", "patchID", "species", "dbh", "trees") %in% all(colnames(obs_df)) is not TRUE')
      obs_array = self$obsDF2arrays(obs_df)
      dbh = obs_array$dbh
      trees = obs_array$trees
      species = obs_array$species
      dims = dim(obs_array$species)
    }

    self$dbh = if (is.null(dbh)) array(0.0, dim = dims) else dbh
    self$trees = if (is.null(trees)) array(0.0, dim = dims) else trees
    self$species = if (is.null(species)) array(sample.int(sp, prod(dims), replace = TRUE), dim = dims) else species
    self$dims = dims
    self$sp = sp
    self$device_r = device
    if (!("torch_device" %in% class(device))) self$device = torch::torch_device(device)

    if (!inherits(self$dbh, "torch_tensor")) self$dbh = torch::torch_tensor(self$dbh, dtype = torch::torch_float32(), device = self$device)
    if (!inherits(self$trees, "torch_tensor")) self$trees = torch::torch_tensor(self$trees, dtype = torch::torch_float32(), device = self$device)
    if (!inherits(self$species, "torch_tensor")) self$species = torch::torch_tensor(self$species, dtype = torch::torch_int64(), device = self$device)

    self$dbh_r = torch::as_array(self$dbh)
    self$trees_r = torch::as_array(self$trees)
    self$species_r = torch::as_array(self$species)
  },

  #' Check and reinitialize tensors if necessary
  check = function() {
    self$device = torch::torch_device(self$device_r)
    self$dbh = check_and_recreate(self$dbh, self$dbh_r, dtype = torch::torch_float32(), device = self$device_r)
    self$trees = check_and_recreate(self$trees, self$trees_r, dtype = torch::torch_float32(), device = self$device_r)
    self$species = check_and_recreate(self$species, self$species_r, dtype = torch::torch_int64(), device = self$device_r)
  },

  #' Transform an observation data frame into three arrays
  #' @param obs_df A data frame containing columns "siteID", "patchID", "species", "dbh", and "trees".
  obsDF2arrays = obsDF2arrays,

  #' Convert the tensors to a data frame format
  asDF = function() {
    array2obsDF(list(dbh = torch::as_array(self$dbh), trees = torch::as_array(self$trees), species = torch::as_array(self$species)))
  }
))
