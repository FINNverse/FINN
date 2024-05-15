library(R6)
library(data.table)

# Function to transform obs_dt into three arrays
obsDF2arrays = function(obs_dt) {
  # Retrieve dimensions from obs_dt
  Nsites <- length(unique(obs_dt$siteID))
  Npatches <- length(unique(obs_dt$patchID))
  maxNcohorts <- max(obs_dt$cohortID)

  # Initialize arrays
  species_array <- array(NA, dim = c(Nsites, Npatches, maxNcohorts))
  dbh_array <- array(NA, dim = c(Nsites, Npatches, maxNcohorts))
  trees_array <- array(NA, dim = c(Nsites, Npatches, maxNcohorts))

  # Populate arrays
  for (i in 1:nrow(obs_dt)) {
    site <- obs_dt$siteID[i]
    patch <- obs_dt$patchID[i]
    cohort <- obs_dt$cohortID[i]

    species_array[site, patch, cohort] <- obs_dt$species[i]
    dbh_array[site, patch, cohort] <- obs_dt$dbh[i]
    trees_array[site, patch, cohort] <- obs_dt$trees[i]
  }

  # Return a list of the arrays
  return(list(
    species = species_array,
    dbh = dbh_array,
    trees = trees_array
  ))
}

# Function to transform arrays back to a data.table
array2obsDF <- function(obs_array) {
  # Retrieve dimensions

  species_array = obs_array$species
  dbh_array = obs_array$dbh
  trees_array = obs_array$trees

  Nsites <- dim(species_array)[1]
  Npatches <- dim(species_array)[2]
  maxNcohorts <- dim(species_array)[3]

  # Initialize an empty data.table
  obs_dt <- data.table()

  # Populate data.table
  for (site in 1:Nsites) {
    for (patch in 1:Npatches) {
      for (cohort in 1:maxNcohorts) {
        species <- species_array[site, patch, cohort]
        dbh <- dbh_array[site, patch, cohort]
        trees <- trees_array[site, patch, cohort]

        # Only add rows where dbh is not NA (assuming NA means no data for that cohort)
        if (!is.na(dbh)) {
          obs_dt <- rbind(obs_dt, data.table(
            siteID = site,
            patchID = patch,
            cohortID = cohort,
            species = species,
            dbh = dbh,
            trees = trees
          ))
        }
      }
    }
  }

  return(obs_dt)
}

CohortMat = R6::R6Class("CohortMat", public = list(
  dbh=NULL,
  trees=NULL,
  species=NULL,
  dims=c(50, 30, 10),
  sp = 10,
  device ='cpu',
  initialize = function(dbh=self$dbh, trees=self$trees, species=self$species, dims=self$dims, sp = self$sp, device = self$device, obs_df = NULL) {
    if(!is.null(obs_df)){
      if(all(colnames(obs_df) %in% c("siteID","patchID","species", "dbh", "trees"))) stop("all(colnames(obs_df) %in% c(\"siteID\",\"patchID\",\"species\",\"dbh\",\"trees\")) is not TRUE")
      obs_array = obsDF2arrays(obs_df)
      dbh = obs_array$dbh
      trees = obs_array$trees
      species = obs_array$species
      dims = dim(obs_array$species)
    }

    self$dbh = if(is.null(dbh)) array(0.0, dim = dims) else dbh
    self$trees = if(is.null(trees)) array(0.0, dim = dims) else trees
    self$species = if(is.null(species)) array(sample.int(sp, prod(dims),replace = TRUE), dim = dims) else species
    self$dims = dims
    self$sp = sp
    self$device = torch::torch_device(device)

    if(!inherits(self$dbh, "torch_tensor")) self$dbh = torch::torch_tensor(self$dbh, dtype=torch_float32(), device=self$device)
    if(!inherits(self$trees, "torch_tensor")) self$trees = torch::torch_tensor(self$trees, dtype=torch_float32(), device=self$device)
    if(!inherits(self$species, "torch_tensor")) self$species = torch::torch_tensor(self$species, dtype=torch_int64(), device=self$device)

  },
  # Function to transform obs_dt into three arrays
  obsDF2arrays = obsDF2arrays,
  asDF = function() {
    array2obsDF(list(dbh = torch::as_array(self$dbh),trees = torch::as_array(self$trees),species = torch::as_array(self$species)))
    }
))

