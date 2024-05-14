library(R6)

# Function to transform obs_dt into three arrays
obsDF2arrays = function(obs_dt) {
  # Retrieve dimensions from obs_dt
  Nsites <- length(unique(obs_dt$siteID))
  Npatches <- length(unique(obs_dt$patchID))
  maxNcohorts <- max(obs_dt$cohortID)

  # Initialize arrays
  species_array <- array(NA, dim = c(Nsites, Npatches, maxNcohorts))
  dbh_array <- array(NA, dim = c(Nsites, Npatches, maxNcohorts))
  nTree_array <- array(NA, dim = c(Nsites, Npatches, maxNcohorts))

  # Populate arrays
  for (i in 1:nrow(obs_dt)) {
    site <- obs_dt$siteID[i]
    patch <- obs_dt$patchID[i]
    cohort <- obs_dt$cohortID[i]

    species_array[site, patch, cohort] <- obs_dt$species[i]
    dbh_array[site, patch, cohort] <- obs_dt$dbh[i]
    nTree_array[site, patch, cohort] <- obs_dt$nTree[i]
  }

  # Return a list of the arrays
  return(list(
    species = species_array,
    dbh = dbh_array,
    nTree = nTree_array
  ))
}

library(data.table)

# Function to transform arrays back to a data.table
array2obsDF <- function(obs_array) {
  # Retrieve dimensions

  species_array = obs_array$species
  dbh_array = obs_array$dbh
  nTree_array = obs_array$nTree

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
        nTree <- nTree_array[site, patch, cohort]

        # Only add rows where dbh is not NA (assuming NA means no data for that cohort)
        if (!is.na(dbh)) {
          obs_dt <- rbind(obs_dt, data.table(
            siteID = site,
            patchID = patch,
            cohortID = cohort,
            species = species,
            dbh = dbh,
            nTree = nTree
          ))
        }
      }
    }
  }

  return(obs_dt)
}


CohortMat = R6::R6Class("CohortMat", public = list(
  dbh=NULL,
  nTree=NULL,
  Species=NULL,
  dims=c(50, 30, 10),
  sp = 10,
  device ='cpu',
  initialize = function(dbh=self$dbh, nTree=self$nTree, Species=self$Species, dims=self$dims, sp = self$sp, device = self$device, obs_df = NULL) {
    if(!is.null(obs_df)){
      if(all(colnames(obs_df) %in% c("siteID","patchID","species", "dbh", "nTree"))) stop("all(colnames(obs_df) %in% c(\"siteID\",\"patchID\",\"species\",\"dbh\",\"nTree\")) is not TRUE")
      obs_array = obsDF2arrays(obs_df)
      dbh = obs_array$dbh
      nTree = obs_array$nTree
      Species = obs_array$species
      dims = dim(obs_array$species)
    }

    self$dbh = if(is.null(dbh)) array(0.0, dim = dims) else dbh
    self$nTree = if(is.null(nTree)) array(0.0, dim = dims) else nTree
    self$Species = if(is.null(Species)) array(sample.int(sp, prod(dims),replace = TRUE), dim = dims) else Species
    self$dims = dims
    self$sp = sp
    self$device = torch::torch_device(device)

    if(!inherits(self$dbh, "torch_tensor")) self$dbh = torch::torch_tensor(self$dbh, dtype=torch_float32(), device=self$device)
    if(!inherits(self$nTree, "torch_tensor")) self$nTree = torch::torch_tensor(self$nTree, dtype=torch_float32(), device=self$device)
    if(!inherits(self$Species, "torch_tensor")) self$Species = torch::torch_tensor(self$Species, dtype=torch_int64(), device=self$device)

  },
  # Function to transform obs_dt into three arrays
  obsDF2arrays = obsDF2arrays
))

