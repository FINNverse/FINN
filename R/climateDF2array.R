#' Transform Climate Data Table to Arrays
#'
#' This function transforms a climate data table into an array with dimensions for sites, year, and Nenv.
#' Optionally, it can create an array with dimensions for sites, year, month, and Nenv.
#'
#' @param climate_dt data.frame The climate data table containing year, month, uniquePLOTid, and environmental variables (e.g., tmp, pre) columns.
#' @param include_month logical Whether to include the month dimension in the output array. Defaults to FALSE.
#'
#' @return An array with dimensions for sites, year, and Nenv (or optionally, sites, year, month, and Nenv).
#'
#' @examples
#' climate_dt <- data.frame(
#'   year = rep(1960:1961, each = 2),
#'   month = rep(1:2, times = 2),
#'   uniquePLOTid = 1:4,
#'   tmp = runif(4, -10, 30),
#'   pre = runif(4, 0, 200)
#' )
#' result <- climateDF2array(climate_dt, include_month = TRUE)
#'
#' @export
climateDF2array <- function(climate_dt, include_month = FALSE) {
  # Retrieve dimensions from climate_dt
  Nsites <- length(unique(climate_dt$uniquePLOTid))
  Nyears <- length(unique(climate_dt$year))
  Nenv <- ncol(climate_dt) - 3  # Subtracting columns for year, month, and uniquePLOTid

  # Get unique ids for dimensions
  site_ids <- unique(climate_dt$uniquePLOTid)
  year_ids <- unique(climate_dt$year)

  if (include_month) {
    Nmonths <- 12
    # Initialize array with dimensions (sites, year, month, Nenv)
    env_array <- array(NA, dim = c(Nsites, Nyears, Nmonths, Nenv))
  } else {
    # Initialize array with dimensions (sites, year, Nenv)
    env_array <- array(NA, dim = c(Nsites, Nyears, Nenv))
  }

  # Populate the array
  for (i in 1:nrow(climate_dt)) {
    site <- which(site_ids == climate_dt$uniquePLOTid[i])
    year <- which(year_ids == climate_dt$year[i])
    env_values <- as.numeric(climate_dt[i, -(1:3)])  # Excluding year, month, and uniquePLOTid columns

    if (include_month) {
      month <- climate_dt$month[i]
      env_array[site, year, month, ] <- env_values
    } else {
      env_array[site, year, ] <- env_values
    }
  }

  return(env_array)
}

