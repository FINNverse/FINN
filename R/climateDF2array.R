#' @useDynLib FINN, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Transform Climate Data Table to Arrays
#'
#' This function transforms a climate data table into an array with dimensions for sites, year, and Nenv.
#' Optionally, it can create an array with dimensions for sites, year, month, and Nenv.
#'
#' @param climate_dt data.frame The climate data table containing year, month (optional), siteID, and environmental variables (e.g., tmp, pre) columns.
#' @param include_month logical Whether to include the month dimension in the output array. Defaults to FALSE.
#' @param env_vars character vector A vector of column names representing environmental variables in the climate data table.
#'
#' @return An array with dimensions for sites, year, and Nenv (or optionally, sites, year, month, and Nenv).
#'
#' @examples
#' climate_dt <- data.frame(
#'   year = rep(1960:1961, each = 2),
#'   month = rep(1:2, times = 2),
#'   siteID = 1:4,
#'   tmp = runif(4, -10, 30),
#'   pre = runif(4, 0, 200)
#' )
#' result <- climateDF2array(climate_dt, include_month = TRUE, env_vars = c("tmp", "pre"))
#'
#' @export
climateDF2array <- function(climate_dt, include_month = FALSE, env_vars) {
  # Ensure correct order in climate_dt
  climate_dt <- data.frame(climate_dt)
  if (include_month) {
    selected_cols <- c("year", "month", "siteID",  env_vars)
    climate_dt <- climate_dt[, selected_cols]
    climate_dt <- climate_dt[order(climate_dt$siteID, climate_dt$year, climate_dt$month), ]
  } else {
    selected_cols <- c("year", "siteID", env_vars)
    climate_dt <- climate_dt[, selected_cols]
    climate_dt <- climate_dt[order(climate_dt$siteID, climate_dt$year), ]
  }

  # Retrieve dimensions from climate_dt
  Nsites <- length(unique(climate_dt$siteID))
  Nyears <- length(unique(climate_dt$year))
  Nenv <- length(env_vars)

  # Get unique ids for dimensions
  site_ids <- unique(climate_dt$siteID)
  year_ids <- unique(climate_dt$year)

  if (include_month) {
    Nmonths <- 12
  } else {
    Nmonths <- 1
  }

  # Convert climate_dt to a matrix for Rcpp
  climate_matrix <- as.matrix(climate_dt)

  # Call the Rcpp function
  env_array <- climateDF2arrayCpp(climate_matrix, as.integer(site_ids), as.integer(year_ids), include_month, Nsites, Nyears, Nmonths, Nenv, env_vars)

  # Reshape the resulting vector into the appropriate array dimensions
  if (include_month) {
    env_array <- array(env_array, dim = c(Nsites, Nyears, Nmonths, Nenv))
  } else {
    env_array <- array(env_array, dim = c(Nsites, Nyears, Nenv))
  }

  return(env_array)
}

# Example usage
# climate_dt <- data.frame(
#   year = rep(1960:1961, each = 2),
#   month = rep(1:2, times = 2),
#   siteID = 1:4,
#   tmp = runif(4, -10, 30),
#   pre = runif(4, 0, 200)
# )
# result <- climateDF2array(climate_dt, include_month = TRUE, env_vars = c("tmp", "pre"))
# print(result)
