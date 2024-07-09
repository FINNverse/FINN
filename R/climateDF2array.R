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
library(Rcpp)
sourceCpp("inst/rcpp_functions.cpp")

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
  } else {
    Nmonths <- 1
  }

  # Start time for ETA calculation
  start_time <- Sys.time()
  total_rows <- nrow(climate_dt)

  # Convert climate_dt to a matrix for Rcpp
  climate_matrix <- as.matrix(climate_dt)

  # Call the Rcpp function
  env_array <- populate_array(climate_matrix, site_ids, year_ids, include_month, Nsites, Nyears, Nmonths, Nenv)

  # Reshape the resulting vector into the appropriate array dimensions
  if (include_month) {
    env_array <- array(env_array, dim = c(Nsites, Nyears, Nmonths, Nenv))
  } else {
    env_array <- array(env_array, dim = c(Nsites, Nyears, Nenv))
  }

  # Calculate elapsed time and estimate remaining time
  elapsed_time <- Sys.time() - start_time
  avg_time_per_iteration <- elapsed_time / total_rows
  remaining_time <- avg_time_per_iteration * (total_rows - total_rows)
  percentage_done <- (total_rows / total_rows) * 100

  # Overwrite previous console output with estimated remaining time and percentage
  cat(sprintf("\r%d of %d - %.2f%% - Elapsed time: %.2f secs, Estimated remaining time: %.2f secs",
              total_rows, total_rows, percentage_done, as.numeric(elapsed_time, units = "secs"), as.numeric(remaining_time, units = "secs")))

  cat("\n")
  return(env_array)
}

# # Example usage
# climate_dt <- data.frame(
#   year = rep(1960:1961, each = 2),
#   month = rep(1:2, times = 2),
#   uniquePLOTid = 1:4,
#   tmp = runif(4, -10, 30),
#   pre = runif(4, 0, 200)
# )
# result <- climateDF2array(climate_dt, include_month = TRUE)
# print(result)
