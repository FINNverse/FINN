#' @export
# Define the climateDF2array function with automatic detection and error handling
climateDF2array <- function(climate_dt, env_vars) {
  # Ensure data.table is used
  if (!data.table::is.data.table(climate_dt)) {
    climate_dt <- data.table::as.data.table(climate_dt)
  }

  # Determine time resolution based on the presence of columns
  if ("day" %in% colnames(climate_dt)) {
    if (!("month" %in% colnames(climate_dt))) {
      stop("The 'day' column is provided, but the 'month' column is missing. Please provide 'month' along with 'day'.")
    }
    if (!("year" %in% colnames(climate_dt))) {
      stop("The 'day' column is provided, but the 'year' column is missing. Please provide 'year' along with 'day'.")
    }
    time_resolution <- "day"
  } else if ("month" %in% colnames(climate_dt)) {
    if (!("year" %in% colnames(climate_dt))) {
      stop("The 'month' column is provided, but the 'year' column is missing. Please provide 'year' along with 'month'.")
    }
    time_resolution <- "month"
  } else if ("year" %in% colnames(climate_dt)) {
    time_resolution <- "year"
  } else {
    stop("Neither 'year', 'month', nor 'day' columns are provided. Please provide at least the 'year' column.")
  }

  # Set the correct order of required columns based on time_resolution
  required_cols <- c("siteID", "year")

  if (time_resolution %in% c("month", "day")) {
    required_cols <- c(required_cols, "month")
  }

  if (time_resolution == "day") {
    required_cols <- c(required_cols, "day")
  }

  # Add environmental variables at the end of the required columns vector
  required_cols <- c(required_cols, env_vars)

  # Check if the required columns are present
  missing_cols <- setdiff(required_cols, colnames(climate_dt))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing from the input data: ", paste(missing_cols, collapse = ", "))
  }

  # Reorder columns to match the expected order using with = FALSE
  climate_dt <- climate_dt[, required_cols, with = FALSE]

  # Ensure the data is sorted correctly based on time_resolution
  if (time_resolution == "day") {
    data.table::setorder(climate_dt, siteID, year, month, day)
  } else if (time_resolution == "month") {
    data.table::setorder(climate_dt, siteID, year, month)
  } else {
    data.table::setorder(climate_dt, siteID, year)
  }

  # Retrieve dimensions from climate_dt
  Nsites <- length(unique(climate_dt$siteID))
  Nyears <- length(unique(climate_dt$year))
  Nenv <- length(env_vars)

  # Get unique ids for dimensions
  site_ids <- sort(unique(climate_dt$siteID))
  year_ids <- sort(unique(climate_dt$year))

  # Set dimensions based on time_resolution
  if (time_resolution == "day") {
    Nmonths <- 12
    Ndays <- 31
  } else if (time_resolution == "month") {
    Nmonths <- 12
    Ndays <- 1
  } else {
    Nmonths <- 1
    Ndays <- 1
  }

  # Convert climate_dt to a matrix for Rcpp
  climate_matrix <- as.matrix(climate_dt)

  # Call the Rcpp function
  env_array <- climateDF2arrayCpp(climate_matrix, as.integer(site_ids), as.integer(year_ids),
                                  time_resolution == "month", time_resolution == "day",
                                  Nsites, Nyears, Nmonths, Ndays, Nenv, env_vars)

  # Reshape the resulting vector into the appropriate array dimensions
  if (time_resolution == "day") {
    env_array <- array(env_array, dim = c(Nsites, Nyears, Nmonths, Ndays, Nenv))
  } else if (time_resolution == "month") {
    env_array <- array(env_array, dim = c(Nsites, Nyears, Nmonths, Nenv))
  } else {
    env_array <- array(env_array, dim = c(Nsites, Nyears, Nenv))
  }

  return(env_array)
}
