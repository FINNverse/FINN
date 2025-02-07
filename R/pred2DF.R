#' Convert Prediction Arrays to Data Frames
#'
#' This function takes prediction arrays from a model output and converts them
#' into data frames. The data frames can be returned in either 'wide' or 'long'
#' format. The function processes site-level, patch-level, and cohort-level predictions.
#'
#' @param pred A list containing prediction arrays for site-level, patch-level,
#' and cohort-level data. Each element of `pred` should have a structure similar
#' to what is outlined in the details section.
#' @param format A character string indicating the desired format of the output
#' data frames. Must be either "wide" (default) or "long".
#'
#' @details
#' The `pred` argument should be a list containing at least a `Predictions` element,
#' which itself is a list of arrays. The arrays represent predictions for different
#' metrics such as dbh/ba, tree counts, AL (aboveground live biomass), growth rates,
#' mortality rates, and regeneration rates for sites, patches, or cohorts. The
#' dimensionality of the arrays should correspond to different factors, such as
#' siteID, year, species, and optionally patch or cohortID.
#'
#' The function first converts each prediction array into a data frame, properly
#' naming and converting the relevant dimensions. It then merges these data frames
#' by common identifiers such as siteID, year, species, patch, and cohortID.
#' Depending on the `format` parameter, the data frames are returned in either a
#' wide format (one row per site/patch/cohort per year with multiple columns for
#' different metrics) or a long format (one row per site/patch/cohort per year
#' per metric).
#'
#' @return A list of data frames. The list may contain up to three elements:
#' `site`, `patch`, and `cohort`, corresponding to the processed site-level,
#' patch-level, and cohort-level predictions, respectively.
#'
#' @examples
#' \dontrun{
#' # Assuming `model_output` is a list with the structure required by `pred2DF`
#' result <- pred2DF(model_output, format = "long")
#' print(result$site)
#' }
#'
#' @export
pred2DF <- function(pred, format = "wide") {
  # Initialize an empty data.table to store the site-level results
  site_dt = data.table()

  # Loop through the site-level prediction metrics
  for (i in c("dbh", "ba", "trees", "growth", "mort", "reg", "r_mean_ha")) {
    # Set dimension names
    dim_names = c("siteID", "year", "species")
    # Convert the prediction array to a data frame
    site_temp_dt <- data.table(
      as.data.frame.table(
        pred$Predictions$Site[[i]],
        responseName = i)
    )
    # Rename dimension columns
    setnames(site_temp_dt, old = paste0("Var", 1:length(dim_names)), new = dim_names)
    # Convert dimension columns to integers
    for(dim_i in dim_names) site_temp_dt[[dim_i]] <- as.integer(site_temp_dt[[dim_i]])

    # Merge data frames if not the first iteration
    if(nrow(site_dt) == 0) {
      site_dt <- site_temp_dt
    } else {
      site_dt <- merge(site_dt, site_temp_dt, by = c("siteID", "species", "year"))
    }
  }

  # Process patch-level predictions if available
  if("Patch" %in% names(pred$Predictions)) {
    patch_dt <- data.table()
    for(i in names(pred$Predictions$Patch[[1]])) {
      patch_temp_dt2 <- data.table()
      for(year_i in seq_along(pred$Predictions$Patch)) {
        # Set dimension names
        dim_names = c("siteID", "patchID", "species")
        # Convert the patch-level array to a data frame
        patch_temp_dt <- data.table(
          as.data.frame.table(
            pred$Predictions$Patch[[year_i]][[i]],
            responseName = i)
        )
        # Rename dimension columns
        setnames(patch_temp_dt, old = paste0("Var", 1:length(dim_names)), new = dim_names)
        # Convert dimension columns to integers
        for(dim_i in dim_names) patch_temp_dt[[dim_i]] <- as.integer(patch_temp_dt[[dim_i]])
        patch_temp_dt$year = year_i
        patch_temp_dt2 <- rbind(patch_temp_dt2, patch_temp_dt)
        setcolorder(patch_temp_dt, c("siteID", "patchID", "species", "year", i))
      }
      if(nrow(patch_dt) == 0) {
        patch_dt <- patch_temp_dt2
      } else {
        patch_dt <- merge(patch_dt, patch_temp_dt, by = c("siteID", "patchID", "species", "year"))
      }
    }
  }

  # Process cohort-level predictions if available
  if("Cohort" %in% names(pred$Predictions)) {
    cohort_dt <- data.table()
    for(year_i in seq_along(pred$Predictions$Cohort$cohortID)) {
      # Select the cohort array for the current year
      cohorts_array_i = pred$Predictions$Cohort$cohortID[[year_i]]
      # Set dimension names
      dim_names = c("siteID", "patchID", "species")
      # Convert the cohort-level array to a data frame
      cohort_temp_dt <- data.table(
        as.data.frame.table(
          cohorts_array_i,
          responseName = "cohortID")
      )
      # Rename dimension columns
      setnames(cohort_temp_dt, old = paste0("Var", 1:length(dim_names)), new = dim_names)
      # Convert dimension columns to integers
      for(dim_i in dim_names) cohort_temp_dt[[dim_i]] <- as.integer(cohort_temp_dt[[dim_i]])

      # Remove rows where cohortID equals 0
      #cohort_temp_dt <- cohort_temp_dt[cohortID != 0]
      cohort_temp_dt <- cohort_temp_dt[cohort_temp_dt$cohortID != 0,]

      # Get the indices of the cohortID in the cohort array
      idx_mat <- matrix(nrow = 0, ncol = 3)
      for(cohortID_i in cohort_temp_dt$cohortID) {
        idx_mat = rbind(idx_mat, which(cohortID_i == cohorts_array_i, arr.ind = T))
      }

      # Add the year column to the data table
      cohort_temp_dt$year = year_i

      # Add the cohort data to cohort_dt
      for (i in names(pred$Predictions$Cohort$cohortStates[[year_i]])) {
        cohort_temp_dt[[i]] = pred$Predictions$Cohort$cohortStates[[year_i]][[i]][idx_mat]
      }
      cohort_dt <- rbind(cohort_dt, cohort_temp_dt)
    }
  }

  # Convert to long format if specified
  if(format == "long") {
    if(exists("site_dt")) site_dt <- melt(site_dt, id.vars = c("siteID", "year", "species"), variable.name = "variable")
    if(exists("patch_dt")) patch_dt <- melt(patch_dt, id.vars = c("siteID", "patchID", "year", "species"), variable.name = "variable")
    if(exists("cohort_dt")) cohort_dt <- melt(cohort_dt, id.vars = c("siteID", "patchID", "year", "species", "cohortID"), variable.name = "variable")
  } else if(format != "wide") {
    stop("Invalid format argument. Use either 'long' or 'wide'")
  }

  # Prepare the output list
  out <- list()
  if(exists("site_dt")) out$site <- site_dt
  if(exists("patch_dt")) out$patch <- patch_dt
  if(exists("cohort_dt")) out$cohort <- cohort_dt

  return(out)
}
