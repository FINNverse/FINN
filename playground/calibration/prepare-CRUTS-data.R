#
# # Load the raster and ncdf4 packages
# library(raster)
# library(ncdf4)
#
# # Function to extract time information from the NetCDF file
# get_nc_time_info <- function(nc_file_path) {
#   # Open the NetCDF file using ncdf4
#   nc <- nc_open(nc_file_path)
#
#   # Extract the time variable
#   time <- ncvar_get(nc, "time")
#
#   # Get the units attribute to understand the time origin
#   time_units <- ncatt_get(nc, "time", "units")$value
#
#   # Close the NetCDF file
#   nc_close(nc)
#
#   # Convert the time variable to POSIXct (assuming units are "days since YYYY-MM-DD")
#   time_origin <- sub("days since ", "", time_units)
#   time_dates <- as.POSIXct(time * 86400, origin = time_origin, tz = "UTC")
#
#   # Extract the year, month, and day for each time point
#   time_years <- as.numeric(format(time_dates, "%Y"))
#   time_months <- as.numeric(format(time_dates, "%m"))
#   time_days <- as.numeric(format(time_dates, "%d"))
#
#   # Combine the band, year, month, and day into a data frame
#   time_info <- data.frame(band = seq_along(time_dates), year = time_years, month = time_months, day = time_days)
#
#   return(time_info)
# }
#
# # Define the root path to the data directory
# root_path <- "data/calibration-data/CRUTS"
#
# # Define the subdirectories for different variables
# # variables <- c("tmp", "pre")
# # variables <- c("tmp")
# variables <- c("pre")
#
# # Create a directory to store the output raster files (if it doesn't already exist)
# output_dir <- file.path(root_path, "tiffs")
# if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#
# # Iterate over each variable
# for (variable in variables) {
#   # Define the path to the NetCDF files for the current variable
#   var_path <- file.path(root_path, variable)
#
#   # Get a list of NetCDF files for the current variable
#   nc_files <- list.files(var_path, pattern = "\\.nc$", full.names = TRUE)
#
#   # Iterate over each NetCDF file
#   for (file_path in nc_files) {
#     # Extract time information from the NetCDF file
#     time_df <- get_nc_time_info(file_path)
#
#     # Read the NetCDF file into a RasterBrick object
#     nc_data <- brick(file_path, varname = variable)
#
#     # Loop through each year
#     for (year in unique(time_df$year)) {
#       # Select the layers for the current year
#       year_layers <- which(time_df$year == year)
#
#       # Create a stack of 12 layers for the current year (one for each month)
#       year_stack <- stack()
#       for (month in 1:12) {
#         month_layer <- which(time_df$year == year & time_df$month == month)
#         if (length(month_layer) > 0) {
#           year_stack <- addLayer(year_stack, nc_data[[month_layer]])
#         } else {
#           # If no data for this month, add an empty layer
#           year_stack <- addLayer(year_stack, raster(nrows=nrow(nc_data), ncols=ncol(nc_data), vals=NA))
#         }
#       }
#
#       # Define the path for the output raster file
#       output_raster_path <- file.path(output_dir, paste0(variable, "_year_", year, ".tif"))
#
#       # Write the raster stack to a file
#       writeRaster(year_stack, output_raster_path, format = "GTiff", overwrite = TRUE)
#
#       # Print a message indicating the process is complete for the current year
#       cat("Raster file for", variable, "year", year, "has been created and saved at:", output_raster_path, "\n")
#     }
#   }
# }
