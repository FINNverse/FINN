library(data.table)
library(raster)

plot_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv")
plot_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT),]

plotid_dt <- fread("data/calibration-data/FIA/plotID_dt.csv")
plot_dt <- merge(plot_dt, plotid_dt, by= "uniquePLOTid_txt")

selected_plots <- plot_dt[
  PLOT_STATUS_CD == 1 & QA_STATUS == 1 & INTENSITY == 1 & DESIGNCD == 1 & INVYR != 9999
]

# Create a data.table to store the extracted values
plot_dt <- unique(selected_plots[, .(uniquePLOTid, LAT, LON)])

tiff_dir <- "data/calibration-data/CRUTS/tiffs"

# Define chunk size
chunk_size <- 2000000

# Specify the minimum and maximum year for extraction
min_year <- 1960  # Set the minimum year
max_year <- 2023  # Set the maximum year


variables <- c("tmp", "pre")

# Measure the time taken for the extraction process in chunks
total_time_taken <- system.time({
  for (variable in variables) {
    output_csv_path <- paste0("data/calibration-data/FIA/plot_cruts_",variable,".csv")
    # Get list of raster files for the current variable
    raster_files <- list.files(tiff_dir, pattern = paste0(variable, "_year_.*\\.tif$"), full.names = TRUE)

    for (raster_file_path in raster_files) {
      # Extract the year from the raster file name
      year <- as.numeric(gsub(paste0(variable, "_year_|\\.tif"), "", basename(raster_file_path)))

      # Check if the year is within the specified range
      if (year >= min_year && year <= max_year) {
        # Load the raster stack for the current year
        year_stack <- stack(raster_file_path)

        # Process in chunks
        for (i in seq(1, nrow(plot_dt), by = chunk_size)) {
          end <- min(i + chunk_size - 1, nrow(plot_dt))
          coords_chunk <- plot_dt[i:end, .(LON, LAT)]

          # Measure time for each chunk
          out_dt <- data.table()
          chunk_time_taken <- system.time({
            extracted_values <- extract(year_stack, coords_chunk)
            for (j in 1:12) {
              temp_plot_dt <- plot_dt[i:end, .(uniquePLOTid)]
              temp_plot_dt$year <- year
              temp_plot_dt$month <- j
              temp_plot_dt$variable <- variable
              temp_plot_dt$value <- extracted_values[, j]
              out_dt <- rbind(out_dt, temp_plot_dt)
            }
            if (year == min_year)
              fwrite(out_dt, output_csv_path) else
                fwrite(out_dt, output_csv_path, append = TRUE)
          })

          # Convert time to minutes and seconds
          chunk_time_minutes <- floor(chunk_time_taken[3] / 60)
          chunk_time_seconds <- round(chunk_time_taken[3] %% 60, 4)

          # Print the time taken for this chunk in minutes and seconds
          cat("Time taken for chunk", variable, "year", year, i, "to", end, ":", chunk_time_minutes, "minutes", chunk_time_seconds, "seconds\n")
        }
      }
    }
  }
})

# Print the total time taken
total_time_minutes <- floor(total_time_taken[3] / 60)
total_time_seconds <- round(total_time_taken[3] %% 60, 4)
cat("Total time taken:", total_time_minutes, "minutes", total_time_seconds, "seconds\n")

tmp <- fread("data/calibration-data/FIA/plot_cruts_tmp.csv")
pre <- fread("data/calibration-data/FIA/plot_cruts_pre.csv")

plot_cruts_climate_dt <- merge(
  tmp[,.(uniquePLOTid,year,month,tmp = value)],
  pre[,.(uniquePLOTid,year,month,pre = value)],
  by = c("year","month","uniquePLOTid"))

fwrite(plot_cruts_climate_dt, "data/calibration-data/FIA/preparedV3/plot_cruts_climate_dt.csv")
