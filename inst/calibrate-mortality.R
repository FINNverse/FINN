#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Project: FINN
## Script purpose:  Explore FIA data for calibration of mortality parameters.
## Date: Mon Jul  8 09:17:04 2024
## Author: Yannek Kaeber <y.kaeber@posteo.de>
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# Conversion functions
circle_area_cm2m2 <- function(d) pi * (d/100 / 2)^2
inches2cm <- function(x) x * 2.54
ft2m <- function(x) x * 0.3048
sqrft2m2 <- function(x) x * 0.092903
acre2ha <- function(x) x * 0.404686
sqrft2acre <- function(x) x * 2.29568e-5
radius2area <- function(x) pi * x^2

# Load required libraries
source("R/calibrate-height.R")
library(data.table)
library(FINN)

# Define constants
plotRadius = 24
sampled_area_acre <- sqrft2acre(radius2area(plotRadius)) * 4
sampled_area_ha <- acre2ha(sampled_area_acre)
representated_area_acre <- 1
representated_area_ha <- acre2ha(representated_area_acre)

# Load data
plot_dt <- fread("data/calibration-data/FIA/preparedV3/plot.csv")
state_dt <- fread("data/calibration-data/FIA/preparedV3/state.csv")
mort_calib_data <- fread("data/calibration-data/FIA/preparedV3/mort_calib_data.csv")
cohort_mort_calib_data_dt <- fread("data/calibration-data/FIA/preparedV3/cohort_mort_calib_data_dt.csv")
climate_dt <- fread("data/calibration-data/FIA/preparedV3/plot_cruts_climate_dt.csv")
species_dt <- fread("data/calibration-data/FIA/FIADB_REFERENCE/REF_SPECIES.csv")

# Merge plot and state data
plot_dt <- merge(plot_dt, state_dt[, .(STATECD = FIPS_Code, State_Name)], by = "STATECD")
mort_calib_data <- merge(mort_calib_data, unique(plot_dt[, .(uniquePLOTid, State_Name)]), by = "uniquePLOTid")

# Select specific states
selected_plot_dt <- plot_dt[STATECD %in% c(41, 12, 06, 08, 23)]

# Filter plots based on criteria
selected_plot_dt <- selected_plot_dt[
  PLOT_STATUS_CD == 1 &
    QA_STATUS == 1 &
    INTENSITY == 1 &
    DESIGNCD == 1 &
    INVYR != 9999
]

# Select relevant climate and mortality data
selected_climate_dt <- climate_dt[uniquePLOTid %in% selected_plot_dt$uniquePLOTid]
mort_calib_data <- mort_calib_data[uniquePLOTid %in% selected_plot_dt$uniquePLOTid]

# Calculate death years and update mortality data
mort_calib_data[, yearDead := min(year[y == 1], na.rm = TRUE), by = uniqueTREEid]

# Filter data based on death year
mort_calib_data2 <- mort_calib_data[year <= yearDead]
mort_calib_data2[, time_of_death_years := yearDead - 1959]
mort_calib_data2[, time_of_death_months := time_of_death_years * 12]

# Calculate patch area in hectares
acre_to_0.067haPatch <- (1 / representated_area_ha) * sampled_area_ha

# Create cohort data table
cohort_dt <- mort_calib_data2[DBH >= 12.7, .(
  uniquePLOTid,
  siteID = as.integer(as.factor(uniquePLOTid)),
  patchID = 1,
  uniqueTREEid,
  species = species,
  dbh = DBH,
  y,
  trees = round(trees_acre * acre_to_0.067haPatch),
  time_of_death_years,
  time_of_death_months,
  year
)]
cohort_dt[, cohortID := as.integer(factor(uniqueTREEid, ordered = TRUE)), by = .(siteID)]

# Convert cohort data to arrays
cohort_array <- obsDF2arrays(cohort_dt, additional_cols = c("time_of_death_months", "y"))
dim(cohort_array$trees)
dim(cohort_array$species)
dim(cohort_array$dbh)
dim(cohort_array$time_of_death_months)
dim(cohort_array$y)

# Calibrate height parameters
heightPars_dt <- calibrate_height(mort_calib_data)

# Filter climate data for selected plots and years
selected_climate_dt <- selected_climate_dt[
  uniquePLOTid %in% unique(cohort_dt$uniquePLOTid) &
    year <= max(cohort_dt$year)
]

# Add running month column
selected_climate_dt[, year2 := as.integer(factor(year, ordered = TRUE))]
selected_climate_dt[, running_month := (12 * year2 - 12) + month]

# Check max running month
max(selected_climate_dt$running_month)
max(cohort_dt$time_of_death_months)

# Convert climate data to array
climate_array <- climateDF2array(
  climate_dt = selected_climate_dt[, .(
    uniquePLOTid = as.integer(factor(uniquePLOTid)),
    year = running_month,
    pre = round(pre, 2),
    tmp = round(tmp),
    2
  )],
  include_month = F,
  env_vars = c("pre", "tmp")
)


# This is the climate array
dim(climate_array)

# Those are the cohort arrays
dim(cohort_array$species)
dim(cohort_array$dbh)
dim(cohort_array$trees)
dim(cohort_array$time_of_death_months)
dim(cohort_array$y)

#Have fun
