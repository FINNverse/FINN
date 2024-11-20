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
# select the states based on the state code
# the selected states are Oregon, Florida, California, Colorado, and Maine
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
mort_calib_data2 <- mort_calib_data[year <= yearDead | y == 0]
table(mort_calib_data2$y)/sum(table(mort_calib_data2$y))


mort_calib_data2[, time_of_death_years := yearDead - 1959]
mort_calib_data2[, time_of_death_months := time_of_death_years * 12]
# mort_calib_data2[, DBH := DBH_before,]

# Calculate patch area in hectares
acre_to_0.067haPatch <- (1 / representated_area_ha) * sampled_area_ha

# Create cohort data table
cohort_dt <- copy(mort_calib_data2[!is.na(DBH_before) & (DBH_before >= 12.7), .(
  uniquePLOTid,
  # siteID = as.integer(as.factor(uniquePLOTid)),
  siteID = as.integer(as.factor(paste0(uniquePLOTid,year))),
  patchID = 1,
  uniqueTREEid,
  SPCD = species,
  species = as.integer(factor(species)),
  dbh = DBH_before,
  y,
  trees = 1,
  time_of_death_years,
  time_of_death_months,
  height = actualHeight,
  ba = BA_total,
  year
)])

table(cohort_dt$y)/sum(table(cohort_dt$y))

cohort_dt <- merge(cohort_dt, species_dt[,.(SPCD, SCIENTIFIC_NAME)], by = "SPCD", all.x = TRUE)
View(cohort_dt[, .(
  N = .N,
  "0" = sum(y == 0),
  "1" = sum(y == 1)
  ), by = .(SPCD, species, SCIENTIFIC_NAME)][order(-N)])

cohort_dt[, cohortID := 1:.N, by = .(siteID)]

# Convert cohort data to arrays
cohort_array <- obsDF2arrays(cohort_dt, additional_cols = c("time_of_death_months", "y", "height", "ba"))
sum(!is.na(cohort_array$trees),na.rm = T)
sum(!is.na(cohort_array$species),na.rm = T)
sum(!is.na(cohort_array$dbh),na.rm = T)
sum(!is.na(cohort_array$time_of_death_months),na.rm = T)
sum(!is.na(cohort_array$height),na.rm = T)
dim(cohort_array$species)
dim(cohort_array$height)
dim(cohort_array$dbh)
hist(cohort_array$dbh)
hist(cohort_array$dbh)
dim(cohort_array$time_of_death_months)
dim(cohort_array$y)

# Calibrate height parameters
heightPars_dt <- calibrate_height(mort_calib_data)

# Filter climate data for selected plots and years
selected_climate_dt <- selected_climate_dt[
  uniquePLOTid %in% unique(cohort_dt$uniquePLOTid) &
    year <= max(cohort_dt$year)
]

selected_climate_dt2 <- merge(selected_climate_dt, unique(cohort_dt[,.(uniquePLOTid, siteID)]), by = "uniquePLOTid", allow.cartesian=TRUE)

# Add running month column
selected_climate_dt2[, year2 := as.integer(factor(year, ordered = TRUE))]
selected_climate_dt2[, running_month := (12 * year2 - 12) + month]

# Check max running month
max(selected_climate_dt2$running_month)
max(cohort_dt$time_of_death_months)

# Convert climate data to array
climate_array <- climateDF2array(
  climate_dt = selected_climate_dt2[, .(
    siteID,
    year = running_month,
    pre = round(pre, 1),
    tmp = round(tmp, 1)
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
