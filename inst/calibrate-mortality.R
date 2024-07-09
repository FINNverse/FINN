#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Project: FINN
## Script purpose:  Explore FIA data for calibration of mortality parameters.
## Date: Mon Jul  8 09:17:04 2024
## Author: Yannek Kaeber <y.kaeber@posteo.de>
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

library(data.table)

# Load data
plot_dt <- fread("data/calibration-data/FIA/prepared/plot_dt.csv")
state_dt <- fread("data/calibration-data/FIA/prepared/state_dt.csv")
treevars_dt <- fread("data/calibration-data/FIA/prepared/treevars_dt.csv")
climate_dt <- fread("data/calibration-data/FIA/prepared/plot_cruts_climate_dt.csv")

plot_dt <- merge(plot_dt, state_dt[,.(STATECD = FIPS_Code, State_Name)], by = "STATECD")
treevars_dt <- merge(treevars_dt, unique(plot_dt[,.(uniquePLOTid,State_Name)]), by = "uniquePLOTid")

# select the states Oregon (41), Florida (12), California (06), Colorado (08), Maine (23)
selected_plot_dt <- plot_dt[STATECD %in% c(41, 12, 06, 08, 23)]
selected_treevars_dt <- treevars_dt[uniquePLOTid %in% selected_plot_dt$uniquePLOTid]
selected_climate_dt <- climate_dt[uniquePLOTid %in% selected_plot_dt$uniquePLOTid]

hist(selected_climate_dt$tmp)
hist(selected_climate_dt$pre)
hist(treevars_dt$DBH)
hist(treevars_dt$ACTUALHTm)

heightPars_dt <- calibrate_height(treevars_dt)

cohort_dt <-
  treevars_dt[,.(
    siteID = uniquePLOTid,
    patchID = 1,
    cohortID = uniqueTREEid,
    species = species,
    dbh = DBH,
    trees = TPA_UNADJ
    )]

climate_array <- climateDF2array(selected_climate_dt)
dim(climate_array)
climate_array2 <- climateDF2array(selected_climate_dt, include_month = T)
dim(climate_array2)
cohort_array <- FINN::obsDF2arrays(cohort_dt)
dim(cohort_array)
