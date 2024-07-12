#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Project: FINN
## Script purpose:  Explore FIA data for calibration of mortality parameters.
## Date: Mon Jul  8 09:17:04 2024
## Author: Yannek Kaeber <y.kaeber@posteo.de>
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
library(data.table)
library(FINN)

plotRadius = 24
sampled_area_acre <- sqrft2acre(radius2area(plotRadius))*4
sampled_area_ha <- acre2ha(sqrft2acre(radius2area(plotRadius))*4)
representated_area_acre <- 1
representated_area_ha <- acre2ha(1)

# Load data
plot_dt <- fread("data/calibration-data/FIA/preparedV3/plot.csv")
state_dt <- fread("data/calibration-data/FIA/preparedV3/state.csv")
mort_calib_data <- fread("data/calibration-data/FIA/preparedV3/mort_calib_data.csv")
cohort_mort_calib_data_dt <- fread("data/calibration-data/FIA/preparedV3/cohort_mort_calib_data_dt.csv")
climate_dt <- fread("data/calibration-data/FIA/preparedV3/plot_cruts_climate_dt.csv")
plotID_dt <- fread("data/calibration-data/FIA/plotID_dt.csv")

plot_dt <- merge(plot_dt, state_dt[,.(STATECD = FIPS_Code, State_Name)], by = "STATECD")
mort_calib_data <- merge(mort_calib_data, unique(plot_dt[,.(uniquePLOTid,State_Name)]), by = "uniquePLOTid")

# select the states Oregon (41), Florida (12), California (06), Colorado (08), Maine (23)
selected_plot_dt <- plot_dt[STATECD %in% c(41, 12, 06, 08, 23)]
# select plots based on PLOT_STATUS_CD, QA_STATUS, INTENSITY, DESIGNCD, INVYR
# PLOTSTATUS_CD = 1 --> Only sampled plots are included
# QA_STATUS = 1 --> Only plots that have passed the quality assurance
# INTENSITY = 1 --> Only plots that are part of the FIA standard grid
# DESIGNCD = 1 --> Only plots that follow the "harmonized" set of the FIA sampling
# INVYR != 9999 --> Only plots that have a valid inventory year
selected_plot_dt <- selected_plot_dt[PLOT_STATUS_CD == 1 & QA_STATUS == 1 & INTENSITY == 1 & DESIGNCD == 1 & INVYR != 9999]
# selected_plot_dt[,uniquePLOTOBSid := paste0(uniquePLOTid,INVYR),]
# treevars_dt[,uniquePLOTOBSid := paste0(uniquePLOTid,year),]
# selected_treevars_dt <- treevars_dt[uniquePLOTOBSid %in% selected_plot_dt$uniquePLOTOBSid]
selected_climate_dt <- climate_dt[uniquePLOTid %in% selected_plot_dt$uniquePLOTid]
mort_calib_data <- mort_calib_data[uniquePLOTid %in% selected_plot_dt$uniquePLOTid]


uniqueN(selected_plot_dt$uniquePLOTid)
uniqueN(selected_treevars_dt$uniquePLOTid)
uniqueN(selected_climate_dt$uniquePLOTid)

heightPars_dt <- calibrate_height(treevars_dt)

cohort_dt <-
  treevars_dt[,.(
    # siteID = uniquePLOTid,
    siteID = as.integer(as.factor(uniquePLOTid)),
    patchID = 1,
    cohortID = uniqueTREEid,
    # cohortID = as.integer(as.factor(uniqueTREEid)),
    species = species,
    dbh = DBH,
    fresh_dead,
    trees = round(TPA_UNADJ),
    year
    )]
uniqueN(treevars_dt$uniquePLOTid)
cohort_dt <- cohort_dt[siteID < 10000]
cohort_dt[,":="(
         periodID = as.integer(factor(year, ordered = T)),
         Nperiods = .N
         ),
         by = .(siteID,
                species = as.integer(as.factor(species)),
                dbh,
                trees)
         ]

cohort_dt[,":="(
         cohortID = as.integer(as.factor(cohortID)),
         patchID = sample.int(10, .N, replace = T)
         ),
         by = .(siteID,
                species = as.integer(as.factor(species)),
                dbh,
                trees)
         ]

# cohort_dt <- cohort_dt[Nperiods > 1]
table(cohort_dt$periodID, cohort_dt$year)
cohort_array <- obsDF2arrays(cohort_dt[periodID == 1])
dim(cohort_array$trees)
dim(cohort_array$dbh)
dim(cohort_array$species)

library(FINN)
climate_array <- climateDF2array(climate_dt = selected_climate_dt[year < 1962], include_month = T)
dim(climate_array)
dim(climate_array)
climate_array2 <- climateDF2array(selected_climate_dt, include_month = T)
dim(climate_array2)
