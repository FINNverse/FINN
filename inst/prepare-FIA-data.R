#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Project: FINN
## Script purpose:  Explore FIA data for calibration of mortality parameters.
## Date: Mon Jul  8 09:17:04 2024
## Author: Yannek Kaeber <y.kaeber@posteo.de>
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

library(data.table)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## conversion tables ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
state_dt <- data.table(
  FIPS_Code = c(1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56),
  State_Name = c('Alabama', 'Alaska', 'Arizona', 'Arkansas', 'California', 'Colorado', 'Connecticut', 'Delaware', 'District of Columbia', 'Florida',
                 'Georgia', 'Hawaii', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Maryland', 'Massachusetts',
                 'Michigan', 'Minnesota', 'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico',
                 'New York', 'North Carolina', 'North Dakota', 'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island', 'South Carolina',
                 'South Dakota', 'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia', 'Washington', 'West Virginia', 'Wisconsin', 'Wyoming')
)
species_dt <- fread("data/calibration-data/FIA/FIADB_REFERENCE/REF_SPECIES.csv")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## select plots based on plot variables ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
plot_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv")
# subplot_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_SUBP_COND.csv")
# cond_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_COND.csv")
# pop_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_POP_ESTN_UNIT.csv")
# subpcond_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_SUBP_COND.csv")

# create uniquePLOTid based on "natural" unique cariables for PLOT TABLE
# see FIA documentation for details
plot_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT), by = INVYR]
plot_dt[, uniquePLOTid := as.integer(factor(uniquePLOTid_txt, levels = unique(uniquePLOTid_txt)))]

plotID_dt <- unique(plot_dt[, .(uniquePLOTid, uniquePLOTid_txt)])

# save table with unique IDs for linking climate variable extraction to the correct PLOT
fwrite(plotID_dt, "data/calibration-data/FIA/plotID_dt.csv")

# select plots based on PLOT_STATUS_CD, QA_STATUS, INTENSITY, DESIGNCD, INVYR
# PLOTSTATUS_CD = 1 --> Only sampled plots are included
# QA_STATUS = 1 --> Only plots that have passed the quality assurance
# INTENSITY = 1 --> Only plots that are part of the FIA standard grid
# DESIGNCD = 1 --> Only plots that follow the "harmonized" set of the FIA sampling
# INVYR != 9999 --> Only plots that have a valid inventory year
selected_plots <- plot_dt[
  PLOT_STATUS_CD == 1 & QA_STATUS == 1 & INTENSITY == 1 & DESIGNCD == 1 & INVYR != 9999
    ]

table(selected_plots$MACRO_BREAKPOINT_DIA, useNA = "always")
selected_plots[is.na(MACRO_BREAKPOINT_DIA)]

tree_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_TREE.csv")
tree_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT),]
tree_dt <- merge(tree_dt, plotID_dt, by = "uniquePLOTid_txt", all.x = TRUE)
tree_dt[, uniqueTREEid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT, "_", SUBP, "_", TREE),]
tree_dt[, uniqueTREEid := as.integer(factor(uniqueTREEid_txt, levels = unique(uniqueTREEid_txt)))]
tree_dt[, uniqueSUBPid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT, "_", SUBP),]
tree_dt[, uniqueSUBPid := as.integer(factor(uniqueSUBPid_txt, levels = unique(uniqueSUBPid_txt)))]
selected_trees_dt <- tree_dt[uniquePLOTid %in% unique(selected_plots$uniquePLOTid)]
selected_trees_dt <- selected_trees_dt[INVYR != 9999]

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## calculate variables ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

circle_area <- function(d) pi * (d/100 / 2)^2
inches2cm <- function(x) x * 2.54
ft2m <- function(x) x * 0.3048
sqrft2m2 <- function(x) x * 0.092903
acre2ha <- function(x) x * 0.404686
sqrft2acre <- function(x) x * 2.29568e-5
radius2area <- function(x) pi * x^2

plotRadius = 24
sampled_area_acre <- sqrft2acre(radius2area(plotRadius))*4
sampled_area_ha <- acre2ha(sqrft2acre(radius2area(plotRadius))*4)
representated_area_acre <- 1
representated_area_ha <- acre2ha(1)

selected_trees_dt[,":="(
  DBHcm = inches2cm(DIA),
  height = ft2m(HT)
),]
selected_trees_dt[,ba := round(circle_area(DBHcm), 4),]

plotvars_dt <- selected_trees_dt[,.(
  N = sum((STATUSCD == 1)*TPA_UNADJ, na.rm = T),
  BA = sum(ba*(STATUSCD == 1)*TPA_UNADJ, na.rm = T)
), by = .(uniquePLOTid, INVYR)]

selected_trees_dt <- merge(selected_trees_dt, plotvars_dt, by = c("uniquePLOTid", "INVYR"))

# selected_trees_dt <- selected_trees_dt[,.(
#   uniquePLOTid,
#   uniqueSUBPid,
#   uniqueTREEid,
#   year = INVYR,
#   species = SPCD,
#   DBH = DBHcm,
#   height = ft2m(HT),
#   DIACHECK = DIACHECK, # 0 = accurate, 1 = estimated, 2 = deltawrong, 5 = modeled
#   HTCD = HTCD, # 1 = measured, 2 = totvisest-actmeas, 3 = visest, 4 = estimated
#   ACTUALHTm = ft2m(ACTUALHT),
#   TPA_UNADJ,
#   status = STATUSCD # 0 = none, 1 = alive, 2 = dead, 3 = cut
# )]

selected_trees_dt <- selected_trees_dt[order(INVYR)]

selected_trees_dt[
  order(INVYR),
  ":="(
    STATUSCD_before = data.table::shift(STATUSCD,1,type = "lag"),
    DBHcm_before = data.table::shift(DBHcm,1,type = "lag")
  ),
  by = .(uniquePLOTid, uniqueSUBPid, uniqueTREEid)]

selected_trees_dt[, fresh_dead := STATUSCD == 2 & STATUSCD_before == 1,]
selected_trees_dt[DIACHECK == 0, DBHchange := DBHcm - DBHcm_before,]

selected_trees_dt[
  order(INVYR),
  ":="(
    DBHchange_before = data.table::shift(DBHchange,1,type = "lag")
  ),
  by = .(uniquePLOTid, uniqueSUBPid, uniqueTREEid)]

selected_trees_dt[fresh_dead == T, mort_year := INVYR,]
selected_trees_dt[, trees_acre := TPA_UNADJ,]

mort_calib_data <- selected_trees_dt[STATUSCD == 1 | STATUSCD_before == 1,.(
  uniquePLOTid,
  uniqueSUBPid,
  uniqueTREEid,
  year = INVYR,
  species = SPCD,
  trees_acre,
  DBH = DBHcm,
  DBH_before = DBHcm_before,
  DBHchange,
  DBHchange_before,
  height = ft2m(HT),
  BA_total = BA,
  N_total = N,
  y = as.integer(fresh_dead)
)]


maxDBH <- ceiling(max(selected_trees_dt$DBHcm,na.rm = T)/10)*10
intervalsize <- 10
DBHclasses <- seq(0, maxDBH, intervalsize)
DBHclassesLabels <- DBHclasses[-1]-intervalsize/2
# table for stand initialization
cohort_mort_calib_data_dt <-
  mort_calib_data[, .(trees_acre = sum(trees_acre),
                      N_dead = sum(trees_acre * (y == 1))),
                  by = .(
                    uniquePLOTid,
                    species,
                    year,
                    BA_total,
                    N_total,
                    DBHclass = as.numeric(as.character(cut(DBH, breaks = DBHclasses, labels = DBHclassesLabels)))
                    )]

outputDir <- "data/calibration-data/FIA/preparedV3"
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## save the full tables ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

fwrite(selected_trees_dt, paste0(outputDir, "selected_trees.csv"))
fwrite(mort_calib_data, paste0(outputDir, "mort_calib_data.csv"))
fwrite(cohort_mort_calib_data_dt, paste0(outputDir, "cohort_mort_calib_data_dt.csv"))
fwrite(state_dt, paste0(outputDir, "state.csv"))
fwrite(plot_dt, paste0(outputDir, "plot.csv"))
fwrite(species_dt, paste0(outputDir, "species.csv"))
