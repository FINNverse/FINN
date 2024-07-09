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

# create uniquePLOTid based on "natural" unique cariables for PLOT TABLE
# see FIA documentation for details
plot_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT),]
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

# create Unique IDs for tree table based on "natural" unique variables for TREE TABLE
# see FIA documentation for details
if(file.exists("data/calibration-data/FIA/ENTIRE_TREE_uniquePLOTid.csv")) {
  tree_dt <- fread("data/calibration-data/FIA/ENTIRE_TREE_uniquePLOTid.csv")
}else{
  tree_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_TREE.csv")
  tree_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT),]
  tree_dt <- merge(tree_dt, plotID_dt, by = "uniquePLOTid_txt", all.x = TRUE)
  tree_dt[, uniqueTREEid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT, "_", SUBP, "_", TREE),]
  tree_dt[, uniqueTREEid := as.integer(factor(uniqueTREEid_txt, levels = unique(uniqueTREEid_txt)))]
  tree_dt[, uniqueSUBPid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT, "_", SUBP),]
  tree_dt[, uniqueSUBPid := as.integer(factor(uniqueSUBPid_txt, levels = unique(uniqueSUBPid_txt)))]
  fwrite(tree_dt, "data/calibration-data/FIA/ENTIRE_TREE_uniquePLOTid.csv")
}


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## calculate variables ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

circle_area <- function(d) pi * (d/100 / 2)^2
inches2cm <- function(x) x * 2.54
ft2m <- function(x) x * 0.3048
sqrft2m2 <- function(x) x * 0.092903

tree_dt[,":="(
  DBHcm = inches2cm(DIA)
),]
tree_dt[,ba := round(circle_area(DBHcm), 4),]

species_subplotvars_dt <- tree_dt[,.(
  N = sum(STATUSCD == 1)*TPA_UNADJ,
  BA = sum(ba*(STATUSCD == 1))*TPA_UNADJ,
  meanDBH = weighted.mean(DBHcm, TPA_UNADJ, na.rm = TRUE)
), by = .(uniquePLOTid, uniqueSUBPid, SPCD, year = INVYR)]

# table for stand initialization
cohort_species_subplotvars_dt <- tree_dt[,.(
  N = sum(STATUSCD == 1)*TPA_UNADJ
), by = .(uniquePLOTid, uniqueSUBPid, SPCD, year = INVYR, DBHclass = cut(DBHcm, breaks = seq(0, 260, 5)))]

allspecies_subplotvars_dt <- tree_dt[,.(
  N = sum(STATUSCD == 1)*TPA_UNADJ,
  BA = sum(ba*(STATUSCD == 1))*TPA_UNADJ,
  meanDBH = weighted.mean(DBHcm, TPA_UNADJ, na.rm = TRUE)
), by = .(uniquePLOTid, uniqueSUBPid, year = INVYR)]

round(quantile(allspecies_subplotvars_dt$BA, seq(0,1,0.01), na.rm = TRUE),4)
round(quantile(allspecies_subplotvars_dt$N, seq(0,1,0.01), na.rm = TRUE),4)

table(tree_dt$INVYR, tree_dt$DIACHECK)
table(tree_dt$INVYR, tree_dt$CLIGHTCD)
table(tree_dt$INVYR, tree_dt$CCLCD)

treevars_dt <- tree_dt[,.(
  uniquePLOTid,
  uniqueSUBPid,
  uniqueTREEid,
  year = INVYR,
  species = SPCD,
  DBH = DBHcm,
  height = ft2m(HT),
  DIACHECK = DIACHECK, # 0 = accurate, 1 = estimated, 2 = deltawrong, 5 = modeled
  CLIGHTCD = CLIGHTCD, # 0 = shaded, 1 = top_or_1_side, 2 = top_and_1_side, 3 = top_and_2_sides, 4 = top_and_3_sides, 5 = full_light
  CCLCD = CCLCD, # 1 = open_grown, 2 = dominant, 3 = codominant, 4 = intermediate, 5 = overtopped
  HTCD = HTCD, # 1 = measured, 2 = totvisest-actmeas, 3 = visest, 4 = estimated
  ACTUALHTm = ft2m(ACTUALHT),
  TPA_UNADJ,
  status = STATUSCD # 0 = none, 1 = alive, 2 = dead, 3 = cut
  )]

treevars_dt <- treevars_dt[order(year)]

treevars_dt[
  order(year),
  ":="(
    status_before = data.table::shift(status,1,type = "lag")
    ),
  by = .(uniquePLOTid, uniqueSUBPid, uniqueTREEid)]
treevars_dt[
  order(year),
  ":="(
    DBH_before = data.table::shift(DBH,1,type = "lag")
    ),
  by = .(uniquePLOTid, uniqueSUBPid, uniqueTREEid)]

treevars_dt[, fresh_dead := status == 2 & status_before == 1,]
treevars_dt[, DBHchange := DBH - DBH_before,]

if(!dir.exists("data/calibration-data/FIA/prepared")) {
  dir.create("data/calibration-data/FIA/prepared")
}

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## save the full tables ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

fwrite(treevars_dt, "data/calibration-data/FIA/prepared/treevars_dt.csv")
fwrite(species_subplotvars_dt, "data/calibration-data/FIA/prepared/species_subplotvars.csv")
fwrite(allspecies_subplotvars_dt, "data/calibration-data/FIA/prepared/allspecies_subplotvars.csv")
fwrite(cohort_species_subplotvars_dt, "data/calibration-data/FIA/prepared/cohort_species_subplotvars.csv")
fwrite(state_dt, "data/calibration-data/FIA/prepared/state_dt.csv")
fwrite(plot_dt, "data/calibration-data/FIA/prepared/plot_dt.csv")
fwrite(species_dt, "data/calibration-data/FIA/prepared/species_dt.csv")


