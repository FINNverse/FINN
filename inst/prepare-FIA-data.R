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

# # Create the data.table
# state_IDs <- data.table(
#   STATECD = c(1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29),
#   state_name = c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", "Delaware", "District of Columbia",
#                  "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana",
#                  "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri")
# )

# state_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_REF_STATE_ELEV.csv")

# STATUSCD_dt <- data.table(
#   0 = "none",
#   1 = "alive",
#   2 = "dead",
#   3 = "cut"
# )
#
# DIAHTCD_dt <- data.table(
#   1 = "breast",
#   2 = "root"
# )
#
# HTCD_dt <- data.table(
#   1 = "measured",
#   2 = "totvisest-actmeas",
#   3 = "visest",
#   4 = "estimed"
# )
# AGENTCD_dt <- data.table(
#   00 = "NoAgent",
#   10 = "Insect",
#   20 = "Disease",
#   30 = "Fire",
#   40 = "Animal",
#   50 = "Weather",
#   60 = "Vegetation",
#   70 = "Unknown",
#   80 = "Silvicultural"
#   )
#
# DIACHECK_dt <- data.table(
#   0 = "accurate",
#   1 = "estimated",
#   2 = "deltawrong",
#   5 = "modeled"
# )


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## select plots based on plot variables ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
plot_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv")
condition_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_COND.csv")

hist(condition_dt$BALIVE)
summary(condition_dt$BALIVE)

plot_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT),]
plot_dt[, uniquePLOTid := as.integer(factor(uniquePLOTid_txt, levels = unique(uniquePLOTid_txt)))]
plotID_dt <- unique(plot_dt[, .(uniquePLOTid, uniquePLOTid_txt)])
fwrite(plotID_dt, "data/calibration-data/FIA/plotID_dt.csv")
selected_plots <- plot_dt[
  PLOT_STATUS_CD == 1 & QA_STATUS == 1 & INTENSITY == 1 & DESIGNCD == 1 & INVYR != 9999
    ]

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

treevars_dt <- tree_dt[,.(
  DBH = DBHcm,
  status = STATUSCD,
  year = INVYR,
  uniquePLOTid,
  uniqueSUBPid,
  uniqueTREEid
  )]


# testIDs <- as.integer(names(table(treevars_dt[year < 1995,]$uniqueTREEid)[table(treevars_dt[year < 1995,]$uniqueTREEid) > 2]))
# treevars_dt[,.(Nyears = uniqueN(year)), .(uniqueTREEid)]
#
#
# treevars_test_dt <- treevars_dt[uniqueTREEid %in% testIDs]
# treevars_test_dt[,.(Nyears = uniqueN(year)),.(uniqueTREEid)]
#
# treevars_test_dt[status == 2]
treevars_dt[
  order(year),
  status_before := data.table::shift(status,1,type = "lag"),
  by = .(uniquePLOTid, uniqueSUBPid, uniqueTREEid)]

treevars_dt[, fresh_dead := status == 2 & status_before == 1,]

if(!dir.exists("data/calibration-data/FIA/prepared")) {
  dir.create("data/calibration-data/FIA/prepared")
}

fwrite(treevars_dt, "data/calibration-data/FIA/prepared/treevars_dt.csv")
fwrite(species_subplotvars_dt, "data/calibration-data/FIA/prepared/species_subplotvars.csv")
fwrite(allspecies_subplotvars_dt, "data/calibration-data/FIA/prepared/allspecies_subplotvars.csv")
fwrite(cohort_species_subplotvars_dt, "data/calibration-data/FIA/prepared/cohort_species_subplotvars.csv")


