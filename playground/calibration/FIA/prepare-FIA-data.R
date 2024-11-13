#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Project: FINN
## Script purpose:  Explore FIA data for calibration of mortality parameters.
## Date: Mon Jul  8 09:17:04 2024
## Author: Yannek Kaeber <y.kaeber@posteo.de>
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

library(data.table)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## functions and constants ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# Constants for TPA_UNADJ for each plot type
TPA_SUBPLOT <- 6.018046
TPA_MICROPLOT <- 74.965282
TPA_MACROPLOT <- 0.999188

# Function to calculate plot size in hectares
calculate_plot_size <- function(tpa_unadj) {
  plot_size_acres <- 1 / tpa_unadj
  plot_size_hectares <- plot_size_acres * 0.404686
  return(plot_size_hectares)
}

# Calculating plot sizes in hectares for subplots and microplots
subplot_size_ha <- calculate_plot_size(TPA_SUBPLOT)
microplot_size_ha <- calculate_plot_size(TPA_MICROPLOT)
macroplot_size_ha <- calculate_plot_size(TPA_MACROPLOT)

circle_area <- function(d) pi * (d/100 / 2)^2
inches2cm <- function(x) x * 2.54
ft2m <- function(x) x * 0.3048
sqrft2m2 <- function(x) x * 0.092903
acre2ha <- function(x) x * 0.404686
sqrft2acre <- function(x) x * 2.29568e-5
radius2area <- function(x) pi * x^2

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
survey_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_SURVEY.csv")
plot_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv")
subplot_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_SUBPLOT.csv")
subp_cond_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_SUBP_COND.csv")
cond_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_COND.csv")
# pop_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_POP_ESTN_UNIT.csv")
# subpcond_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_SUBP_COND.csv")

# create uniquePLOTid based on "natural" unique variables for PLOT TABLE
# see FIA documentation for details
plot_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT), by = INVYR]
plot_dt[, uniquePLOTid := as.integer(factor(uniquePLOTid_txt, levels = unique(uniquePLOTid_txt)))]
subplot_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT),]
subplot_dt[, uniqueSUBPid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT, "_", SUBP),]
cond_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT),]
subp_cond_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT),]
subp_cond_dt[, uniqueSUBPid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT, "_", SUBP),]

# save table with unique IDs for linking climate variable extraction to the correct PLOT
plotID_dt <- unique(plot_dt[, .(uniquePLOTid, uniquePLOTid_txt)])

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## select inventories based on attributes in INVENTORY ####
# here the unique link to PLOT is SRV_CN
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

selected_survey_dt <- survey_dt[
  INVYR != 9999
  ]

table(selected_survey_dt$CYCLE, selected_survey_dt$INVYR)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## select in PLOT ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# select plots based on PLOT_STATUS_CD, QA_STATUS, INTENSITY, DESIGNCD, INVYR
# PLOTSTATUS_CD = 1 --> Only sampled plots are included
# QA_STATUS = 1 --> Only plots that have passed the quality assurance. 1 = Standard production plot.
# INTENSITY = 1 --> Only plots that are part of the FIA standard grid
# DESIGNCD = 1 --> Only plots that follow the "harmonized" set of the FIA sampling
# INVYR != 9999 --> Only plots that have a valid inventory year
selected_plots_dt <- plot_dt[
  SRV_CN %in% selected_survey_dt$CN &
  PLOT_STATUS_CD == 1 &
    QA_STATUS == 1 &
    INTENSITY == 1 &
    DESIGNCD == 1 &
    INVYR != 9999
    ,
  .(
    uniquePLOTid,
    uniquePLOTid_txt,
    INVYR,
    MEASYEAR,
    MEASMON,
    MEASDAY,
    PREV_PLT_CN, # Link to CN of perivous inventory
    LAT,
    LON,
    ELEV,
    ECOSUBCD, # use for selecting ecoregions of interest (see "data/calibration-data/FIA/73326-wo-gtr-76d-cleland2007.pdf)
    MANUAL # Version number of field guide
    )
  ]

# test_dt = selected_plots_dt[, .(
#   minyear = min(INVYR),
#   maxyear = max(INVYR),
#   all_years = paste0(sort(unique(INVYR)), collapse = ","),
#   Nyears = uniqueN(INVYR)
# ) ,.(uniquePLOTid_txt)]
# test_dt[Nyears > 1]
# table(test_dt$Nyears)


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## select in COND####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# COND_STATUS_CD = 1 --> Sampled - at least one accessible forest land condition present on subplot.
# RESERVCD = 1 -->  Reserved land is permanently prohibited from being managed for the production of
#                   wood products through statute or agency mandate; the prohibition cannot be changed
#                   through decision of the land manager. Logging may occur to meet protected area
#                   objectives. Examples include designated Federal wilderness areas, national parks
#                   and monuments, and most State parks.
# TRTCD1 = 0 --> No treatment

cond_dt[,CONDID2 := as.integer(factor(INVYR, ordered = T)), by = uniquePLOTid_txt]
cond_dt <- cond_dt[CONDID == CONDID2]

cond_dt[DSTRBCD1 != DSTRBCD1_P2A & DSTRBCD1 == 0, ":="(
  DSTRBCD1 = DSTRBCD1_P2A,
  DSTRBYR1 = DSTRBYR1_P2A
)]
cond_dt[DSTRBCD2 != DSTRBCD1_P2A & DSTRBCD2 == 0, ":="(
  DSTRBCD2 = DSTRBCD2_P2A,
  DSTRBYR2 = DSTRBYR2_P2A
)]
cond_dt[DSTRBCD3 != DSTRBCD1_P2A & DSTRBCD3 == 0, ":="(
  DSTRBCD3 = DSTRBCD3_P2A,
  DSTRBYR3 = DSTRBYR3_P2A
)]

table(cond_dt$STUMP_CD_PNWRS, useNA = "always")

table(cond_dt$TRTCD1, useNA = "always")

cond_dt[,minyear := min(INVYR), by = uniquePLOTid_txt]
cond_dt[,perdiodID := as.integer(factor(INVYR,ordered = T)), by = uniquePLOTid_txt]


cond_dt[is.na(TRTCD1), TRTCD1 := 999,]

cond_dt[perdiodID == 1, TRTCD1 := 0,]
cond_dt[perdiodID == 1, TRTCD1_P2A := 0,]
cond_dt[, anyTRTCD1not0 := any(TRTCD1 != 0), by = uniquePLOTid_txt]
cond_dt[, anyTRTCD1_P2Anot0 := any(TRTCD1_P2A != 0), by = uniquePLOTid_txt]

selected_cond_dt <- cond_dt[
  uniquePLOTid_txt %in% unique(selected_plots_dt$uniquePLOTid_txt) &
    COND_STATUS_CD == 1 &
    # RESERVCD == 1 &
    anyTRTCD1not0 == F &
    STUMP_CD_PNWRS != "N"
  ,
  .(
    uniquePLOTid_txt,
    RESERVCD,
    CONDID,
    CONDID2,
    INVYR,
    FORTYPCD, # Forest type code derived (see appendix D in FIA documentation)
    FLDTYPCD, # Field type code from field team (see appendix D in FIA documentation)
    PROP_BASIS, # Proportion basis. A value indicating what type of fixed-size subplots were installed when this plot was sampled. This information is needed to use the proper adjustment factor for the stratum in which the plot occurs (see POP_STRATUM.ADJ_FACTOR_SUBP and POP_STRATUM.ADJ_FACTOR_MACR).
    SLOPE, # percentage
    ASPECT, # degree
    PHYSCLCD, # Physiographic class code. A code indicating the general effect of land form, topographical position, and soil on moisture available to trees. (see FIA documentation)
    DSTRBCD1, # Disturbance code 1. A code indicating the kind of disturbance occurring since the last measurement or within the last 5 years for new plots. The area affected by the disturbance must be at least 1 acre in size.
    DSTRBYR1,
    DSTRBCD2, # Disturbance code 2. A code indicating the kind of disturbance occurring since the last measurement or within the last 5 years for new plots. The area affected by the disturbance must be at least 1 acre in size.
    DSTRBYR2,
    DSTRBCD3, # Disturbance code 3. A code indicating the kind of disturbance occurring since the last measurement or within the last 5 years for new plots. The area affected by the disturbance must be at least 1 acre in size.
    DSTRBYR3,
    BALIVE # Basal area per acre of live trees. Basal area in square feet per acre of all live trees 1.0 inch d.b.h./d.r.c. sampled in the condition. (ft2/ac)
  )
]

# selected_cond_dt <- selected_cond_dt[,-"V1"]
selected_cond_dt <- merge(selected_plots_dt, selected_cond_dt, by = c("uniquePLOTid_txt", "INVYR"))

# test_dt = selected_cond_dt[, .(
#   minyear = min(INVYR),
#   maxyear = max(INVYR),
#   all_years = paste0(sort(unique(INVYR)), collapse = ","),
#   Nyears = uniqueN(INVYR)
# ) ,.(uniquePLOTid_txt)]
# test_dt[Nyears > 1]
# table(test_dt$Nyears, useNA = "always")
#
# table(cond_dt$perdiodID, cond_dt$TRTCD1, useNA = "always")
# round(table(cond_dt$perdiodID, cond_dt$TRTCD1, useNA = "always")/sum(table(cond_dt$perdiodID, cond_dt$TRTCD1, useNA = "always")),2)
# round(table(cond_dt$perdiodID, cond_dt$TRTCD1_P2A, useNA = "always")/sum(table(cond_dt$perdiodID, cond_dt$TRTCD1_P2A, useNA = "always")),2)
# round(table(cond_dt$TRTCD1, cond_dt$TRTCD1_P2A, useNA = "always")/sum(table(cond_dt$TRTCD1, cond_dt$TRTCD1_P2A, useNA = "always")),2)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## select SUBPLOT####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# SUBP = 1,2,3,4 --> Only subplots 1,2,3,4 are included. These correspond to the FIA standard grid (see FIA documentation)
# SUBP_STATUS_CD = 1 --> Sampled - at least one accessible forest land condition present on subplot.
selected_subplots_dt <- subplot_dt[
  uniquePLOTid_txt %in% unique(selected_plots_dt$uniquePLOTid_txt) &
    SUBP %in% c(1,2,3,4) &
    SUBP_STATUS_CD == 1 &
    MICRCOND == 1,
    SUBPCOND == 1,.(
      uniquePLOTid_txt,
      uniqueSUBPid_txt,
      SUBP,
      INVYR
    )
  ]

selected_subplots_dt <- merge(selected_cond_dt, selected_subplots_dt, by = c("uniquePLOTid_txt", "INVYR"))

# test_dt = selected_subplots_dt[, .(
#   minyear = min(INVYR),
#   maxyear = max(INVYR),
#   all_years = paste0(sort(unique(INVYR)), collapse = ","),
#   Nyears = uniqueN(INVYR)
# ) ,.(uniquePLOTid_txt)]
# test_dt[Nyears > 1]
# table(test_dt$Nyears, useNA = "always")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## select SUBP_COND####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
subp_cond_dt[,CONDID2 := as.integer(factor(INVYR, ordered = T)), by = uniqueSUBPid_txt]
subp_cond_dt[CONDID2 > CONDID, CONDID2 := max(CONDID), by = uniqueSUBPid_txt]
subp_cond_dt <- subp_cond_dt[CONDID == CONDID2]

# SUBPCOND_PROP = 1 --> Only subplots with a proportion of 1 are included. This means that the subplot fully meets the conditions defined in SUBPLOT table (which is SUBCOND == 1). In simple words, only plots that are fully covered by forest are included.
selected_subp_cond_dt <- subp_cond_dt[
  uniqueSUBPid_txt %in% unique(selected_subplots_dt$uniqueSUBPid_txt) &
    MICRCOND_PROP == 1 &
    SUBPCOND_PROP == 1,.(
      uniquePLOTid_txt,
      uniqueSUBPid_txt,
      INVYR
    )
  ]

final_selected_plots <- merge(selected_subplots_dt, selected_subp_cond_dt, by = c("uniquePLOTid_txt", "uniqueSUBPid_txt", "INVYR"))

# SUBP_COND_CHNG_MTRX is not relevant here because only subplots with SUBPCOND_PROP == 1 are included.

if(!dir.exists("data/calibration-data/FIA_V2")) {
  dir.create("data/calibration-data/FIA_V2")
}
fwrite(final_selected_plots, "data/calibration-data/FIA_V2/selected_plots.csv")

# remove all objects but the final_selected_plots
rm(list = setdiff(ls(), "final_selected_plots"))
gc()
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## select in TREE ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# read TREE table
tree_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_TREE.csv")
tree_thresh_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_TREE_GRM_THRESHOLD.csv")
tree_dt <- merge(tree_dt, tree_thresh_dt[,.(TRE_CN, THRESHOLD_TYPE, DIA_THRESHOLD = DIA)], by.x = c("CN"), by.y = c("TRE_CN"), all.x = T)
tree_dt[, uniquePLOTid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT),]
tree_dt[, uniqueSUBPid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT, "_", SUBP),]
tree_dt[, uniqueTREEid_txt := paste0(STATECD, "_", UNITCD, "_", COUNTYCD, "_", PLOT, "_", SUBP, "_", TREE),]
tree_dt[, uniqueTREEid := as.integer(factor(uniqueTREEid_txt, levels = unique(uniqueTREEid_txt)))]

## STATUSCD
# 0 = no status
# 1 = alive
# 2 = dead
# 3 = removed (cut)
## PREV_STATUS_CD
# 1 = live tree
# 2 = dead tree
## AGENTCD
# 00 = no agent recorded (only allowed on live trees in data prior to 1999)
# 10 = insect
# 20 = disease
# 30 = fire
# 40 = animal
# 50 = weather
# 60 = vegetation (e.g., suppression, competition, vines/kudzu)
# 70 = unknown/not sure/other (includes death from human activity not related to silvicultural or landclearing activity)
# 80 = silvicultural or landclearing activity (death caused by harvesting or other activity like girdling, chaining, etc.)
## MORTCD
# 0 = does not qualify as mortality
# 1 = does qualify as mortality
## STANDING_DEAD_CD
# 0 = no tree does not qualify as standing dead
# 1 = yes tree does qualify as standing dead
## RECONCILECD
# 1 = ingrowth (new tally tree not qualifying as through growth or a new tree on previously nonforest land)
# 2 = through growth (new tally tree ≥5 inches d.b.h. that was not missed at the previous inventory)
# 3 = retired code for missed live (live tree missed at previous inventory, now invalid for PLOT.MANUAL ≥ 9.0)
# 4 = retired code for missed dead (dead tree missed at previous inventory, now invalid for PLOT.MANUAL ≥ 9.0)
# 5 = shrank (tree that shrunk below threshold diameter, only currently live)
# 6 = physical movement (tree moved beyond plot radius or previously outside radius but now within due to natural causes)
# 7 = cruiser error (tree erroneously tallied or missed at the previous inventory)
# 8 = procedural change (tree previously tallied in error or missed, now tallied due to procedural change)
# 9 = nonsampled area (tree tallied in a sampled condition at previous inventory, now in a nonsampled condition)
## DIA [inch]
## PREVDIA [inch]
# Previous diameter. The previous diameter, in inches, of the sample tree at the point of diameter measurement. Populated for remeasured trees.
## DIAHTCD
# 1 = breast height
# 2 = root collar
## DIACHECK
# 0 = diameter accurately measured
# 1 = diameter estimated
# 2 = diameter measured at different location than previous measurement (remeasurement trees only)
# 5 = diameter modeled in the office (used with periodic inventories)
## HT [inch]
# Total height of the tree. Total height refers to the height the tree would have if broken
# parts of the ground were added.
## ACTUALHT
# Actual height. (All live and standing dead tally trees 1.0 inch d.b.h./d.r.c.)
# The length (height) of the tree to the nearest foot from ground level to the highest
# remaining portion of the tree still present and attached to the bole. If ACTUALHT = HT,
# then the tree does not have a broken top. If ACTUALHT <HT, then the tree does have a
# broken or missing top.
## HEIGHTCD
# 1 = field measurement
# 2 = total length visible estimated - actual length measured
# 3 = visual estimate
# 4 = estimated with a model
## CCLCD
# Crown class code. A code indicating the amount of sunlight received and the crown
# position within the canopy.
# 1 = open crown
# 2 = dominant
# 3 = codominant
# 4 = intermediate
# 5 = overtopped
## TPA_UNADJ
# Trees per acre unadjusted. The number of trees per acre that the sample tree theoretically
# represents based on the sample design. For fixed-radius plots taken with the mapped plot
# design (PLOT.DESIGNCD = 1), TPA_UNADJ is set to a constant derived from the plot size and
# equals 6.018046 for trees sampled on subplots, 74.965282 for trees sampled on microplots,
# and 0.999188 for trees sampled on macroplots. Variable-radius plots were often used in
# earlier inventories, so the value in TPA_UNADJ decreases as the tree diameter increases.
# Based on the procedures described in Bechtold and Patterson (2005), this attribute must
# be adjusted using factors stored in the POP_STRATUM table to derive population estimates.
# Examples of estimating population totals are shown in The Forest Inventory and Analysis
# Database: Population Estimation User Guide.
selected_tree_dt1 <- tree_dt[
  TPA_UNADJ == 6.018046 # only pick trees sampled on the subplot
  ,.(
    uniquePLOTid_txt,
    uniqueSUBPid_txt,
    uniqueTREEid,
    uniqueTREEid_txt,
    INVYR,
    CN,
    PREV_TRE_CN,
    STATUSCD,
    MORTCD,
    MORTYR,
    AGENTCD,
    STANDING_DEAD_CD,
    SPCD,
    DIA,
    PREVDIA,
    DIAHTCD,
    DIACHECK,
    THRESHOLD_TYPE,
    DIA_THRESHOLD,
    HT,
    ACTUALHT,
    HTCD,
    CCLCD,
    RECONCILECD,
    TPA_UNADJ
  )]

selected_tree_dt2 <- merge(
  selected_tree_dt1, final_selected_plots,
  by = c("uniquePLOTid_txt", "uniqueSUBPid_txt", "INVYR")
)


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## select in TREE_GRM_COMPONENT ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# in this table we select the derived growth and mortality

## ANN_DIA_GROWTH [inch]
# Computed annual diameter growth. The annual diameter growth for the tree expressed as inches per year.
## ANN_HT_GROWTH [inch]
# Computed annual height growth. The annual height growth for a tree expressed as feet per year.
## MICR_COMPONENT_AL_FOREST
# CUT0 = tree killed due to harvesting by T2, applicable in periodic-to-periodic, periodic-to-annual, and modeled GRM estimates
# CUT1 = tree previously in estimate at T1 and killed due to harvesting by T2
# CUT2 = tree grew across threshold diameter since T1 and was killed due to harvesting by T2
# INGROWTH = tree grew across threshold diameter for the estimate since T1
# MORTALITY0 = tree died of natural causes by T2, applicable in periodic-to-periodic, periodic-to-annual, and modeled GRM estimates
# MORTALITY1 = tree was previously in estimate at T1 and died of natural causes by T2
# NOT USED = tree was either live or dead at T1 and has no status at T2
# SURVIVOR = tree remained live and in estimate from T1 through T2
# UNKNOWN = tree lacks information required to classify, usually due to procedural changes
# REVERSION1 = tree grew across threshold diameter for the estimate by midpoint of measurement interval and condition reverted to land basis by T2
# REVERSION2 = tree grew across threshold diameter after midpoint and condition reverted to land basis by T2
# DIVERSION0 = tree removed from estimate by non-harvesting action by T2, in periodic-to-periodic, periodic-to-annual, and modeled GRM estimates
# DIVERSION1 = tree previously in estimate at T1 and condition diverted from land basis by T2
# DIVERSION2 = tree grew across threshold diameter since T1 and condition diverted from land basis by T2
# CULLINC = not used at this time
# CULLDEC = not used at this time
# N/A - A2A = component not defined or does not exist, only in annual-to-annual GRM estimates
# N/A - A2A SOON = component change not defined or does not exist, only in annual-to-annual GRM estimates
# N/A - MODELED = component change not defined or does not exist, only in annual-to-annual GRM estimates
# N/A - P2A = component change not defined or does not exist, only in periodic-to-annual GRM estimates
# N/A - P2P = component change not defined or does not exist, only in periodic-to-periodic GRM estimates
# N/A - PERIODIC = component change not defined or does not exist, only in periodic-to-periodic GRM estimates
## SUBP_COMPONENT_AL_FOREST
# see MICR_COMPONENT_AL_FOREST

# tree_grm_component_dt <- fread("data/calibration-data/FIA/CSV_FIADB_ENTIRE/ENTIRE_TREE_GRM_COMPONENT.csv")
# selected_tree_dt3 <- merge(
#   selected_tree_dt2,
#   tree_grm_component_dt[,.(
#     CN = TRE_CN,
#     ANN_DIA_GROWTH,
#     ANN_HT_GROWTH,
#     MICR_COMPONENT_AL_FOREST,
#     SUBP_COMPONENT_AL_FOREST
#     )], by = c("CN"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## calculate variables ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# complete time series for each tree
plot_eval_dt = selected_tree_dt1[, .(
  minyear = min(INVYR),
  maxyear = max(INVYR),
  all_years = paste0(sort(unique(INVYR)), collapse = ","),
  Nyears = uniqueN(INVYR)
  ) ,.(uniqueSUBPid_txt)]

# only keep plots with at least 2 years of data
plot_eval_dt <- plot_eval_dt[Nyears > 1]
# select plots
selected_tree_dt3 <- selected_tree_dt2[uniqueSUBPid_txt %in% unique(plot_eval_dt$uniqueSUBPid_txt),]

# create a table with all combinations of SUBPid and INVYR
all_SUBPid_comb_dt <- unique(selected_tree_dt3[,.(uniqueSUBPid_txt, uniquePLOTid, uniquePLOTid_txt, INVYR)])

# create a table with all combinations of TREEid and SUBPid
all_TREEid_dt <- unique(selected_tree_dt3[,.(uniqueTREEid_txt, uniqueSUBPid_txt, uniquePLOTid, uniquePLOTid_txt, uniqueTREEid)])

# merge all combinations of SUBPid and INVYR with all combinations of TREEid and SUBPid
all_SUBPid_comb_dt <- merge(all_SUBPid_comb_dt, all_TREEid_dt, by = c("uniqueSUBPid_txt","uniquePLOTid","uniquePLOTid_txt"), allow.cartesian = T)

# add all combinations of SUBPid and INVYR to selected_tree_dt3 so that each uniqueTREEid_text has each INVYR and SUBPid combination
selected_tree_dt4 <- merge(all_SUBPid_comb_dt, selected_tree_dt3, by = c("uniqueSUBPid_txt", "uniqueTREEid_txt", "uniquePLOTid", "uniquePLOTid_txt", "uniqueTREEid", "INVYR"), all.x = T)

# label new rows
selected_tree_dt4[is.na(CN), missingrow := T,]
selected_tree_dt4[!is.na(CN), missingrow := F,]

# only keep plots with all trees DIACHECK == 1
selected_tree_dt4[missingrow == F, DIACHECKisNot0 := sum(DIACHECK != 0 & STATUSCD == 1), by = uniqueTREEid_txt]
selected_tree_dt4[missingrow == F, STATUSCDis3 := sum(STATUSCD == 3), by = uniqueTREEid_txt]
selected_tree_dt5 <- selected_tree_dt4[
  !(uniqueSUBPid_txt %in% selected_tree_dt4[
    DIACHECKisNot0 > 0 | STATUSCDis3 > 0 | DIA_THRESHOLD > 5
    ]$uniqueSUBPid_txt),
  ]

## label time series
# derive full time series for each plot. fill trees with no information in first year
# add INVYR before and after
selected_tree_dt5[, INVYR_before := data.table::shift(INVYR, 1, type = "lag"), by = uniqueTREEid_txt]
selected_tree_dt5[, INVYR_after := data.table::shift(INVYR, 1, type = "lead"), by = uniqueTREEid_txt]
selected_tree_dt5[, INVYR_diff := INVYR - INVYR_before, ]
# add status for not existing trees with STATUSCD = 9
selected_tree_dt5[missingrow == T & is.na(STATUSCD), STATUSCD := 9,]
# add STATUSCD before and after
selected_tree_dt5[, STATUSCD_before := data.table::shift(STATUSCD, 1, type = "lag"), by = uniqueTREEid_txt]
selected_tree_dt5[, STATUSCD_after := data.table::shift(STATUSCD, 1, type = "lead"), by = uniqueTREEid_txt]
# add DIA before and after
selected_tree_dt5[, DIA_before := data.table::shift(DIA, 1, type = "lag"), by = uniqueTREEid_txt]
selected_tree_dt5[, DIA_after := data.table::shift(DIA, 1, type = "lead"), by = uniqueTREEid_txt]
# calculate absolute growth difference
selected_tree_dt5[, DIA_DIFF := DIA - DIA_before, ]
# calculate annual growth
selected_tree_dt5[STATUSCD == 1, DIA_GROWTH := DIA_DIFF/INVYR_diff, ]

# revive "dead" trees i.e. tree that become alive after beeing labeled dead
selected_tree_dt5[STATUSCD == 2 & STATUSCD_after == 1, STATUSCD := 1,]
# update STATUSCD before and after including the corrections
selected_tree_dt5[, STATUSCD_before := data.table::shift(STATUSCD, 1, type = "lag"), by = uniqueTREEid_txt]
selected_tree_dt5[, STATUSCD_after := data.table::shift(STATUSCD, 1, type = "lead"), by = uniqueTREEid_txt]

# define recruitment events
selected_tree_dt5[, recruitment := STATUSCD == 1 & !(RECONCILECD %in% c(3,7,8)) & STATUSCD_before == 9,]

# define mortality events
selected_tree_dt5[, mortality := (STATUSCD == 2 & STATUSCD_before == 1),]

table(selected_tree_dt5$mortality, useNA = "always")
table(selected_tree_dt5$recruitment, useNA = "always")
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## check errors ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
selected_tree_dt5[recruitment == T, RECR_GROWTH := (DIA-5)/INVYR_diff, ]

# 0.1 % of plots with trees that have extreme growth are removed
min_growth = quantile(selected_tree_dt5$DIA_GROWTH, c(0.005), na.rm = T)
max_growth = quantile(selected_tree_dt5$DIA_GROWTH, c(0.995), na.rm = T)
selected_tree_dt6 <- selected_tree_dt5[
  !(uniqueSUBPid_txt %in% unique(selected_tree_dt5[DIA_GROWTH < min_growth & DIA_GROWTH > max_growth & RECR_GROWTH > max_growth,]$uniqueSUBPid_txt)),
]
# remaining trees with growth < 0 are set to 0
selected_tree_dt6[DIA_GROWTH < 0, DIA_GROWTH := 0,]



## height
# HEIGHTCD == 1

selected_tree_dt6$plot_size = subplot_size_ha

selected_tree_dt6[,":="(
  DBHcm = inches2cm(DIA),
  height = ft2m(HT)
),]
selected_tree_dt6[,ba := round(circle_area(DBHcm), 4),]

plotvars_dt <- selected_tree_dt6[,.(
  N = sum((STATUSCD == 1)/plot_size, na.rm = T),
  BA = sum(ba*(STATUSCD == 1)/plot_size, na.rm = T)
), by = .(uniqueSUBPid_txt, INVYR)]

plot(plotvars_dt$N, plotvars_dt$BA)

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
  actualHeight = ft2m(ACTUALHT),
  HTCD = HTCD,
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

outputDir <- "data/calibration-data/FIA/preparedV3/"
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
