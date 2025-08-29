## Preparing Oregon FIA data for ingestion into FINN model
## 08 August 2025
## D. Perret

## (1) Download FIA data for the state of Oregon

### using rFIA

library(rFIA)
library(dplyr)

dir.path = "data/FIA" # <- change this to match your file structure

# options(timeout=3600)
# rFIA::getFIA(states = "OR",
#              dir = dir.path,
#              common = TRUE)

### using the FIA DataMart

#### This is the way that I tend to access data. https://apps.fs.usda.gov/fia/datamart/datamart.html
#### Navigate to the state you want to download, and there's an option to download all state tables as .csv files in a zip archive.

## (2) Load FIA data into data object

fia <- readFIA(dir = dir.path,
               common = T,
               states = "OR")

## (3) Putting together tree data

tree.dat <- fia$TREE %>%
  filter(INVYR > 2000, # only take trees from annual inventories
         DIA > 5) %>%  # only count trees on subplots
  left_join(fia$PLOT %>%
              select(CN, MEASYEAR),
            by = c("PLT_CN" = "CN")) %>%
  mutate(treeID = paste(COUNTYCD, PLOT, SUBP, TREE, sep="_"), # unique IDs that are temporally consistent
         siteID = paste(COUNTYCD, PLOT,sep="_"),
         patchID = paste(COUNTYCD, PLOT,SUBP,sep="_"),
         dbh = DIA*2.54, # converting dbh from in to cm
         status = case_when(STATUSCD == 1 & is.na(PREV_TRE_CN) ~ "new", # trees that are alive and not recorded prior
                            STATUSCD == 1 & !is.na(PREV_TRE_CN) ~ "alive",
                            STATUSCD == 2 ~ "dead"),
         mort_cause = case_when(STATUSCD==2 & AGENTCD %in% c(00,70) ~ "unknown",
                                STATUSCD==2 & AGENTCD == 10 ~ "insect",
                                STATUSCD==2 & AGENTCD == 20 ~ "disease",
                                STATUSCD==2 & AGENTCD == 30 ~ "fire",
                                STATUSCD==2 & AGENTCD == 40 ~ "animal",
                                STATUSCD==2 & AGENTCD == 50 ~ "weather",
                                STATUSCD==2 & AGENTCD == 60 ~ "competition",
                                STATUSCD==2 & AGENTCD == 80 ~ "land use",
                                STATUSCD==2 & is.na(AGENTCD) &
                                  (PREV_STATUS_CD==1 |
                                     is.na(PREV_STATUS_CD)) ~ "unknown")) %>%
  select(siteID, patchID, treeID, year = MEASYEAR,
         SPCD, dbh, status, mort_cause) #note: species codes can be converted to full names using the species reference table found here: https://usfs-public.app.box.com/v/FIA-TreeSpeciesList/folder/290601914132

## (4) site data

site.dat <- fia$PLOT %>%
  mutate(siteID = paste(COUNTYCD, PLOT, sep="_")) %>%
  filter(INVYR > 2000,
         PLOT_STATUS_CD == 1, # plots that were actually measured
         siteID %in% tree.dat$siteID) %>% # just making sure data are consistent between tables
  select(siteID,
         year = MEASYEAR,
         x = LON,
         y = LAT)

library(terra)
source("R/createInputs.R")
library(FINN)
library(data.table)

species_dt <- as.data.table(readxl::read_xlsx("data/FIA/v9-5_2024-10_Natl_MasterTreeSpeciesList.xlsx"))
species_dt <- species_dt[,.(SPCD = `FIA Code`, species_name = SCIENTIFIC_NAME)]

tree.dat <- merge(tree.dat, species_dt, by = "SPCD", all.x = T)
fwrite(site.dat, "vignettes/FIAsubset_sites.csv")
fwrite(tree.dat, "vignettes/FIAsubset_trees.csv")

site_dt <- fread("vignettes/FIAsubset_sites.csv")
tree_dt <- fread("vignettes/FIAsubset_trees.csv")


# create env_dt

env_dt <- data.table(site.dat)
# env_dt[,year := as.integer(factor(year)), by = siteID]
site_pts <- vect(env_dt, geom = c("x", "y"), crs = "EPSG:4326")  # SpatVector points

temp <- terra::rast("../worldclim/wc2.1_10m_bio/wc2.1_10m_bio_1.tif")
prec <- terra::rast("../worldclim/wc2.1_10m_bio/wc2.1_10m_bio_12.tif")

v_proj <- project(site_pts, crs(temp))

env_dt$temp <- extract(temp, v_proj, method = "near")[,2]
env_dt$prec <- extract(prec, v_proj, method = "near")[,2]

# create obs_dt
# you must add flags to each tree observation for its status
# tree_dt <- data.table(tree.dat)
tree_dt[, status_before := data.table::shift(status,1,type = "lag"), by = treeID]
tree_dt[status == "new", reg := T,]
tree_dt[status == "dead" & status_before == "alive", mort := T,]
tree_dt[status == "alive" | status == "new", living := T,]
obs_dt <- makeObsData(tree_dt, plotsize = 0.06)

# select sites with 15 most common species
obs_dt[,Nspecies := uniqueN(species_name), by = siteID]
tree_dt[,Nspecies := uniqueN(species_name), by = siteID]

Nspecies = 5
# round(sort(table(tree_dt[status == "alive"]$species_name, useNA = "always"))/length(tree_dt[status == "alive"]$species_name),2)
# cumsum(round(sort(table(tree_dt[status == "alive"]$species_name), decreasing = T)/length(tree_dt[status == "alive"]$species_name),2))
selected_species <- names(sort(table(tree_dt[status == "alive"]$species_name, useNA = "always"), decreasing = T)[1:Nspecies])
# valid_species <- as.integer(names(sort(table(tree_dt[status == "alive"]$species), useNA = "always", decreasing = T)[1:Nspecies]))

obs_dt[, Nspecies_selected := sum(species_name %in% selected_species), by = siteID]
obs_dt[, Nspecies_notselected := sum(!(species_name %in% selected_species)), by = siteID]
obs_dt <- obs_dt[Nspecies_notselected == 0]

uniqueN(obs_dt$siteID)
selected_patchIDs <- obs_dt[,.(NyearsPerPatch = uniqueN(year)), by = .(siteID,patchID)][NyearsPerPatch == 3]$patchID
# selected_siteIDs <- sample(unique(obs_dt[Npatches == 4]$siteID),500)

obs_dt = obs_dt[patchID %in% selected_patchIDs]

inputDat <- resolveSiteIDs(site_dt = site_dt, tree_dt = tree_dt, obs_dt = obs_dt[patchID %in% selected_patchIDs], env_dt = env_dt[,.(siteID, year, temp, prec)], years2periods = T, Npatches_fixed = 4, selection_priority = "fixed patches", all = F, createInitCohorts = T)


table(inputDat$siteID_dt$year)
table(inputDat$tree_dt$year)
table(inputDat$site_dt$year)
table(inputDat$obs_dt$year)
# table(inputDat$env_dt)
# inputDat$initCohorts

m1 = finn(
  N_species = uniqueN(inputDat$obs_dt$species),recruits_dbh = 12.9,
  competition_process = createProcess(~0, func = FINN::competition, optimizeSpecies = TRUE),
  growth_process = createProcess(~temp+prec, FINN::growth, optimizeSpecies = TRUE, optimizeEnv = T),
  regeneration_process = createProcess(~temp+prec, FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = T),
  mortality_process = createProcess(~temp+prec, func = FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = T)
)

uniqueN(inputDat$obs_dt$siteID)
uniqueN(inputDat$env_dt$siteID)


# TODO put this shit into a function

obs_dt <- inputDat$obs_dt

# expand all combinations of siteID, patchID, year, species
# fill missing values for ba, trees, and reg with 0
# fill missing values for growth, dbh, mort with NA
obs_dt <- merge(
  obs_dt,
  data.table(expand.grid(
    siteID = unique(obs_dt$siteID),
    patchID = unique(obs_dt$patchID),
    year = unique(obs_dt$year),
    species = unique(obs_dt$species)
  )),
  by = c("siteID", "patchID", "year", "species"),
  all = TRUE
)
obs_dt[is.na(ba), ba := 0]
obs_dt[is.na(trees), trees := 0]
obs_dt[is.na(reg), reg := 0]
obs_dt[is.na(growth), growth := NA_real_]
obs_dt[is.na(mort), mort := NA_real_]
obs_dt[is.na(dbh), dbh := NA_real_]
obs_dt = obs_dt[order(siteID, patchID, species, year)]

obs_dt <- obs_dt[,.(
  ba = mean(ba, na.rm = T),
  trees = mean(trees, na.rm = T),
  dbh = mean(dbh, na.rm = T),
  growth = mean(growth, na.rm = T),
  mort = mean(mort, na.rm = T),
  reg = mean(reg, na.rm = T)
), by = .(siteID, year, species)]



m1$fit(
  env = inputDat$env_dt[,.(siteID,year,temp = as.numeric(scale(temp)),prec = as.numeric(scale(prec)))],
  data = obs_dt,
  init = inputDat$initCohorts,
  device = "cpu",
  epochs = 10,
  batchsize = 10L,
  patch_size = 0.06,
  patches = 4,
  lr = 0.01
)


