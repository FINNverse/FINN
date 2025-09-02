############################################################
# Preparing Oregon FIA data for ingestion into FINN model
# Restructured to follow the workflow diagram
# 08 Aug 2025 — D. Perret (refactor)
############################################################

## 0) Setup ----
library(rFIA)
library(dplyr)
library(data.table)
library(terra)
library(readxl)
library(FINN)

source("R/createInputs.R")   # supplies makeObsData(), resolveSiteIDs(), createInitCohorts() if present

dir.path <- "data/FIA"       # <— adjust to your layout

## 1) Download / Load FIA (state = OR) ----
# rFIA::getFIA(states = "OR", dir = dir.path, common = TRUE)
fia <- readFIA(dir = dir.path, common = TRUE, states = "OR")

## 2) Build site_dt (from PLOT)  [site_dt] ----
site_dt <- fia$PLOT %>%
  mutate(siteName = paste(COUNTYCD, PLOT, sep = "_")) %>%
  filter(INVYR > 2000, PLOT_STATUS_CD == 1) %>%
  transmute(
    siteName,
    year = MEASYEAR,
    x = LON,
    y = LAT
  ) %>%
  distinct() %>%
  as.data.table()

## 3) Build raw tree_dt (from TREE + PLOT year)  [tree_dt] ----
tree_dt <- fia$TREE %>%
  filter(INVYR > 2000, DIA > 5) %>%            # annual program, subplots only
  left_join(fia$PLOT %>% select(CN, MEASYEAR),
            by = c("PLT_CN" = "CN")) %>%
  mutate(
    treeName  = paste(COUNTYCD, PLOT, SUBP, TREE, sep = "_"),
    siteName  = paste(COUNTYCD, PLOT, sep = "_"),
    patchName = paste(COUNTYCD, PLOT, SUBP, sep = "_"),
    year    = MEASYEAR,
    dbh     = DIA * 2.54                         # in → cm
  ) %>%
  select(siteName, patchName, treeName, year, SPCD, dbh, STATUSCD,
         PREV_TRE_CN, PREV_STATUS_CD, AGENTCD) %>%
  as.data.table()

## 3a) Attach species names (optional but handy) ----
species_dt <- as.data.table(read_xlsx("data/FIA/v9-5_2024-10_Natl_MasterTreeSpeciesList.xlsx"))[
  , .(SPCD = `FIA Code`, species_name = SCIENTIFIC_NAME)
]
tree_dt <- merge(tree_dt, species_dt, by = "SPCD", all.x = TRUE)

## 4) Define status flags on trees  [define status] ----
# (new/alive/dead + mort cause) — strictly on tree records
tree_dt[, status :=
          fifelse(STATUSCD == 1 & is.na(PREV_TRE_CN), "new",
                  fifelse(STATUSCD == 1 & !is.na(PREV_TRE_CN), "alive",
                          fifelse(STATUSCD == 2, "dead", NA_character_)))]

tree_dt[, mort_cause := fifelse(STATUSCD == 2 & AGENTCD %in% c(0, 70), "unknown",
                                fifelse(STATUSCD == 2 & AGENTCD == 10, "insect",
                                        fifelse(STATUSCD == 2 & AGENTCD == 20, "disease",
                                                fifelse(STATUSCD == 2 & AGENTCD == 30, "fire",
                                                        fifelse(STATUSCD == 2 & AGENTCD == 40, "animal",
                                                                fifelse(STATUSCD == 2 & AGENTCD == 50, "weather",
                                                                        fifelse(STATUSCD == 2 & AGENTCD == 60, "competition",
                                                                                fifelse(STATUSCD == 2 & AGENTCD == 80, "land use", NA_character_))))))))]

# convenience lag for obs building
setorder(tree_dt, treeName, year)
tree_dt[, status_before := shift(status, 1, type = "lag"), by = treeName]

## 5) Derive environmental covariates from site_dt  [site_dt ➜ get env ➜ env_dt] ----
# NOTE: env is computed per (siteID, year) as in the diagram.
env_dt <- copy(site_dt)   # keep siteID, year, x, y

site_pts <- vect(env_dt, geom = c("x", "y"), crs = "EPSG:4326")

# Example: WorldClim BIO1 (temp) & BIO12 (prec)
temp <- terra::rast("../worldclim/wc2.1_10m_bio/wc2.1_10m_bio_1.tif")
prec <- terra::rast("../worldclim/wc2.1_10m_bio/wc2.1_10m_bio_12.tif")

v_proj <- project(site_pts, crs(temp))
env_dt[, temp := extract(temp, v_proj, method = "near")[, 2] ]
env_dt[, prec := extract(prec, v_proj, method = "near")[, 2] ]
env_dt <- env_dt[, .(siteName, year, temp, prec)]

## 6) Create observation table  [tree_dt ➜ create obs_dt] ----
# flags required by makeObsData(): reg/mort/living per (tree, year)
td <- copy(tree_dt)
td[status == "new",                          reg := TRUE]
td[status == "dead" & status_before == "alive", mort := TRUE]
td[status %chin% c("alive", "new"),          living := TRUE]

# Aggregate to FINN obs structure (ba/trees/dbh/growth/mort/reg per site/patch/species/year)
sampled_siteNames <- sample(unique(td$siteName), 100)
obs_list <- makeObsData(tree_dt = td[siteName %in% sampled_siteNames], plotsize = 0.06, aggregate_by_site = F, Nspecies = 0, Npatches = 4)

obs_dt = obs_list$obs_dt
tree_dt = obs_list$tree_dt

## 7) Select sites (species & patch consistency)  [select sites] ----
# Example selection: keep top N species and patches with full time coverage

## 9) Resolve IDs & finalize model inputs  [resolve IDs ➜ {env_dt, obs_dt, init_dt}] ----
# This step aligns site/patch/species IDs across env_dt, obs_dt, init_dt and optionally
# converts years→periods and fixes the number of patches, as in your original call.
inputs <- resolveSiteIDs(
  tree_dt            = as.data.table(tree_dt),
  obs_dt             = as.data.table(obs_dt),
  env_dt             = as.data.table(env_dt),
  createInitCohorts = F
)

tree_dt <- inputs$tree_dt[year != 0]
obs_dt <- inputs$obs_dt[year != 0]
env_dt <- inputs$env_dt[year != 0, .(siteID, year, temp = as.numeric(scale(temp)),prec = as.numeric(scale(prec)))]
init_trees <- inputs$tree_dt[year == 0 & living == T]
init_cohorts <- makeInitCohorts(init_trees, Nspecies = max(init_trees$species))

m1 <- finn(
  N_species            = uniqueN(obs_dt$species),
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0,     FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createProcess(~1+temp, FINN::growth,   optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~1+temp, FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~1+temp, FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)

m1$fit(
  env        = env_dt[,-c("prec")],
  data       = obs_dt,
  init       = init_cohorts,
  device     = "cpu",
  epochs     = 500,
  batchsize  = 10L,
  patch_size = 0.06,
  patches    = 4,
  lr         = 0.0000001
)

str(m1$raw_g)


## (Optional) Minimal model to debug growth/mortality structure ----
# m1 <- finn(
#   N_species = uniqueN(obs_dt$species), recruits_dbh = 12.9,
#   competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
#   growth_process       = createProcess(~1, FINN::growth,      optimizeSpecies = TRUE, optimizeEnv = TRUE),
#   regeneration_process = createProcess(~1, FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
#   mortality_process    = createProcess(~1, FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
# )
# m1$fit(env = env_dt[, .(siteID, year, temp, prec)], data = obs_dt, init = init_dt, ...)

## 12) Save (optional) ----
fwrite(site_dt, "vignettes/FIAsubset_sites.csv")
fwrite(tree_dt, "vignettes/FIAsubset_trees.csv")
fwrite(env_dt,  "vignettes/FIAsubset_env.csv")
fwrite(obs_dt,  "vignettes/FIAsubset_obs.csv")
if (!is.null(init_dt)) fwrite(init_dt, "vignettes/FIAsubset_init.csv")
