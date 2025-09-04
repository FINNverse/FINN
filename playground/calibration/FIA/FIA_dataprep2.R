############################################################
# Preparing Oregon FIA data for ingestion into FINN model
# Restructured to follow the workflow diagram
# 08 Aug 2025 — D. Perret (refactor)
############################################################

## 0) Setup ----
library(dplyr)
library(data.table)
library(terra)
library(readxl)
library(FINN)

source("R/createInputs.R")   # supplies makeObsData(), resolveSiteIDs(), createInitCohorts() if present

if(!file.exists("data/FIA/FIA_trees_raw.csv")){

  library(rFIA)
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

  #write output so that the data can be read without rFIA package
  fwrite(tree_dt, "data/FIA/FIA_trees_raw.csv")
  fwrite(site_dt, "data/FIA/FIA_sites_raw.csv")
}else{
  tree_dt <- fread("data/FIA/FIA_trees_raw.csv")
  sites_dt <- fread("data/FIA/FIA_sites_raw.csv")
}

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
# sampled_siteNames <- sample(unique(td$siteName), 500)
add_site_completeness_flags(env_dt, td)

obs_list <- makeObsData(
  tree_dt = td[complete == T],
  plotsize = 0.06,
  aggregate_by_site = F,
  # NspeciesQuantile = 0.9,
  Nspecies = 5,
  Npatches = 4, minNyears = 3)

# obs_list$obs_dt[,.(Nyear = uniqueN(year)), by = siteName]
set.seed(2)
obs_dt = obs_list$obs_dt[siteName %in% sample(unique(obs_list$obs_dt$siteName), 50)]
# obs_dt[,.(Nyear = uniqueN(year)), by = siteName]
tree_dt = obs_list$tree_dt[siteName %in% unique(obs_dt$siteName)]
# tree_dt[,.(Nyear = uniqueN(year)), by = siteName]

uniqueN(tree_dt$species_name)
unique(tree_dt$species_name)
unique(obs_dt$species_name)

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
inputs$obs_dt_patches[siteID == 1]
tree_dt <- inputs$tree_dt[year != 0]
obs_dt <- inputs$obs_dt[year != 0]
uniqueN(inputs$tree_dt[year != 0]$species)
uniqueN(inputs$obs_dt[year != 0]$species)
env_dt <- unique(inputs$env_dt[year != 0, .(siteID, year, temp = as.numeric(scale(temp)),prec = as.numeric(scale(prec)))])
init_trees <- inputs$tree_dt[year == 0]
init_trees[living == F, trees := 0,]
init_cohorts <- makeInitCohorts(init_trees, Nspecies = max(obs_dt$species))

m1 <- finn(
  N_species            = uniqueN(obs_dt$species),
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0,     FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createProcess(~1+temp+prec, FINN::growth,   optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~1+temp+prec, FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~1+temp+prec, FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)

# first fit ####
m1$fit(
  env        = env_dt,
  data       = unique(obs_dt),
  init_cohort = init_cohorts,
  device     = "cpu",
  epochs     = 100,
  batchsize  = 50L,
  patch_size = 0.06,
  patches    = 4, weights = c(0.1, 10, 1.0, 10.0, 1, 1),
  lr         = 0.01#, loss = c("mse","mse","mse","mse","mse","mse","mse")
)

par1 = m1$parameters_r

# simulate ####
env_dt2 <- env_dt[year == 1]
for(i in 2:200){
  env_dt2 <- rbind(env_dt2, env_dt[year == 1,.(siteID, year = i, prec, temp)])
}
species_dt <- unique(inputs$obs_dt[,.(species, species_name)])
sim2 <- m1$simulate(env_dt2, init_cohort = NULL)
p_dat <- sim2$long$site[, .(value = mean(value)), by = .(year, species, variable)]
p_dat <- merge(p_dat, species_dt, by = "species", all.x = TRUE)
p_dat[variable %in% c("ba","trees"), value := value/0.06]
wide_dt <- sim2$wide$site
wide_dt <- merge(wide_dt, species_dt, by = "species", all.x = TRUE)
wide_dt[,.(ba = mean(ba)), by = species_name]
library(ggplot2)
ggplot(p_dat[year <= 200],
       aes(x = year, y = value, color = factor(species_name))) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y")


pred <- m1$simulate(env_dt, init_cohort = init_cohorts)
pred_dt <- pred$wide$site

obs_dt_long <- melt(obs_dt, id.vars = c("siteID","year","species"), measure.vars = c("ba","trees","growth","mort","reg"), value.name = "Obs")
# obs_dt_long$SimObs <- "obs"
pred_dt_long <- melt(pred_dt, id.vars = c("siteID","year","species"), measure.vars = c("ba","trees","growth","mort","reg"), value.name = "Sim")
# pred_dt_long$SimObs <- "sim"

comp_dt_long <- merge(obs_dt_long, pred_dt_long, by = c("siteID","year","species","variable"), all = TRUE)

#calculate cors
comp_dt_long[,.(cor = cor(Sim, Obs, use = "complete.obs")), by = .(variable)]

ggplot(comp_dt_long,
       aes(x = Sim, y = Obs, color = factor(species))) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~variable, scales = "free_y")


par1afterSimulate = m1$parameters_r

# fit again ####
m1$fit(
  env        = env_dt,
  data       = unique(obs_dt),
  init_cohort = init_cohorts,
  device     = "cpu",
  epochs     = 500,
  batchsize  = 10L,
  patch_size = 0.06,
  patches    = 4, weights = c(0.1, 10, 1.0, 10.0, 1, 1),
  lr         = 0.01#, loss = c("mse","mse","mse","mse","mse","mse","mse")
)
# sim0 <- m1$simulate(env_dt, init_cohort = NULL)
par2 <- m1$parameters_r
par1$par_theta_recruits_raw
par2$par_theta_recruits_raw

# print difference between each element in par1 and par2
for(i in names(par1)){
  print(i)
  print(par1[[i]] - par2[[i]])
}

# simulate again ####
sim2 <- m1$simulate(env_dt2, init_cohort = NULL)
p_dat <- sim2$long$site[, .(value = mean(value)), by = .(year, species, variable)]
p_dat <- merge(p_dat, species_dt, by = "species", all.x = TRUE)
p_dat[variable %in% c("ba","trees"), value := value/0.06]
wide_dt <- sim2$wide$site
wide_dt <- merge(wide_dt, species_dt, by = "species", all.x = TRUE)
wide_dt[,.(ba = mean(ba)), by = species_name]
library(ggplot2)
ggplot(p_dat[year <= 200],
       aes(x = year, y = value, color = factor(species_name))) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y")



par1$par_theta_recruits_raw
par2$par_theta_recruits_raw


par1$par_theta_recruits_raw
par2$par_theta_recruits_raw
par1$par_competition_unconstrained
par2$par_competition_unconstrained
m1$fit(
  env        = env_dt,
  data       = unique(obs_dt),
  init       = init_cohorts,
  device     = "cpu",
  epochs     = 100,
  batchsize  = 10L,
  patch_size = 0.06,
  patches    = 4, weights = c(0.1, 10, 1.0, 10.0, 1, 1),
  lr         = 0.000000000000000000000001#, loss = c("mse","mse","mse","mse","mse","mse","mse")
)
par3 <- m1$parameters_r
par3$par_theta_recruits_raw
identical(par1,par3)


env_dt2 <- env_dt[year == 1]
for(i in 2:200){
  env_dt2 <- rbind(env_dt2, env_dt[year == 1,.(siteID, year = i, prec, temp)])
}


species_dt <- unique(inputs$obs_dt[,.(species, species_name)])
sim2 <- m1$simulate(env_dt2, init_cohort = NULL)
p_dat <- sim2$long$site[, .(value = mean(value)), by = .(year, species, variable)]
p_dat <- merge(p_dat, species_dt, by = "species", all.x = TRUE)
p_dat[variable %in% c("ba","trees"), value := value/0.06]

wide_dt <- sim2$wide$site
wide_dt <- merge(wide_dt, species_dt, by = "species", all.x = TRUE)
wide_dt[,.(ba = mean(ba)), by = species_name]

library(ggplot2)
ggplot(p_dat[year <= 200],
       aes(x = year, y = value, color = factor(species_name))) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y")

p_dat2 <- sim2$long$site[, .(value = mean(value)), by = .(year, species, variable,siteID)]
p_dat2 <- merge(env_dt, p_dat2, by = c("year","siteID"), all.x = TRUE)
p_dat2[variable %in% c("ba","trees"), value := value/0.06]
p_dat2 <- merge(p_dat2, species_dt, by = "species", all.x = TRUE)

# add prec and value to variable column by making it long

p1 = ggplot(p_dat2,
       aes(x = prec, y = value, color = factor(species_name))) +
  geom_point() +
  geom_smooth()+
  facet_wrap(~variable, scales = "free_y", ncol = 1)+
  theme_minimal()+
  scale_color_discrete(name = "Species")

p2 = ggplot(p_dat2,
       aes(x = temp, y = value, color = factor(species_name))) +
  geom_point() +
  geom_smooth()+
  facet_wrap(~variable, scales = "free_y", ncol = 1)+
  theme_minimal()+
  scale_color_discrete(name = "Species")

#arrange with common legend
library(patchwork)
(p1 | p2) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')


## 12) Save (optional) ----
fwrite(site_dt, "vignettes/FIAsubset_sites.csv")
fwrite(tree_dt, "vignettes/FIAsubset_trees.csv")
fwrite(env_dt,  "vignettes/FIAsubset_env.csv")
fwrite(obs_dt,  "vignettes/FIAsubset_obs.csv")
if (!is.null(init_dt)) fwrite(init_dt, "vignettes/FIAsubset_init.csv")
