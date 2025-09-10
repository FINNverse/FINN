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
Sys.setenv(CUDA_VISIBLE_DEVICES=0)
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
  site_dt <- fread("data/FIA/FIA_sites_raw.csv")
}

hist(tree_dt[,.(Ntrees = uniqueN(treeName)), by = .(siteName, patchName)]$Ntrees)

# convenience lag for obs building
setorder(tree_dt, treeName, year)
tree_dt[, status_before := shift(status, 1, type = "lag"), by = treeName]

## 5) Derive environmental covariates from site_dt  [site_dt ➜ get env ➜ env_dt] ----
# NOTE: env is computed per (siteID, year) as in the diagram.
env_dt <- copy(site_dt)   # keep siteID, year, x, y

site_pts <- vect(env_dt, geom = c("x", "y"), crs = "EPSG:4326")

# Example: WorldClim BIO1 (temp) & BIO12 (prec)
temp <- terra::rast("../worldclim/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")
tempmax <- terra::rast("../worldclim/wc2.1_30s_bio/wc2.1_30s_bio_5.tif")
tempmin <- terra::rast("../worldclim/wc2.1_30s_bio/wc2.1_30s_bio_6.tif")
prec <- terra::rast("../worldclim/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")
precseas <- terra::rast("../worldclim/wc2.1_30s_bio/wc2.1_30s_bio_15.tif")
precwarmq <- terra::rast("../worldclim/wc2.1_30s_bio/wc2.1_30s_bio_18.tif")

v_proj <- project(site_pts, crs(temp))
# env_dt[, temp := extract(temp, v_proj, method = "near")[, 2] ]
# env_dt[, prec := extract(prec, v_proj, method = "near")[, 2] ]
env_dt[, ":="(
  temp = extract(temp, v_proj, method = "near")[, 2],
  tempmax = extract(tempmax, v_proj, method = "near")[, 2],
  tempmin = extract(tempmin, v_proj, method = "near")[, 2],
  prec = extract(prec, v_proj, method = "near")[, 2],
  precseas = extract(precseas, v_proj, method = "near")[, 2],
  precwarmq = extract(precwarmq, v_proj, method = "near")[, 2]
),]

env_dt <- env_dt[, -c("x","y")]

## 6) Create observation table  [tree_dt ➜ create obs_dt] ----
# flags required by makeObsData(): reg/mort/living per (tree, year)
td <- copy(tree_dt)
td[status == "new",                          reg := TRUE]
td[status == "dead" & status_before == "alive", mort := TRUE]
td[status %chin% c("alive", "new"),          living := TRUE]

td <- td[!(siteName %in% unique(td[mort_cause %in% c("land use", "fire")]$siteName))]

# Aggregate to FINN obs structure (ba/trees/dbh/growth/mort/reg per site/patch/species/year)
# sampled_siteNames <- sample(unique(td$siteName), 500)
add_site_completeness_flags(env_dt, td)

td[,.(Npatches = uniqueN(patchName)),.(siteName)]$Npatches |> table()
td[,.(Nyears = uniqueN(year)),.(siteName)]$Nyears |> table()

obs_list <- makeObsData(
  tree_dt = td[complete == T],
  plotsize = 0.06,
  aggregate_by_site = F,
  NspeciesQuantile = 0.99,
  #Nspecies = 10,
  Npatches = 4, minNyears = 3)

# select sites from obs_list in sites_pts and plot
selected_sites_pts <- site_pts[site_pts$siteName %in% unique(obs_list$obs_dt$siteName),]
plot(site_pts)
# add selected_sites_pts to plot
points(selected_sites_pts, col = "red", pch = 19)

# obs_list$obs_dt[,.(Nyear = uniqueN(year)), by = siteName]
set.seed(2)
obs_dt = obs_list$obs_dt
# obs_dt[,.(Nyear = uniqueN(year)), by = siteName]
tree_dt = obs_list$tree_dt[siteName %in% unique(obs_dt$siteName)]
# tree_dt[,.(Nyear = uniqueN(year)), by = siteName]

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
env_dt <- inputs$env_dt[,-c("complete","OrigYear","siteName")]
#scale env_dt
env_dt_names <- c("siteID","year", names(env_dt)[!(names(env_dt) %in% c("siteID","year"))])
env_dt <- env_dt[, ..env_dt_names]

scale_vals <- list()
par(mfrow = c(3,3))
for(i in names(env_dt)[!(names(env_dt) %in% c("siteID","year"))]){
  hist(env_dt[[i]], main = i)
  scale_vals[[i]]$mean <- mean(env_dt[[i]])
  scale_vals[[i]]$sd <- sd(env_dt[[i]])
  env_dt[[i]] <- as.numeric(scale(env_dt[[i]]))
}

init_trees <- inputs$tree_dt[year == 0]
init_trees[living == F, trees := 0,]
init_cohorts <- makeInitCohorts(init_trees, Nspecies = max(obs_dt$species))

m1 <- finn(
  N_species            = uniqueN(obs_dt$species),
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0,     FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createProcess(~., FINN::growth,   optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~., FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~., FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)

# m1 <- torch::torch_load("data/fia_v5_process_finn.pt")

# first fit ####
m1$fit(
  env        = env_dt,
  data       = unique(obs_dt),
  init_cohort = init_cohorts,
  device     = "gpu",
  epochs     = 10000,
  batchsize  = 200L,
  patch_size = 0.06,
  patches    = 4, weights = c(0.1, 10, 1.0, 10.0, 1, 1),
  lr         = 0.01#, loss = c("mse","mse","mse","mse","mse","mse","mse")
)

torch::torch_save(m1, "data/fia_v7_process_finn.pt")

m1$par_theta_recruits
# simulate ####
sampled_sites <- sample(unique(env_dt$siteID),100)
env_dt2 <- env_dt[siteID %in% sampled_sites & year == 1]
for(i in 1:1500){
  env_dt_temp <- env_dt[siteID %in% sampled_sites & year == 1]
  env_dt_temp$year = i
  env_dt2 <- rbind(env_dt2, env_dt_temp)
}

species_dt <- unique(inputs$obs_dt[,.(species, species_name)])
sim2 <- m1$simulate(env_dt2, init_cohort = NULL, device = "gpu", patches = 100, patch_size = 0.06)
p_dat <- sim2$long$site[, .(value = mean(value)), by = .(year, species, variable)]
p_dat <- merge(p_dat, species_dt, by = "species", all.x = TRUE)
p_dat[variable %in% c("ba","trees"), value := value/0.06]

wide_dt <- sim2$wide$site
wide_dt <- merge(wide_dt, species_dt, by = "species", all.x = TRUE)
wide_dt[,.(ba = mean(ba)), by = species_name][order(ba)]
wide_dt[,.(mort = mean(mort)), by = species_name][order(mort)]

m1$par_theta_recruits
library(ggplot2)
ggplot(p_dat[year <= 1500],
       aes(x = year, y = value, color = factor(species_name))) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y")

p_dat2 <- sim2$long$site[, .(value = mean(value)), by = .(year, species, variable,siteID)]
# p_dat2 <- merge(env_dt[,-"year"], p_dat2, by = c("siteID"), all.x = TRUE)
p_dat2[variable %in% c("ba","trees"), value := value/0.06]
p_dat2 <- merge(p_dat2, species_dt, by = "species", all.x = TRUE)

site_dt2 <- fread("data/FIA/FIA_sites_raw.csv")

site_dt2 <- merge(site_dt2, unique(inputs$siteID_dt[,.(siteName, siteID)]), by = "siteName")

p_dat3 <- merge(p_dat2, unique(site_dt2[,.(siteID,x,y)]), by = "siteID")

p_dat4 <- p_dat3[variable == "ba" & year == 1500, .(ba_share = mean(value, na.rm = T)), by = .(lon = cut(x, breaks= 50), species_name)]
ggplot(p_dat4, aes(x = lon, y = ba_share, fill = species_name))+
  geom_bar(stat = "identity")


# add prec and value to variable column by making it long


p_dat2_env <- sim2$long$site[, .(value = mean(value)), by = .(year, species, variable,siteID)]
p_dat2_env <- merge(unique(env_dt[,-"year"]), p_dat2_env, by = c("siteID"), all.x = TRUE)
p_dat2_env[variable %in% c("ba","trees"), value := value/0.06]
p_dat2_env <- merge(p_dat2_env, species_dt, by = "species", all.x = TRUE)

p1 = ggplot(p_dat2_env[variable == "reg"],
            aes(x = prec, y = value)) +
  geom_point() +
  geom_smooth()+
  facet_wrap(~species_name, scales = "free_y", ncol = 2)+
  theme_minimal()+
  scale_color_discrete(name = "Species")

p2 = ggplot(p_dat2_env[variable == "reg"],
            aes(x = temp, y = value)) +
  geom_point() +
  geom_smooth()+
  facet_wrap(~species_name, scales = "free_y", ncol = 2)+
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

growth = m1$parameters_r$nn_growth.0.weight

growth <- cbind(species_dt[,-"species"], growth)
growth_dt <- data.table::melt(growth, id.vars = "species_name")
ggplot(growth_dt, aes(x = value, y = species_name,))+
  geom_bar(stat = "identity")+
  facet_wrap(~variable)

mort = m1$parameters_r$nn_mortality.0.weight
mort <- cbind(species_dt[,-"species"], mort)
mort_dt <- data.table::melt(mort, id.vars = "species_name")
ggplot(mort_dt, aes(x = value, y = species_name,))+
  geom_bar(stat = "identity")+
  facet_wrap(~variable)

reg = m1$parameters_r$nn_regeneration.0.weight
reg <- cbind(species_dt[,-"species"], reg)
reg_dt <- data.table::melt(reg, id.vars = "species_name")
ggplot(reg_dt, aes(x = value, y = species_name,))+
  geom_bar(stat = "identity")+
  facet_wrap(~variable)


comp1_dt <- rbindlist(list(
  data.table(growth_dt, process = "growth"),
  data.table(mort_dt, process = "mort"),
  data.table(reg_dt, process = "reg")
))
fwrite(comp1_dt, "data/fia_v7_process_finn_envpars.csv")

print.data.frame(comp1_dt,max = 500)

comp1_dt[, variable := factor(
  variable,
  levels = paste0("V", 1:(length(scale_vals)+1)),
  labels = c("intercept", names(scale_vals))),]


ggplot(comp1_dt, aes(x = value, y = species_name, fill = process))+
  geom_bar(stat = "identity", position = position_dodge2())+
  facet_wrap(~variable,nrow = 1)

disp_dt <- data.table(
  disp_recr = torch::as_array(m1$par_theta_recruits),
  species_name = species_dt$species_name
)
ggplot(disp_dt, aes(x = disp_recr, y = forcats::fct_reorder(species_name, disp_recr)))+
  geom_bar(stat = "identity")

m1$par_theta_recruit

par_mort <- m1$parameters_r$par_mortality_unconstrained
colnames(par_mort) <- c("mortLight", "mortSize", "mortGrowth")

par_growth <- m1$parameters_r$par_mortality_unconstrained
colnames(par_mort) <- c("mortLight", "mortSize", "mortGrowth")



