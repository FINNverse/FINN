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
library(ggplot2)
Sys.setenv(CUDA_VISIBLE_DEVICES=0)
source("R/createInputs.R")   # supplies makeObsData(), resolveSiteIDs(), createInitCohorts() if present

if(!file.exists("data/FIA/FIA_trees_raw.csv")){
# if(T){

  library(rFIA)
  dir.path <- "data/FIA"       # <— adjust to your layout

  ## 1) Download / Load FIA (state = OR) ----
  # options(timeout = 999999)
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
temp <- terra::rast("~/Downloads/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")
tempmax <- terra::rast("~/Downloads/wc2.1_30s_bio/wc2.1_30s_bio_5.tif")
tempmin <- terra::rast("~/Downloads/wc2.1_30s_bio/wc2.1_30s_bio_6.tif")
prec <- terra::rast("~/Downloads/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")
precseas <- terra::rast("~/Downloads/wc2.1_30s_bio/wc2.1_30s_bio_15.tif")
precwarmq <- terra::rast("~/Downloads/wc2.1_30s_bio/wc2.1_30s_bio_18.tif")

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
td[status == "dead" & status_before %in% c("alive","new"), mort := TRUE]
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
  NspeciesQuantile = 0.98,
  fix_period_length = 10,
  dbh_growth_thresh = c(-10,50),
  #Nspecies = 5,
  Npatches = 4, minNyears = 3)

trees_dt <- obs_list$tree_dt

hist(trees_dt$dbh_growth)
hist(trees_dt$rel_growth)

# select sites from obs_list in sites_pts and plot
selected_sites_pts <- site_pts[site_pts$siteName %in% unique(obs_list$obs_dt$siteName),]
plot(site_pts)
# add selected_sites_pts to plot
points(selected_sites_pts, col = "red", pch = 19)

# obs_list$obs_dt[,.(Nyear = uniqueN(year)), by = siteName]
set.seed(2)
# uniqueN(obs_list$obs_dt$siteName)
# uniqueN(excluded_sites)
# excluded_sites <- obs_list$obs_dt[,.(ba = sum(ba,na.rm = T)), by = .(siteName,year)][ba == 0]$siteName
# dev.off()
obs_dt = obs_list$obs_dt
# obs_dt = obs_list$obs_dt[!(siteName %in% excluded_sites)]
# obs_dt[,.(Nyear = uniqueN(year)), by = siteName]
tree_dt = obs_list$tree_dt[siteName %in% unique(obs_dt$siteName)]
# tree_dt[,.(Nyear = uniqueN(year)), by = siteName]

## 7) Select sites (species & patch consistency)  [select sites] ----
# Example selection: keep top N species and patches with full time coverage

## 9) Resolve IDs & finalize model inputs  [resolve IDs ➜ {env_dt, obs_dt, init_dt}] ----
# This step aligns site/patch/species IDs across env_dt, obs_dt, init_dt and optionally
# converts years→periods and fixes the number of patches, as in your original call.

# obs_dt[, period := as.integer(as.factor(year)), by = .(siteName)]

inputs <- resolveSiteIDs(
  tree_dt            = as.data.table(tree_dt),
  obs_dt             = as.data.table(obs_dt),
  env_dt             = as.data.table(env_dt),
  createInitCohorts = F
)

tree_dt <- inputs$tree_dt[year != 0]
tree_dt_complete <- inputs$tree_dt
obs_dt <- inputs$obs_dt[year != 0]
env_dt <- inputs$env_dt[,-c("complete","OrigYear","siteName")]
summary(obs_dt)
#scale env_dt
env_dt_names <- c("siteID","year", names(env_dt)[!(names(env_dt) %in% c("siteID","year"))])
env_unscaled_dt <- copy(env_dt[, ..env_dt_names])
env_dt <- env_dt[, ..env_dt_names]

env_scales_dt <- data.table()
par(mfrow = c(3,3))
for(i in names(env_dt)[!(names(env_dt) %in% c("siteID","year"))]){
  hist(env_dt[[i]], main = i)
  env_scales_dt <- rbind(
    env_scales_dt,
    data.table(variable = i,mean = mean(env_dt[[i]]),sd = sd(env_dt[[i]]))
    )
  env_dt[[i]] <- as.numeric(scale(env_dt[[i]]))
}


init_trees <- inputs$tree_dt[year == 0]
init_trees[living == F, trees := 0,]

if(!dir.exists("vignettes/fit-fia-data")){
  dir.create("vignettes/fit-fia-data")
}

summary(obs_dt)
fwrite(env_scales_dt, "vignettes/fit-fia-data/env_scales_dt.csv")
fwrite(env_unscaled_dt, "vignettes/fit-fia-data/env_unscaled_dt.csv")
fwrite(env_dt, "vignettes/fit-fia-data/env_dt.csv")
fwrite(obs_dt, "vignettes/fit-fia-data/obs_dt.csv")
fwrite(inputs$obs_dt, "vignettes/fit-fia-data/full_obs_dt.csv")
fwrite(inputs$tree_dt, "vignettes/fit-fia-data/full_tree_dt.csv")
fwrite(init_trees, "vignettes/fit-fia-data/init_trees.csv")
species_dt <- unique(obs_dt[,.(species,species_name)])
fwrite(species_dt, "vignettes/fit-fia-data/species_dt.csv")

init_cohorts <- makeInitCohorts(init_trees, Nspecies = max(obs_dt$species))

m1 <- finn(
  N_species            = uniqueN(obs_dt$species),
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0,     FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createProcess(~., FINN::growth,   optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~., FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~., FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)

# first fit ####
# m1$fit(
#   env        = env_dt,
#   data       = obs_dt,
#   init_cohort = init_cohorts,
#   device     = "cpu",
#   epochs     = 5000,
#   batchsize  = 500L,
#   patch_size = 0.06,
#   patches    = 4, weights = c(0.1, 10, 1.0, 10.0, 1, 1),
#   lr         = 0.01#, loss = c("mse","mse","mse","mse","mse","mse","mse")
# )

# torch::torch_save(m1, "data/fia_v15_10k_process_finn.pt")
m1 <- torch::torch_load("data/fia_v15_10k_process_finn.pt")


## 1. numeric coordinates from the figure
lon_w <- c(124.00, 123.75, 123.50, 123.25, 123.00,
           122.75, 122.50, 122.37, 122.25, 122.12,
           122.00, 121.96, 121.92, 121.87, 121.83,
           121.79, 121.75, 121.71, 121.67, 121.62,
           121.58, 121.54, 121.50, 121.46, 121.37,
           121.25, 121.00)

lat_const <- 44.13  # 44.13°N from the caption

## 2. build env_dt with lon = negative west, lat = constant
env_dt_BugSol2000 <- data.table(
  x = -lon_w,                 # longitude in decimal deg East (west is negative)
  y = rep(lat_const, length(lon_w))  # latitude
)

## 3. create SpatVector in EPSG:4326
site_pts_BugSol2000 <- vect(env_dt_BugSol2000, geom = c("x", "y"), crs = "EPSG:4326")

## 5. project points to raster CRS (WorldClim is already EPSG:4326 usually,
##    but this is still valid and keeps your workflow unchanged)
v_proj_BugSol2000 <- project(site_pts_BugSol2000, crs(temp))

# install.packages("mapview")
# library(mapview)
# mapview(v_proj_BugSol2000)

# env_dt[, temp := extract(temp, v_proj, method = "near")[, 2] ]
# env_dt[, prec := extract(prec, v_proj, method = "near")[, 2] ]
env_dt_BugSol2000[, ":="(
  temp = extract(temp, v_proj_BugSol2000, method = "near")[, 2],
  tempmax = extract(tempmax, v_proj_BugSol2000, method = "near")[, 2],
  tempmin = extract(tempmin, v_proj_BugSol2000, method = "near")[, 2],
  prec = extract(prec, v_proj_BugSol2000, method = "near")[, 2],
  precseas = extract(precseas, v_proj_BugSol2000, method = "near")[, 2],
  precwarmq = extract(precwarmq, v_proj_BugSol2000, method = "near")[, 2]
),]


env_dt_BugSol2000scaled <- copy(env_dt_BugSol2000)

for(i in names(env_dt_BugSol2000)[!(names(env_dt_BugSol2000) %in% c("x","y"))]){
  scale_mean <- env_scales_dt[variable == i]$mean
  scale_sd <- env_scales_dt[variable == i]$sd
  env_dt_BugSol2000scaled[[i]] <- (env_dt_BugSol2000[[i]] - scale_mean) / scale_sd
}

# m1 <- torch::torch_load("data/fia_v15_10k_process_finn.pt")

env_dt_BugSol2000scaled[, siteID := 1:.N]

Nyears = 1000
sim_env_dt_BugSol2000scaled <- data.table()
for(i in 1:Nyears){
  sim_env_dt_BugSol2000scaled <- rbind(
    sim_env_dt_BugSol2000scaled,
    data.table(env_dt_BugSol2000scaled[,-c("x","y")], year = i)
    )
}

dist_dt <- sim_env_dt_BugSol2000scaled[,.(siteID,year)]

dist_dt$intensity <- rbinom(nrow(dist_dt),1,0.01)*runif(nrow(dist_dt),0.1,0.2)


dist_dt[, intensity := rbinom(.N,1,0.01)*runif(.N,0.1,0.2), by = siteID]

hist(dist_dt$intensity)
table(dist_dt$intensity)

sim1 <- m1$simulate(sim_env_dt_BugSol2000scaled,disturbance = dist_dt, init_cohort = NULL, device = "cpu", patches = 100, patch_size = 0.06)

sim_long_out <- copy(sim1$long$site)
sim_wide_out <- copy(sim1$wide$site)


# sim_long_out <- merge(sim_long_out, species_dt, by = "species", all.x = TRUE)
# sim_long_out[variable %in% c("ba","trees"), value := value/0.06]
# pdf("succesions.pdf",width = 15,height = 10)
# for(i in unique(sim_long_out$siteID)){
#   p =   ggplot(sim_long_out[species %in% selected_species & siteID == i],
#            aes(x = year, y = value, color = species_name)) +
#       geom_line() +
#       facet_wrap(~variable, scales = "free")+
#       ggtitle(paste0(i))
#   print(p)
# }
# dev.off()
# dev.off()
par(mfrow = c(2,1), mar = c(5, 4, 4, 1))
plot(temp~x,data = env_dt_BugSol2000, type = "b", las = 1, ylab = "Temperature (°C)", xlab = "Longitude (°W)")
plot(prec/10~x,data = env_dt_BugSol2000, type = "b", las = 1, ylab = "Precipitation (cm)", xlab = "Longitude (°W)")

p_dat <- sim_long_out[, .(value = mean(value)), by = .(year, species, variable)]
# p_dat <- sim_long_out
p_dat <- merge(p_dat, species_dt, by = "species", all.x = TRUE)
p_dat[variable %in% c("ba","trees"), value := value/0.06]


# sim_wide_out[year == 1000, .(ba = sum(ba)), .(species)][order(ba)]$species
# sim_wide_out[, ba_rank := as.integer(factor(ba,ordered = T)), .(siteID,year)]

ggplot(p_dat,
       aes(x = year, y = value, color = species_name)) +
  geom_line() +
  facet_wrap(~variable, scales = "free")

selected_species = unique(p_dat[variable == "ba" & value > 3]$species_name)
sim_long_out <- merge(sim_long_out, species_dt, by = "species", all.x = TRUE)
sim_long_out[variable %in% c("ba","trees"), value := value/0.06]
pdf("succesions.pdf",width = 15,height = 10)
for(i in unique(sim_long_out$siteID)){
  p =   ggplot(sim_long_out[species_name %in% selected_species & siteID == i],
           aes(x = year, y = value, color = species_name)) +
      geom_line() +
      facet_wrap(~variable, scales = "free")+
      ggtitle(paste0(i))
  print(p)
}
dev.off()


p_dat2 <- sim1$long$site[, .(value = mean(value)), by = .(year, species, variable,siteID)]
# p_dat2 <- merge(env_dt[,-"year"], p_dat2, by = c("siteID"), all.x = TRUE)
p_dat2[variable %in% c("ba","trees"), value := value/0.06]
p_dat2 <- merge(p_dat2, species_dt, by = "species", all.x = TRUE)

Bug_species <- c(
  "Abies amabilis",
  "Abies grandis",
  "Chamaecyparis nootkatensis",
  "Picea sitchensis",
  "Pinus ponderosa",
  "Pseudotsuga menziesii",
  "Thuja plicata",
  "Tsuga heterophylla",
  "Tsuga mertensiana"
  )
p_dat2[!(species_name %in% Bug_species), species_name := "Other species",]

p_dat4 <- p_dat2[variable == "ba" & year > 999, .(ba = sum(value, na.rm = T)/.N), by = .(siteID, species_name)]

p_dat5 <- merge(unique(env_dt_BugSol2000scaled[, .(siteID, x, y)]), p_dat4, by = "siteID", all.x = TRUE)

p_dat5$lon <- format(abs(p_dat5$x), 2)
p_dat5$lon <- factor(p_dat5$lon, levels = sort(unique(p_dat5$lon), decreasing = T))
p_dat5[,species_name := factor(species_name, levels = c(Bug_species,"Other species")),]

ggplot(p_dat5, aes(x = lon, y = ba, fill = species_name))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  xlab("Longitude (°W)")


ggplot(p_dat5[species_name != "Other species"], aes(x = lon, y = ba, fill = as.character(species_name)))+
  geom_bar(stat = "identity")

library(ggplot2)
library(patchwork)


unique(p_dat2$variable)
p_dat4 <- p_dat2[year == 1000, .(
  ba = sum(value[variable == "ba"], na.rm = T)/.N,
  trees = sum(value[variable == "trees"], na.rm = T)/.N
  ), by = .(siteID, species_name)]

p_dat5 <- merge(unique(env_dt_BugSol2000scaled[, .(siteID, x, y)]), p_dat4, by = "siteID", all.x = TRUE)

p_dat5$lon <- format(abs(p_dat5$x), 2)
p_dat5$lon <- factor(p_dat5$lon, levels = sort(unique(p_dat5$lon), decreasing = T))
p_dat5[,species_name := factor(species_name, levels = c(Bug_species,"Other species")),]


# Temperaturplot
p_temp <- ggplot(env_dt_BugSol2000, aes(x = x, y = temp)) +
  geom_line(colour = "firebrick", linewidth = 1) +
  geom_point(colour = "firebrick") +
  theme_bw() +
  labs(x = "Longitude (°W)", y = "Temperature (°C)") +
  theme(legend.position = "none")

# Niederschlagsplot
p_prec <- ggplot(env_dt_BugSol2000, aes(x = x, y = prec / 10)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_point(colour = "steelblue") +
  theme_bw() +
  labs(x = "Longitude (°W)", y = "Precipitation (cm)") +
  theme(legend.position = "none")

# Balkendiagramm der Arten
hist(p_dat5$trees)
p_species <- ggplot(p_dat5, aes(x = as.integer(lon), y = ba, fill = species_name)) +
  geom_bar(stat = "identity") +
  # geom_line(linewidth = 2) +
  theme_bw() +
  labs(x = "Longitude (°W)", y = "Basal area") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank())

# Kombination: Temperatur oben, Niederschlag darunter, Balken unten
combined_plot <- p_temp / p_prec / p_species +
  plot_layout(heights = c(1, 1, 1.2))

combined_plot


bug_pars <- fread("data/SpecParsPNW.DAT", sep = "\t", skip = 6, encoding = "Latin-1")
chr <- names(bug_pars)[vapply(bug_pars, is.character, logical(1))]
bug_pars[, (chr) := lapply(.SD, function(x) chartr("\u00CA", " ", x)), .SDcols = chr]
# bug_pars[,species_name := gsub("(g)","",,Species),]
bug_pars[Species != "Pinus contorta(l)"]
bug_pars[Species == "Pinus contorta(c)", Species == "Pinus contorta" ,]
bug_pars[Species != "Pseudotsuga menziesii(g)"]
bug_pars[Species == "Pseudotsuga menziesii(m)", Species == "Pseudotsuga menziesii" ,]

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

env_names = names(env_dt)[-c(1:2)]
comp1_dt[, variable := factor(
  variable,
  levels = paste0("V", 1:(length(env_names)+1)),
  labels = c("intercept", env_names)),]


bug_pars[,species_name := Species,]
bug_pars[,WiTN := as.numeric(WiTN),]



disp_dt <- data.table(
  disp_recr = torch::as_array(m1$par_theta_recruits),
  species_name = species_dt$species_name
)
ggplot(disp_dt, aes(x = disp_recr, y = forcats::fct_reorder(species_name, disp_recr)))+
  geom_bar(stat = "identity")


par_mort <- data.table(m1$parameters_r$par_mortality_unconstrained)
par_mort$species_name <- species_dt$species_name
colnames(par_mort) <- c("mortLight", "mortSize", "mortGrowth","species")

par_growth <- data.table(m1$parameters_r$par_growth_unconstrained)
par_growth$species_name <- species_dt$species_name
colnames(par_growth) <- c("growthLight", "growthSize","species")

par_reg <- data.table(m1$parameters_r$par_regeneration_unconstrained)
par_reg$species_name <- species_dt$species_name
colnames(par_reg) <- c("regLight","species")

par_mort <- melt(data.table(par_mort), measure.vars = colnames(par_mort)[-length(colnames(par_mort))])
par_growth <- melt(data.table(par_growth), measure.vars = colnames(par_growth)[-length(colnames(par_growth))])
par_reg <- melt(data.table(par_reg), measure.vars = colnames(par_reg)[-length(colnames(par_reg))])

process_par_dt <- rbindlist(list(
  data.table(par_mort, process = "mort"),
  data.table(par_growth, process = "growth"),
  data.table(par_reg, process = "reg"),
  data.table(species_name = disp_dt$species_name, variable = "regDispersion", value = disp_dt$disp_recr, process = "reg")
))

names(process_par_dt)[1] <- "species_name"
comp1_dt <- rbind(process_par_dt,comp1_dt)

comp_p_dt <- merge(bug_pars,comp1_dt, by = "species_name")
c("Dm" = "evaporative demand")
forclimpars = names(bug_pars[,-c("Species","species_name")])

names(bug_pars[,-c("Species","species_name")])
param_defs <- c(
  sType = "species type",
  Dm    = "max diameter",
  Hm    = "max height",
  Am    = "max age",
  G     = "growth rate",
  DDMin = "min degree days",
  WiTN  = "winter Tmin",
  WiTX  = "winter Tmax",
  DrT   = "drought tolerance",
  kFlT  = "foliage turnover",
  NTol  = "nutrient tolerance",
  Brow  = "browsing tolerance",
  Ly    = "light requirements young trees",
  La    = "light requirements old trees",
  LQ    = "litter quality",
  Imm   = "immaturity age"
)


pdf("comp_pars.pdf",width = 10, height = 10)
p = ggplot(comp_p_dt, aes(y = value, x = sType, color = process))+
  geom_boxplot()+
  facet_wrap(~variable)+
  ggtitle(param_defs[1])
print(p)
for(i in 2:length(param_defs)){
  p = ggplot(comp_p_dt, aes_string(y = "value", x = names(param_defs)[i], color = "process"))+
    geom_point()+
    geom_smooth(method = "lm")+
    facet_wrap(~variable, scales = "free")+
    ggtitle(param_defs[i])
  print(p)
}
dev.off()
cors_dt <- data.table()
summary(comp_p_dt)
cor(comp_p_dt[variable == "mortLight"]$Dm,comp_p_dt[variable == "mortLight"]$value)
cor2 <- function(x,y,...){
  if(length(x) > 0) cor(x,y,...) else NA
}

cors_temp <- comp_p_dt[,.(
  Dm = cor(Dm, value,method = "spearman", use = "complete.obs"),
  Hm = cor(Hm, value,method = "spearman", use = "complete.obs"),
  Am = cor(Am, value,method = "spearman", use = "complete.obs"),
  G = cor(G, value,method = "spearman", use = "complete.obs"),
  DDMin = cor(DDMin, value,method = "spearman", use = "complete.obs"),
  WiTN = cor(WiTN, value,method = "spearman", use = "complete.obs"),
  WiTX = cor(WiTX, value,method = "spearman", use = "complete.obs"),
  DrT = cor(DrT, value,method = "spearman", use = "complete.obs"),
  kFlT = cor(kFlT, value,method = "spearman", use = "complete.obs"),
  NTol = cor(NTol, value,method = "spearman", use = "complete.obs"),
  Brow = cor(Brow, value,method = "spearman", use = "complete.obs"),
  Ly = cor(Ly, value,method = "spearman", use = "complete.obs"),
  La = cor(La, value,method = "spearman", use = "complete.obs"),
  LQ = cor(LQ, value,method = "spearman", use = "complete.obs")
  ), by = .(variable,process)]
p_cors <- melt(cors_temp, id.vars = c("variable","process"))
ggplot(p_cors, aes(y = variable, x = variable.1, fill = value))+
  geom_tile()+
  facet_wrap(~process, scales = "free", ncol = 2)+
  scale_fill_gradient2(limits = c(-1,1))+
  geom_text(aes(label = round(value, 2)))

p = ggplot(comp_p_dt[Brow < 2], aes_string(y = "value", x = "Brow", color = "process"))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_grid(process~variable)+
  ggtitle(param_defs[i])
print(p)

ggplot(comp_p_dt, aes(y = value, x = Hm, color = process))+
  geom_point()+
  facet_wrap(process~variable)



