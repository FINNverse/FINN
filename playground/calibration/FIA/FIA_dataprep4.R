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
  #Nspecies = 5,
  Npatches = 4, minNyears = 3)

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
  device     = "cpu",
  epochs     = 2000,
  batchsize  = 200L,
  patch_size = 0.06,
  patches    = 4, weights = c(0.1, 10, 1.0, 10.0, 1, 1),
  lr         = 0.01#, loss = c("mse","mse","mse","mse","mse","mse","mse")
)
#
torch::torch_save(m1, paste0("data/fia_v10_process_finn.pt"))
# m1 <- torch::torch_load("data/fia_v9(proposal_noOther)_process_finn.pt")

# simulate for fit evaluation
sim2 <- m1$simulate(env_dt, init_cohort = init_cohorts, device = "cpu", patches = 4, patch_size = 0.06)
species_dt <- unique(inputs$obs_dt[,.(species, species_name)])
sim_dt <- sim2$wide$site


# create comparison table to check fit
comp_dt <- merge(
  sim_dt[,.(
    ba = sum(ba, na.rm = T)/0.06,
    trees = sum(trees, na.rm = T)/0.06,
    growth = mean(growth, na.rm = T),
    mort = mean(mort, na.rm = T),
    r_mean_ha = mean(r_mean_ha, na.rm = T)
  ), by = .(siteID, year, species)],
  obs_dt[,.(
    ba = sum(ba, na.rm = T)/0.06,
    trees = sum(trees, na.rm = T)/0.06,
    growth = mean(growth, na.rm = T),
    mort = mean(mort, na.rm = T),
    r_mean_ha = mean(reg, na.rm = T)/0.06
  ), by = .(siteID, year, species)],
  by = c("siteID", "year","species"), suffixes = c(".sim", ".obs"))


safe_cor <- function(x, y, method = "pearson", use = "complete.obs") {
  # both vectors entirely zeros → return 1
  if (all(x == 0, na.rm = TRUE) && all(y == 0, na.rm = TRUE)) {
    return(1)
  }
  # if one has all zeros and the other doesn't → correlation undefined → return NA
  if ((all(x == 0, na.rm = TRUE) && !all(y == 0, na.rm = TRUE)) ||
      (all(y == 0, na.rm = TRUE) && !all(x == 0, na.rm = TRUE))) {
    return(NA_real_)
  }
  # otherwise, return normal correlation
  suppressWarnings(cor(x, y, method = method, use = use))
}

comp_dt <- merge(species_dt, comp_dt, by = c("species"))

cors_dt <-
comp_dt[, .(
  r_ba = safe_cor(ba.sim, ba.obs, method = "spearman"),
  r_trees = safe_cor(trees.sim, trees.obs, method = "spearman"),
  r_growth = safe_cor(growth.sim, growth.obs, method = "spearman"),
  r_mort = safe_cor(mort.sim, mort.obs, method = "spearman"),
  r_r_mean_ha = safe_cor(r_mean_ha.sim, r_mean_ha.obs, method = "spearman")
), by = species_name]


cors_long <- melt(cors_dt, id.vars = "species_name", variable.name = "metric", value.name = "correlation")

library(ggplot2)
ggplot(cors_long, aes(x = correlation, y = species_name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Spearman Correlations by Species",
    x = "Correlation (ρ)",
    y = "Species",
    fill = "Variable"
  )+
  facet_wrap(~metric, nrow = 1)

comp_dt_long <- melt(comp_dt[,-"species"], id.vars = c("species_name", "siteID", "year"), variable.name = "metric")

comp_dt_long[, variable := tstrsplit(metric,"\\.")[[1]],]
comp_dt_long[, simobs := tstrsplit(metric,"\\.")[[2]],]

comp_dt_long <- dcast(comp_dt_long, species_name+siteID+year+variable~simobs)

axis_range = range(c(comp_dt_long[variable == "ba"]$sim,comp_dt_long[variable == "ba"]$obs))
ggplot(comp_dt_long[variable == "ba"], aes(x = sim, y = obs))+
  geom_point()+
  geom_smooth(method = "lm", se = F)+
  facet_wrap(species_name~variable)+
  geom_abline(slope = 1)+
  coord_equal()

axis_range = range(c(comp_dt_long[variable == "trees"]$sim,comp_dt_long[variable == "trees"]$obs))
ggplot(comp_dt_long[variable == "trees"], aes(x = sim, y = obs))+
  geom_point()+
  geom_smooth(method = "lm", se = F)+
  facet_wrap(species_name~variable)+
  geom_abline(slope = 1)+
  coord_equal(
    xlim = axis_range,
    ylim = axis_range
    )

# axis_range = range(c(comp_dt_long[variable == "growth"]$sim,comp_dt_long[variable == "growth"]$obs))
ggplot(comp_dt_long[variable == "growth"], aes(x = sim, y = obs))+
  geom_point()+
  geom_smooth(method = "lm", se = F)+
  facet_wrap(species_name~variable)+
  geom_abline(slope = 1)+
  coord_equal(xlim = c(0,1), ylim = c(0,1))# viele ausreißer

ggplot(comp_dt_long[variable == "mort"], aes(x = sim, y = obs))+
  geom_point()+
  geom_smooth(method = "lm", se = F)+
  facet_wrap(species_name~variable)+
  geom_abline(slope = 1)
  # coord_cartesian(xlim = c(0,1), ylim = c(0,1))


axis_range = range(c(comp_dt_long[variable == "r_mean_ha"]$sim,comp_dt_long[variable == "r_mean_ha"]$obs))
ggplot(comp_dt_long[variable == "r_mean_ha"], aes(x = sim, y = obs))+
  geom_point()+
  geom_smooth(method = "lm", se = F)+
  facet_wrap(species_name~variable)+
  geom_abline(slope = 1)+
  coord_equal(
    xlim = axis_range,
    ylim = axis_range
  )

env_dt2 <- data.table(
  siteID = 1,
  year = 1:2000,
  temp = 0,
  tempmax = 0,
  tempmin = 0,
  prec = 0,
  precseas = 0,
  precwarmq = 0
)

# par(mfrow = c(2,3))
# for(i in names(env_dt)[!(names(env_dt) %in% c("siteID","year"))]){
#   hist(env_dt[[i]], main = i)
# }
# summary(env_dt)
# env_dt2[year == 700, ":="(
#   temp = 2,
#   tempmax = 2.7,
#   tempmin = tempmin,
#   prec = -1.2,
#   precseas = precseas+0.5,
#   precwarmq = -1.3
# ),]

species_dt <- unique(inputs$obs_dt[,.(species, species_name)])

# Grab the leaf parameter
w <- m1$nn_regeneration[[1]]$weight  # this is the 1st layer's weight

# Work on CPU/R array, change it, then copy back on the right device/dtype
parmat <- torch::as_array(w)
parmat[28, 1] <- -50

torch::with_no_grad({
  w$copy_(torch::torch_tensor(
    parmat,
    dtype  = w$dtype,
    device = w$device
  ))
})

sim2 <- m1$simulate(env_dt2, init_cohort = NULL, device = "cpu", patches = 100, patch_size = 0.06)
# sim2 <- m1$simulate(env = env_dt2, patch_size = 0.1, patches = 4, init_cohort = init_cohorts, device = "cpu")
sim_long_out <- copy(sim2$long$site)
sim_wide_out <- copy(sim2$wide$site)


p_dat <- sim_long_out[, .(value = mean(value)), by = .(year, species, variable)]
# p_dat <- sim_long_out
p_dat <- merge(p_dat, species_dt, by = "species", all.x = TRUE)
p_dat[variable %in% c("ba","trees"), value := value/0.06]


# sim_wide_out[year == 1000, .(ba = sum(ba)), .(species)][order(ba)]$species
sim_wide_out[, ba_rank := as.integer(factor(ba,ordered = T)), .(siteID,year)]
selected_species <- unique(sim_wide_out[ba > 2.5*0.06,]$species)
hist(sim_wide_out$ba)
ggplot(p_dat[species %in% selected_species],
       aes(x = year, y = value, color = species_name)) +
  geom_line() +
  facet_wrap(~variable, scales = "free")

ggplot(p_dat[species %in% selected_species & variable == "ba"],
       aes(x = year, y = value, color = species_name)) +
  geom_line() +
  facet_wrap(~variable, scales = "free_x")
ggplot(p_dat[variable == "ba"],
       aes(x = year, y = value)) +
  geom_line() +
  ggplot2::ylab("Basal Area (m2/ha)")+
  facet_wrap(~species_name, scales = "free_x")

ggplot(p_dat,
       aes(x = year, y = value, color = as.character(species_name))) +
  geom_line() +
  facet_wrap(~variable, scales = "free")

# ## 12) Save (optional) ----
# fwrite(site_dt, "vignettes/FIAsubset_sites.csv")
# fwrite(tree_dt, "vignettes/FIAsubset_trees.csv")
# fwrite(env_dt,  "vignettes/FIAsubset_env.csv")
# fwrite(obs_dt,  "vignettes/FIAsubset_obs.csv")
# if (!is.null(init_dt)) fwrite(init_dt, "vignettes/FIAsubset_init.csv")

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

par_mort <- data.table(m1$parameters_r$par_mortality_unconstrained)
par_mort$species_name <- species_dt$species_name[1:15]
colnames(par_mort) <- c("mortLight", "mortSize", "mortGrowth","species")

par_growth <- data.table(m1$parameters_r$par_growth_unconstrained)
par_growth$species_name <- species_dt$species_name[1:15]
colnames(par_growth) <- c("growthLight", "growthSize","species")

par_reg <- data.table(m1$parameters_r$par_regeneration_unconstrained)
par_reg$species_name <- species_dt$species_name[1:15]
colnames(par_reg) <- c("regLight","species")

par_mort <- melt(data.table(par_mort), measure.vars = colnames(par_mort)[-length(colnames(par_mort))])
par_growth <- melt(data.table(par_growth), measure.vars = colnames(par_growth)[-length(colnames(par_growth))])
par_reg <- melt(data.table(par_reg), measure.vars = colnames(par_reg)[-length(colnames(par_reg))])

process_par_dt <- rbindlist(list(
  data.table(par_mort, process = "mortality"),
  data.table(par_growth, process = "growth"),
  data.table(par_reg, process = "regeneration")
))

ggplot(process_par_dt, aes(x = value, y = species, fill = process))+
  geom_bar(stat = "identity", position = position_dodge2())+
  facet_wrap(~variable,nrow = 1)





