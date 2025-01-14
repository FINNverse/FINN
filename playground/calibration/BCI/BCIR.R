library(data.table)
library(FINN)
library(torch)

stand_dt = fread("processed_data/calibration-data/BCI-stand-data/stand_dt.csv")
meteo_dt = fread("processed_data/calibration-data/BCI-stand-data/meteo_dt.csv")
swp_dt = fread("processed_data/calibration-data/BCI-stand-data/swp_dt.csv")
cohorts_dt = fread("processed_data/calibration-data/BCI-stand-data/initial_cohorts.csv")

obs_dt = stand_dt[,.(
  # year = census - 1984,
  year = census,
  siteID, species,
  ba = ba,
  dbh = dbh_cm,
  trees = trees,
  growth = g,
  mort = m,
  reg = r,
  AL = NA_real_)]

# fill sites
# Jedes jahr alle sites!
empty_res =
  lapply(unique(obs_dt$species), function(sp) {
    empty_sites =
      lapply(unique(obs_dt$year), function(yy) {
        expand.grid(siteID = setdiff(obs_dt$siteID |> unique(), obs_dt[species==sp & year == yy]$siteID),
                    year = yy,
                    species = sp
        )
      })
    rbindlist(empty_sites, fill = TRUE)
  })
obs_dt2 = rbindlist(list(obs_dt, rbindlist(empty_res, fill = TRUE)), fill = TRUE)
obs_dt = obs_dt2

# nrow(obs_dt[species == 1])
# nrow(obs_dt[species == 2])

all(table(obs_dt$species) == unique(table(obs_dt$species)))
all(table(obs_dt$siteID) == unique(table(obs_dt$siteID)))
all(table(obs_dt$year) == unique(table(obs_dt$year)))


env_dt =
  lapply(unique(stand_dt$siteID), function(site) {
    tmp = meteo_dt
    tmp$siteID = site
    return(tmp)
  }) |> rbindlist()

env_dt <- merge(env_dt, swp_dt[,.(siteID, swp)], by = "siteID")

env_dt$Prec = scale(env_dt$Prec)
env_dt$T_mean = scale(env_dt$T_mean)
env_dt$T_max = scale(env_dt$T_max)
env_dt$T_min = scale(env_dt$T_min)
env_dt$SR_kW_m2 = scale(env_dt$SR_kW_m2)
env_dt$RH_prc = scale(env_dt$RH_prc)
env_dt$swp = scale(env_dt$swp)

str(env_dt)
str(obs_dt)

# Preparation of cohorts! BUG!!!!
cohorts_dt <- cohorts_dt[,.(cohortID, species, dbh = round(dbh_cm,4), trees)]
cohort1 <- FINN::CohortMat$new(dbh = cohorts_dt$dbh, trees = cohorts_dt$trees, species = cohorts_dt$species, sp = uniqueN(obs_dt$species))
# cohort2 = FINN::CohortMat$new(dbh = cohort1$dbh_r, trees = cohort1$trees_r, species = sp, sp = uniqueN(stand_dt$species))
# sp zum indexen, aber wenn 0 exisitiert -> seg fault / index fault

# Option A) Train model
m = finn(env = env_dt,
         data = obs_dt,
         patches = 1L,
         mortalityProcess = createProcess(~., func = mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
         growthProcess = createProcess(~., func = growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
         regenerationProcess = createProcess(~., func = regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
         device = "gpu",
         optimizeHeight = TRUE,
         init = cohort1,
         lr=0.01,
         epochs = 1,
         patch_size = 0.1,
         batchsize = 350L,sp = uniqueN(obs_dt$species),
         weights = c(1, 0.1, 3, 1.5, 3, 1.0))

# m2 = finn(env = env_dt, data = obs_dt,
#           mortalityProcess = createProcess(~., func = mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE, initSpecies = m$model$parMortTR, initEnv = m$model$parMortEnv_r),
#           growthProcess = createProcess(~., func = growth, optimizeSpecies = TRUE, optimizeEnv = TRUE, initSpecies = m$model$parGrowthTR, initEnv = m$model$parGrowthEnv_r),
#           regenerationProcess = createProcess(~., func = regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE, initSpecies = m$model$parRegTR, initEnv = m$model$parRegEnv_r),
#           device = "gpu",
#           height = as.vector(m$model$parHeightTR),
#           optimizeHeight = TRUE,
#           init = cohort2,
#           lr=0.01,
#           epochs = 1000L,
#           patch_size = 1,
#           batchsize = 350L,
#           weights = c(1, 0.1, 3, 1.5, 3, 1.0),
#           file = "BCI_parameters_1ha.RDS")

# m = readRDS(file = "gfoe2024/bci_model.RDS")
# m$model$check(device = "cpu")
# m$model$device = "cpu"
# m$model$device_r = "cpu"

# pred = predict(m)

library(ggplot2)
# pred$species = as.factor(pred$species)
# ggplot(pred[siteID == 25 & (species %in% 120:195)], aes(x = year, y = value,  color = (species))) +
#   geom_line() +
#   labs(x = "Year",
#        y = "value") +
#   theme_minimal()+
#   facet_wrap(~variable, scales = "free_y")


pars = readRDS("data/calibration-output/BCI_parameters_1ha_1_9_thinned.RDS")

pars_sim = pars[length(pars)]

par_backtr = list()
# i=names(pars_sim[[1]])[5]
for(i in names(pars_sim[[1]])) {
  # check if i beginns with "nn" as characetrs
  if(substr(i, 1, 2) == "nn") {
    par_backtr[[i]] = pars_sim[[1]][[i]]
  }else{
    par_backtr[[i]] = torch::as_array(m$model$getPars(internalPar = pars_sim[[1]][[i]], parRange = m$model$speciesPars_ranges[[i]]))
  }
}

uniqueN(cohort2$species_r)

# obs_df2 = rbindlist(lapply(1:10, function(i) {
#   obs_df$patchID = i
#   return(obs_df)
# }))

cohort1 <- FINN::CohortMat$new(obs_df, sp = 195)
# sp= cohort1$species_r
# sp[sp==0] = 1L
# sp[is.na(sp)] = 1L
# cohort2 = FINN::CohortMat$new(dbh = cohort1$dbh_r, trees = cohort1$trees_r, species = sp, sp = 195)
pdat = data.table()
for(i in 1:10) {
predictions =
  simulateForest(env = env_dt,
                 height = as.vector(par_backtr$parHeight),
                 sp = nrow(par_backtr$parReg),
                 init = cohort1,
                 patch_size = 1,
                 patches = 1,
                 growthProcess = createProcess(~., func = growth, initEnv = par_backtr$nnGrowth.0.weight,initSpecies = par_backtr$parGrowth),
                 mortalityProcess = createProcess(~., func = mortality, initEnv = par_backtr$nnMort.0.weight,initSpecies = par_backtr$parMort),
                 regenerationProcess = createProcess(~., func = regeneration, initEnv = par_backtr$nnReg.0.weight,initSpecies = as.vector(par_backtr$parReg)),
                 device = "cpu")

tmp_dt <- predictions$wide$site
tmp_dt$rep = i
pdat = rbind(pdat, tmp_dt)
}
pdat[,year := year+1984,]
pdat <- merge(stand_dt, pdat, by = c("siteID", "year", "species"))
library(ggplot2)
ggplot() +
  geom_line(aes(x = year, y = ba, group = interaction(rep,species), color = factor(species)), data = pdat[, .(ba = sum(ba.y)), by = .(year, rep, species, siteID)], alpha = 0.3) +
  geom_point(aes(x = year, y = ba, color = factor(species)), data = pdat[, .(ba = sum(ba.x)/100), by = .(year, siteID, species)])+
  labs(x = "Year",
       y = "value") +
  theme_minimal()+
  theme(legend.position = "none")+
  facet_wrap(~siteID)






