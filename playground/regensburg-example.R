library(FINN)
library(torch)
library(data.table)
library(ggplot2)

# the data read in here is created by the script
# "fit-fia-data/FIA_dataprep_vignette.R"
env_scales_dt <- fread("vignettes/fit-fia-data/env_scales_dt.csv")
tree_dt <- fread("vignettes/fit-fia-data/tree_dt.csv")
env_unscaled_dt <- fread("vignettes/fit-fia-data/env_unscaled_dt.csv")
env_dt <- fread("vignettes/fit-fia-data/env_dt.csv")
env_scales_dt <- fread("vignettes/fit-fia-data/env_scales_dt.csv")
obs_dt <- fread("vignettes/fit-fia-data/obs_dt.csv")
full_obs_dt <- fread("vignettes/fit-fia-data/full_obs_dt.csv")
full_tree_dt <- fread("vignettes/fit-fia-data/full_tree_dt.csv")
species_dt <- fread("vignettes/fit-fia-data/species_dt.csv")


init_trees <- full_tree_dt[year == 0,.(siteID, patchID, year, species, species_name, dbh, trees, living)]
init_trees[living != T | is.na(living), trees := 0,] # set number of trees of all non living trees to 0
summary(init_trees)
init_cohorts <- FINN::makeInitCohorts(init_trees, Nspecies = max(obs_dt$species))

m1 <- finn(
  N_species            = uniqueN(obs_dt$species),
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0, FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createProcess(~., FINN::growth,   optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~., FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~., FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)

m1 |> fit(
  env        = env_dt,
  data       = obs_dt,
  init_cohort = init_cohorts,
  device     = "cpu",
  epochs     = 20,
  batchsize  = 500L,
  patch_size = 0.06,
  patches    = 4, weights = c(0.1, 10, 1.0, 10.0, 1, 1),
  lr         = 0.01#, loss = c("mse","mse","mse","mse","mse","mse","mse")
)

predictions =
  m1 |> simulateForest(env = env_dt)

plot(m1, pars = "env", env_names = env_scales_dt$variable, species_names = species_dt$species_name)
plot(m1, pars = "process", species_names = species_dt$species_name)



