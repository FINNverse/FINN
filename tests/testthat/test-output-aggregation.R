library(data.table)
library(ggplot2)
library(FINN)

Ntimesteps = 500  # number of timesteps
Nsites = 1 # number of sites
patch_size = 0.1
# one species
Nsp = 10 # number of species

FINN.seed(1234)
torch::torch_manual_seed(1234)
# we draw the same shade parameters for each process for simplicity
# shade parameters correspond to the fraction of light a species needs to succesfully grow, regenerate, or survive.
shadeSP = runif(Nsp, 0.01, 0.99)

# regeneration parameters
parReg = shadeSP # regeneration is only dependent on shade and environment
parRegEnv = list(matrix(c(
  runif(Nsp, 1, 4), # intercept regulating the overall effect size
  runif(Nsp, -2, 2) # the second parameter modulates the effect of the environmental variable
  ),Nsp, 2))

# growth parameters
parGrowth = matrix(c(
  shadeSP, # see above
  runif(Nsp, 0, 0.1) # the second growth parameter modulates the size dependent growth
  ),Nsp, 2)

parGrowthEnv = list(matrix(c(
  runif(Nsp, -2, 1), # intercept regulating the overall effect size
  runif(Nsp, -1, 2) # the second parameter modulates the effect of the environmental variable
  ),Nsp, 2))

# mortality parameters
parMort = matrix(c(
  shadeSP, # see above
  runif(Nsp, 1, 4) # the second growth parameter modulates the size dependent mortality
  ),Nsp, 2)
parMortEnv = list(matrix(c(
  runif(Nsp, 0, 0.5), # intercept regulating the overall effect size
  runif(Nsp, 0, .5) # the second parameter modulates the effect of the environmental variable
  ), Nsp, 2))

# allometric parameters for the calculation of tree height from a trees diameter
parHeight = runif(Nsp, 0.3, 0.7)

# we first generate a data.table with all combinations of site and timestep.
env_dt <- data.table(
  expand.grid(
    list(
      siteID = 1:Nsites,
      year = 1:Ntimesteps
    )
  )
)

dist_dt <- env_dt

# for this very simple model we will have a constant environment for all sites and timesteps
env_dt$env1 = rep(0, Ntimesteps)

# we can also specify the intensity of disturbances for each timestep

# here we specify a disturbance frequency of 1%, which means that there is a 1% chance each year that a disturbance occurs
disturbance_frequency = 0.2

# the disturbance intensity at each timestep is the fraction of patches that is disturbed at that time step
# here we specify a disturbance frequency of 1%, which means that there is a 1% chance each year that a disturbance occurs
disturbance_intensity = runif(Ntimesteps*Nsites, 0, 0.05)

# this will result in 0 to 20 % of the patches being disturbed at each timestep with a change of 1% that a disturbance occurs at the timestep
dist_dt$intensity = dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, disturbance_frequency)*disturbance_intensity

predictions <- list()
predictions[["patches_1"]] =
  simulateForest(env = env_dt,
                 sp = Nsp,
                 # disturbance = dist_dt, this argument is optional and we will not use it for one patch
                 patches = 1,
                 patch_size = patch_size,
                 height = parHeight,
                 growthProcess = createProcess(~ 1 + env1, func = growth, initEnv = parGrowthEnv, initSpecies = parGrowth),
                 mortalityProcess = createProcess(~ 1 + env1, func = mortality, initEnv = parMortEnv, initSpecies = parMort),
                 regenerationProcess = createProcess(~ 1 + env1, func = regeneration, initEnv = parRegEnv, initSpecies = parReg),
                 debug = T
                 )

predictions[["patches_2"]] =
  simulateForest(env = env_dt,
                 sp = Nsp,
                 # disturbance = dist_dt, this argument is optional and we will not use it for one patch
                 patches = 2,
                 patch_size = patch_size,
                 height = parHeight,
                 growthProcess = createProcess(~ 1 + env1, func = growth, initEnv = parGrowthEnv, initSpecies = parGrowth),
                 mortalityProcess = createProcess(~ 1 + env1, func = mortality, initEnv = parMortEnv, initSpecies = parMort),
                 regenerationProcess = createProcess(~ 1 + env1, func = regeneration, initEnv = parRegEnv, initSpecies = parReg),
                 debug = T
                 )



predictions[["patches_100"]] =
  simulateForest(env = env_dt,
                 sp = Nsp,
                 patches = 100,
                 patch_size = patch_size,
                 disturbance = dist_dt,
                 height = parHeight,
                 growthProcess = createProcess(~ 1 + env1, func = growth, initEnv = parGrowthEnv, initSpecies = parGrowth),
                 mortalityProcess = createProcess(~ 1 + env1, func = mortality, initEnv = parMortEnv, initSpecies = parMort),
                 regenerationProcess = createProcess(~ 1 + env1, func = regeneration, initEnv = parRegEnv, initSpecies = parReg),
                 debug = T
                 )

for(i in c("patches_1", "patches_2", "patches_100")){
  # p_dat <- predictions[[i]]$long$site[variable == "ba", .(value = mean(value)), by = .(year, species, variable)]
  p_dat <- predictions[[i]]$long$site[, .(value = mean(value)), by = .(year, species, variable)]
  p_dat[, variable2 := factor(
    variable,
    levels = c("dbh", "ba", "trees", "AL", "growth", "mort", "reg", "r_mean_ha"),
    labels =  c("avg. DBH [cm]", "Basal Area [m²/ha]", "Trees [N/ha]",
                "Available Light [%]", "Growth [cm/yr]", "Mortality [%]",
                "Reg. Count [N/ha]", "Reg. Mean [N/ha]")
    ),]
  p_dat[variable %in% c("ba", "trees", "reg"), value := value/patch_size,]
  p <- ggplot(p_dat, aes(x = year, y = value, color = factor(species))) +
  geom_line() +
  theme_minimal() +
  labs(x = "Year", y = "Value") +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~variable2, scales = "free_y", ncol = 2, strip.position = "left") +  # Remove label_parsed
  theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90)
  ) +
  guides(color = guide_legend(title = "Species", override.aes = list(linewidth = 5), ncol = 2, title.position = "top")) +
  scale_color_discrete(name = "Species") +
  ggtitle(paste0("patches = ",gsub("patches_", "", i)))

print(p)
}


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## test output aggregation ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
out <- predictions[[2]]
Npatches = max(out$wide$patch$patchID)
cohort_dt <- out$wide$cohort
# patch_dt = cohort_dt[,.(
#   # dbh = sum(dbh*(trees-m*trees), na.rm = T)/sum(trees-m*trees, na.rm = T),
#   dbh = sum(dbh*trees, na.rm = T),
#   ba = sum(BA_stem(dbh) * trees, na.rm = T),
#   trees = sum(trees, na.rm = T),
#   mort = sum(m*trees, na.rm = T)/sum(trees, na.rm = T),
#   growth = sum(g*trees, na.rm = T)/sum(trees, na.rm = T)
# ), by = .(siteID, patchID, species, year)]
#
# site_dt = patch_dt[,.(
#   dbh = sum(dbh, na.rm = T)/sum(trees),
#   ba = sum(ba, na.rm = T)/Npatches,
#   trees = sum(trees, na.rm = T)/Npatches,
#   mort = sum(mort, na.rm = T)/Npatches,
#   growth = sum(growth, na.rm = T)/Npatches
# ), by = .(siteID, year, species)]

site_dt = cohort_dt[,.(
  dbh = sum(dbh*(trees), na.rm = F)/sum(trees, na.rm = F)/Npatches,
  ba = sum(BA_stem(dbh)*trees, na.rm = T)/Npatches,
  trees = sum(trees, na.rm = T)/Npatches,
  # ba = sum(mort_count, na.rm = T),
  #mort = sum(m*trees, na.rm = T)/sum(trees, na.rm = T)/Npatches,
  mort = sum(m*trees, na.rm = T)/sum(trees, na.rm = T)/Npatches,

  # mort_count = sum(mort_count, na.rm = T),
  growth = sum(g*trees, na.rm = T)/sum(trees, na.rm = T)/Npatches,
  Ncohorts = sum(trees > 0, na.rm = T)/Npatches
), by = .(siteID, species, year)]


comp_dt <- merge(site_dt, patch_dt, by = c("siteID", "species", "year"), suffixes = c(".site", ".patch"))

ggplot(comp_dt, aes(x = dbh.site, y = dbh.patch)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "site dbh mean", y = "patch dbh mean") +
  ggtitle("dbh mean site vs patch")
ggplot(comp_dt, aes(x = ba.site, y = ba.patch)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "site ba", y = "patch ba") +
  ggtitle("ba site vs patch")
ggplot(comp_dt, aes(x = trees.site, y = trees.patch)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "site trees", y = "patch trees") +
  ggtitle("trees site vs patch")
ggplot(comp_dt, aes(x = mort.site, y = mort.patch)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "site m", y = "patch m") +
  ggtitle("m site vs patch")
ggplot(comp_dt, aes(x = growth.site, y = growth.patch)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "site g", y = "patch g") +
  ggtitle("g site vs patch")


finn_site_dt <- out$wide$site

comp_out_dt <- merge(
  site_dt, finn_site_dt,
  by = c("siteID", "species", "year"),
  suffixes = c(".site", ".finn")
  )

ggplot(comp_out_dt, aes(x = dbh.site, y = dbh.finn, color = year)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "site dbh mean", y = "finn dbh mean") +
  ggtitle("dbh mean site vs finn")+
  facet_wrap(~species)
ggplot(comp_out_dt, aes(x = ba.site, y = ba.finn, color = year)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "site ba", y = "finn ba") +
  ggtitle("ba site vs finn")+
  facet_wrap(~species)
ggplot(comp_out_dt, aes(x = trees.site, y = trees.finn, color = year)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "site trees", y = "finn trees") +
  ggtitle("trees site vs finn")+
  facet_wrap(~species)


comp_out_dt[,":="(
  diff_dbh = dbh.site - dbh.finn,
  diff_ba = ba.site - ba.finn,
  diff_trees = trees.site - trees.finn,
  diff_mort = mort.site - mort.finn,
  diff_growth = growth.site - growth.finn
  )]

fm = lm(diff_trees ~ year+mort.site+growth.site, data = comp_out_dt)
summary(fm)

