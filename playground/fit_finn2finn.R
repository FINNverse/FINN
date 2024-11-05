library(FINN)
library(data.table)
library(ggplot2)

Ntimesteps = 300  # number of timesteps
Nsites = 1 # number of sites
patch_size = 0.1
# one species
Nsp = 10 # number of species

## Species parameters

FINN.seed(1234)
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
  runif(Nsp, 1, 2) # the second growth parameter modulates the size dependent growth
),Nsp, 2)

parGrowthEnv = list(matrix(c(
  runif(Nsp, -1, 1), # intercept regulating the overall effect size
  runif(Nsp, -1, 1) # the second parameter modulates the effect of the environmental variable
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

## Environment and disturbances {#sec-env-dist}

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

## Simulate

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
                 debug = F
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
                 debug = F
  )

for(i in c("patches_1", "patches_100")){
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


## Fit model

data = predictions[["patches_1"]]$wide$site

# we will fit a model with the following structure
fit <- list()
fit[["patches_100"]] =
  finn(env = env_dt,
       data = data,
       sp = Nsp,
       patches = 100,
       patch_size = patch_size,
       disturbance = dist_dt,
       growthProcess = createProcess(~ 1 + env1, func = growth, optimizeSpecies = T),
       mortalityProcess = createProcess(~ 1 + env1, func = mortality, optimizeSpecies = T),
       regenerationProcess = createProcess(~ 1 + env1, func = regeneration, optimizeSpecies = T),
       optimizeHeight = T, epochs = 100
  )

# compare parameters from fit and original parameters
true_pars = list(
  parHeight = parHeight,
  parGrowth = parGrowth,
  parGrowthEnv = parGrowthEnv[[1]],
  parMort = parMort,
  parMortEnv = parMortEnv[[1]],
  parReg = parReg,
  parRegEnv = parRegEnv[[1]]
)

# make data.table from true_pars
true_pars = lapply(true_pars, function(x) as.data.table(x))
true_pars_dt <- data.table()
i = names(true_pars)[1]
for(i in names(true_pars)) {
  tmp_dt <- data.table(true_pars[[i]])
  names(tmp_dt) <- paste0(i, "_", 1:ncol(tmp_dt))
  true_pars_dt = cbind(true_pars_dt,tmp_dt)
}

estimated_pars = list(
  parHeight = fit[["patches_100"]]$model$parHeightTR,
  parGrowth = fit[["patches_100"]]$model$parGrowthTR,
  parGrowthEnv = fit[["patches_100"]]$model$parGrowthEnv_r[[1]],
  parMort = fit[["patches_100"]]$model$parMortTR,
  parMortEnv = fit[["patches_100"]]$model$parMortEnv_r[[1]],
  parReg = fit[["patches_100"]]$model$parRegTR,
  parRegEnv = fit[["patches_100"]]$model$parRegEnv_r[[1]]
)

# make data.table from estimated_pars
estimated_pars = lapply(estimated_pars, function(x) as.data.table(x))
estimated_pars_dt <- data.table()
i = names(estimated_pars)[1]
for(i in names(estimated_pars)) {
  tmp_dt <- data.table(estimated_pars[[i]])
  names(tmp_dt) <- paste0(i, "_", 1:ncol(tmp_dt))
  estimated_pars_dt = cbind(estimated_pars_dt,tmp_dt)
}


par(mfrow = c(3,3))
for(i in names(true_pars_dt)){
  {
    value_range = range(c(true_pars_dt[[i]], estimated_pars_dt[[i]]))
    # add plot as perfect square scaled to the range of both variables
    plot(true_pars_dt[[i]], estimated_pars_dt[[i]], xlab = "True", ylab = "Estimated",
         main = paste0(i, "R² = ", round(cor(true_pars_dt[[i]], estimated_pars_dt[[i]])^2, 2)),
         xlim = value_range, ylim = value_range)
    # fit model
    fm = lm(estimated_pars_dt[[i]] ~ true_pars_dt[[i]])
    # add line
    abline(fm, col = "red")
  }

}
