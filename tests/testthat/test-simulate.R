# This script tests if the simulation function works as expected
library(testthat)
library(FINN)
library(data.table)
library(ggplot2)

Ntimesteps = 500  # number of timesteps
Nsites = 1 # number of sites
patch_size = 0.1
# one species
Nsp = 10 # number of species

# set the seed for reproducibility
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
  runif(Nsp, 0, 0.1) # the second growth parameter modulates the size dependent growth
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
parAllo = matrix(c(
  runif(Nsp, 0.3, 0.7), # height parameters
  runif(Nsp, 0.9, 1.1) # competition strength parameters
),Nsp, 2)

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
                 parAllo = parAllo,
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


