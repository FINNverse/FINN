library(FINN)
library(data.table)

# one species, 100 patches, 50 sites, 1 env

Nsp = 1
shadeSP = runif(Nsp, 0.01, 0.99)
Npatches = 50
Nsites = 1
Ntimesteps = 1000

speciesPars = list(
  speciesID = 1:Nsp,
  parGrowth = matrix(c(
    shadeSP,
    runif(Nsp, 1, 2)
    ),Nsp, 2),
  parMort = matrix(c(
    shadeSP,
    runif(Nsp, 1, 4)
    ),Nsp, 2),
  parReg = shadeSP,
  parHeight = runif(Nsp, 0.3, 0.7),
  parGrowthEnv = list(matrix(c(
    runif(Nsp, -1, 1),
    runif(Nsp, -1, 1)
    ),Nsp, 2)),
  parMortEnv = list(matrix(c(
    runif(Nsp, 0, 0.5),
    runif(Nsp, 0, .5)
    ), Nsp, 2)),
  parRegEnv = list(matrix(c(
    runif(Nsp, 0, 2),
    runif(Nsp, -2, 2)
    ),Nsp, 2))
  )

site_dt <- data.table(
  expand.grid(
    list(
      siteID = 1:Nsites,
      year = 1:Ntimesteps
    )
  )
)

# hist(rbeta(Ntimesteps*Nsites, 1, 10), xlim = c(0,1))
dist_dt <- site_dt
# dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, 0.001)*1
dist_dt$intensity = 0
dist_dt$intensity[250] = 1
dist_dt$intensity[500] = 0.5
dist_dt$intensity[750] = 0.25
# dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, 0.01)*rbeta(Ntimesteps*Nsites, 1, 5)

env_dt <- site_dt
env_dt$env1 = rep(0, Ntimesteps)
# env_dt$env1 = rep(seq(-2,2,length.out = Nsites), Ntimesteps)

system.time({
  predictions =
    simulateForest(env = env_dt,
                   disturbance = dist_dt,
                   sp = Nsp,
                   patches=Npatches,
                   height = speciesPars$parHeight,speciesPars_ranges = speciesPars_ranges,
                   # growthProcess = createProcess(~1+env1, func = growth, initEnv = list(matrix(c(runif(),3,4), Nsp, 2)), initSpecies = matrix(c(0.2, 3.9), Nsp, 2)),
                   # mortalityProcess = createProcess(~1+env1, func = mortality, initEnv = list(matrix(c(1, -3), Nsp, 2)), initSpecies = matrix(c(0.2, 3.9), Nsp, 2)),
                   # regenerationProcess = createProcess(~1+env1, func = regeneration, initEnv = list(matrix(c(4, 0), 1, 2)), initSpecies = matrix(c(0.2), Nsp)),
                   growthProcess = createProcess(~1+env1, func = growth, initEnv = speciesPars$parGrowthEnv, initSpecies = speciesPars$parGrowth),
                   mortalityProcess = createProcess(~1+env1, func = mortality, initEnv = speciesPars$parMortEnv, initSpecies = speciesPars$parMort),
                   regenerationProcess = createProcess(~1+env1, func = regeneration, initEnv = speciesPars$parRegEnv, initSpecies = speciesPars$parReg),
                   device = "cpu", debug = F)
})

pred = pred2DF(predictions, format = "long")

library(ggplot2)
ggplot(pred$site[, .(value = mean(value) ), by = .(year, species, variable)], aes(x = year, y = value, color = factor(species))) +
  geom_line() +
  labs(x = "Year",
       y = "value") +
  theme_minimal()+
  coord_cartesian(ylim = c(0,NA))+
  facet_wrap(~variable, scales = "free_y")





