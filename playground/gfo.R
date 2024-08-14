library(FINN)
library(data.table)

# one species, 100 patches, 50 sites, 1 env

Nsp = 1
Npatches = 100
Nsites = 50
Ntimesteps = 100

site_dt <-
  data.table(
    expand.grid(
      list(
        species = 1:Nsp,
        siteID = 1:Nsites,
        year = 1:Ntimesteps,
        intensity = rbinom(Ntimesteps, 1, 0.1)*rbeta(Ntimesteps, 2, 2)
        )
      )
    )

env_dt <- data.table(
  expand.grid(
    list(
      env1 = seq(-2,2,length.out = Nsites)
    )
  )
)
env_dt$siteID <- 1:nrow(env_dt)

sim_dt <- merge(site_dt, env_dt, by = "siteID")
dist_dt <- site_dt

disturbance = rbinom(Ntimesteps, 1, 0.1)*rbeta(Ntimesteps, 2, 2)

system.time({
  predictions =
    simulateForest(env = sim_dt,
                   disturbance = dist_dt,
                   sp = 1L,
                   patches=100L
                   ,device = "cpu")
})


FINN::simulateForest()



finn = f(sp = sp, env = 2L, device = "cuda:0", which = "all" ,
                # parGrowth = matrix(c(0.8, 15), sp, 2, byrow = TRUE),
                # parMort = matrix(c(runif(sp), runif(sp, 1, 4)), sp, 2, byrow = FALSE),
                # parReg = runif(sp, 0.8, 0.9), # any value between 0 and 1. 0 = species needs no light for regeneration, 1 = species needs full light for regeneration
                # parHeight = runif(sp, 0.3, 0.7), # plausible range is 0 to 1. default should be between 0.3 and 0.7
                # parGrowthEnv = list(matrix(c(10, 10, -5), sp, 2)),
                # parMortEnv = list(matrix(c(-1, -1, 1)*0, sp, 2)),
                # parRegEnv = list(matrix(c(1, 2, 3), sp, 2)),
                hidden_growth = c(5L, 5L),
                patch_size_ha = 0.1)
env = torch::torch_randn(size = c(sites, 60L, 2))
# env = torch::torch_zeros(size = c(sites, 100, 2))

system.time({

  pred = finn$predict(dbh = initCohort$dbh[,,],
                      trees = initCohort$trees[,,],
                      species = initCohort$species[,,], env = env[,,], patches = patches, debug = FALSE, verbose = TRUE)
})







