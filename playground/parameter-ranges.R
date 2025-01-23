library(FINN)
library(data.table)
library(ggplot2)
library(torch)

## Create test cohorts ####
dbh_shape = 5
trees_vec = dbh = c(0, 1, 10, 100)
dbh_vec = c(1,seq(20,100,20), 200, 300)
Nspecies_vec = c(1, 3, 6)
patches_vec = c(1, 10, 100)

# test base configuration for Nspecies and Npatches
test_cases_dt <-
  data.table(
    expand.grid(
      Nspecies = Nspecies_vec,
      Npatches = patches_vec
    )
  )

# create cohorts that are each repeated within each test case
test_cohorts_dt <-
  data.table(
    expand.grid(
      trees = trees_vec,
      dbh = dbh_vec
    )
  )
test_cohorts_dt$siteID = 1:nrow(test_cohorts_dt)

# combine test cases and cohorts
cohorts_dt = data.table()
for(i in 1:nrow(test_cases_dt)){
  Nspecies = test_cases_dt[i,]$Nspecies
  Npatches = test_cases_dt[i,]$Npatches
  test_cohort_dt_i = data.table()
  for(j in 1:nrow(test_cohorts_dt)){
    i_trees = test_cohorts_dt[j,]$trees
    i_dbh = test_cohorts_dt[j,]$dbh
    if(i_trees == 0 | i_dbh == 0){
      test_cohort_dt_j = data.table(
        trees = 0,
        dbh = 0,
        patchID = 1,
        cohortID = 1
      )
    }else{
      test_cohort_dt_j = rweibull_cohorts(
        trees = i_trees, dbh_shape = dbh_shape, dbh_scale = i_dbh,
        dbh_class_range = 1
        )
    }
    test_cohort_dt_j$species = sample(1:Nspecies, nrow(test_cohort_dt_j), replace = T)
    test_cohort_dt_j$siteID = j
    test_cohort_dt_i <- rbind(test_cohort_dt_i, test_cohort_dt_j, fill = T)
  }
  for(p in 1:Npatches){
    test_cohort_dt_i$patchID = p
    cohorts_dt =
      rbind(
        cohorts_dt,
        data.frame(
          simID = i,
          Npatches = Npatches,
          Nspecies = Nspecies,
          Weib_trees = i_trees,
          dbh_shape = dbh_shape,
          dbh_scale = i_dbh,
          test_cohort_dt_i
        ))
    }
  }

# create list with cohorts and site data
test_set_list = list()
for(i in 1:nrow(test_cases_dt)){
  i_test_set = cohorts_dt[simID == i]
  test_set_list[[i]] = list(
    "base_config" =  test_cases_dt[i,],
    "site_dt" = unique(i_test_set[,c("siteID", "Nspecies", "Npatches", "Weib_trees", "dbh_shape", "dbh_scale")]),
    "cohorts_dt" = i_test_set[,c("siteID", "patchID", "cohortID", "species", "trees", "dbh")],
    "cohorts" = CohortMat$new(i_test_set, sp = Nspecies)
    )
}

## Testing individual functions ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
### competition  ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

out_dt = data.table()
for(i in 1:length(test_set_list)){
  i_test_set = test_set_list[[i]]
  Nsp = i_test_set$base_config$Nspecies

  parComp = matrix(c(
    runif(Nsp, 0.3, 0.7), # parHeight
    runif(Nsp, 0.2, 0.2) # Competition strength
  ),Nsp, 2)

  parComp = torch_tensor(
    cbind(
      parComp[,1],parComp[,2]
    ), requires_grad=TRUE, dtype=torch_float32()
  )
  # cat(i, "started competition function\n")
  comp_out = competition(
    dbh = i_test_set$cohorts$dbh,
    species = i_test_set$cohorts$species,
    trees= i_test_set$cohorts$trees,
    parComp = parComp,
    h = NULL,
    patch_size_ha = 0.1
  )
  cohortHeights = height(i_test_set$cohorts$dbh, parComp[,1][i_test_set$cohorts$species])$unsqueeze(4)
  ba = BA_stand(i_test_set$cohorts$dbh, i_test_set$cohorts$trees, patch_size_ha = 0.1)
  # cat(i, "finished competition function\n")

  # cat("looping through results\n")
  comp_dt_temp <- data.table()
  for(i_site in unique(i_test_set$cohorts_dt$siteID)){
    # for(i_patch in unique(i_test_set$cohorts_dt[siteID == i_site]$patchID)){
    for(i_cohort in unique(i_test_set$cohorts_dt[siteID == i_site]$cohortID)){
      # cat("looping over site", i_site, "patch", i_patch, "cohort", i_cohort, "\n")
      comp_dt_temp <- rbind(
        comp_dt_temp,
        data.table(
          simID = i,
          siteID = i_site,
          patchID = 1:dim(comp_out)[2],
          cohortID = i_cohort,
          height = as.vector(torch::as_array(cohortHeights)[i_site,,i_cohort,]),
          ba = as.vector(torch::as_array(ba)[i_site,,i_cohort]),
          comp = as.vector(torch::as_array(comp_out)[i_site,,i_cohort])
        )
      )
    }
    # }
  }
  out_dt = rbind(
    out_dt,
    comp_dt_temp
    )
  cat(i, "finished\n")
}
out_dt2 <- merge(test_cohorts_dt, out_dt, by = c("siteID"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## Growth ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
out_dt_growth = data.table()
i=1
for(i in 1:length(test_set_list)){
  i_test_set = test_set_list[[i]]
  Nsp = i_test_set$base_config$Nspecies

  parGrowth = matrix(c(
    runif(Nsp, 0.3, 0.7), # parHeight
    runif(Nsp, 0, 0.05) # Competition strength
  ),Nsp, 2)

  envpred = seq(-2, 2, 1)

  parGrowth_dt <- data.table(
    species = 1:Nsp,
    parGrowth1 = parGrowth[,1],
    parGrowth2 = parGrowth[,2]
  )

  parGrowth = torch_tensor(
    cbind(
      parGrowth[,1],parGrowth[,2]
    ), requires_grad=TRUE, dtype=torch_float32()
  )
  # cat(i, "started competition function\n")
  light = competition(
    dbh = i_test_set$cohorts$dbh,
    species = i_test_set$cohorts$species,
    trees= i_test_set$cohorts$trees,
    parComp = parComp,
    h = NULL,
    patch_size_ha = 0.1
  )


  growth_dt_temp <- data.table()
  for(i_env in 1:length(envpred)){
    pred = FINN::index_species(torch_ones(i_test_set$cohorts$species$shape[1],Nsp), i_test_set$cohorts$species)*envpred[i_env]
    cat("envpred = ",envpred[i_env], "\n")
    growth_out = growth(
      dbh = i_test_set$cohorts$dbh,
      species = i_test_set$cohorts$species,
      parGrowth = parGrowth,
      light = light,
      pred = pred
    )
    cohortHeights = height(i_test_set$cohorts$dbh, parComp[,1][i_test_set$cohorts$species])$unsqueeze(4)
    ba = BA_stand(i_test_set$cohorts$dbh, i_test_set$cohorts$trees, patch_size_ha = 0.1)
    # cat(i, "finished competition function\n")

    # cat("looping through results\n")
    for(i_site in unique(i_test_set$cohorts_dt$siteID)){
      # for(i_patch in unique(i_test_set$cohorts_dt[siteID == i_site]$patchID)){
      for(i_cohort in unique(i_test_set$cohorts_dt[siteID == i_site]$cohortID)){
        # cat("looping over site", i_site, "patch", i_patch, "cohort", i_cohort, "\n")
        growth_dt_temp <- rbind(
          growth_dt_temp,
          data.table(
            simID = i,
            siteID = i_site,
            patchID = 1:dim(growth_out)[2],
            cohortID = i_cohort,
            height = as.vector(torch::as_array(cohortHeights)[i_site,,i_cohort,]),
            growth = as.vector(torch::as_array(growth_out)[i_site,,i_cohort]),
            ba = as.vector(torch::as_array(ba)[i_site,,i_cohort]),
            light = as.vector(torch::as_array(light)[i_site,,i_cohort]),
            species = as.vector(torch::as_array(i_test_set$cohorts$species)[i_site,,i_cohort]),
            envpred = envpred[i_env]
          )
        )
      }
    }
    # }
  }

  growth_dt_temp <- merge(growth_dt_temp, parGrowth_dt, by = c("species"))

  out_dt_growth = rbind(
    out_dt_growth,
    growth_dt_temp
  )
  cat(i, "finished\n")
}
out_dt_growth2 <- merge(test_cohorts_dt, out_dt_growth, by = c("siteID"))

ggplot(out_dt_growth2[light == 1 & species == 6 & growth < 1], aes(x = as.factor(parGrowth2), y = as.factor(envpred)))+
  geom_tile(aes(fill = growth))

ggplot(out_dt_growth2[ba < 50], aes(x = cut(parGrowth2, breaks = 10), y = growth))+
  geom_boxplot()+
  coord_cartesian(ylim = c(0,0.5))

summary(out_dt_growth2[growth < 0.2 & light > parGrowth1])
nrow(out_dt_growth2[growth < 0.2 & light > parGrowth1 & parGrowth2 < 1 & envpred == 0])
summary(out_dt_growth2[growth > 1 & light > parGrowth1 & parGrowth2 < 0.5])
nrow(out_dt_growth2[growth > 1 & light > parGrowth1])
nrow(out_dt_growth2[growth > 1 & light > parGrowth1 & parGrowth2 < 0.5 & envpred == 0])

ggplot(out_dt_growth2, aes(x = parGrowth2, y = growth, color = species))+
  geom_point()


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
### Regeneration  ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=


out_dt = data.table()
for(i in 1:length(test_set_list)){
  i_test_set = test_set_list[[i]]
  Nsp = i_test_set$base_config$Nspecies

  parReg = runif(Nsp, 0.01, 0.99)

  parComp = torch_tensor(parReg, requires_grad=TRUE, dtype=torch_float32())

  # cat(i, "started competition function\n")
  comp_out = regeneration(species =
  )
  cohortHeights = height(i_test_set$cohorts$dbh, parComp[,1][i_test_set$cohorts$species])$unsqueeze(4)
  ba = BA_stand(i_test_set$cohorts$dbh, i_test_set$cohorts$trees, patch_size_ha = 0.1)
  # cat(i, "finished competition function\n")

  # cat("looping through results\n")
  comp_dt_temp <- data.table()
  for(i_site in unique(i_test_set$cohorts_dt$siteID)){
    # for(i_patch in unique(i_test_set$cohorts_dt[siteID == i_site]$patchID)){
    for(i_cohort in unique(i_test_set$cohorts_dt[siteID == i_site]$cohortID)){
      # cat("looping over site", i_site, "patch", i_patch, "cohort", i_cohort, "\n")
      comp_dt_temp <- rbind(
        comp_dt_temp,
        data.table(
          simID = i,
          siteID = i_site,
          patchID = 1:dim(comp_out)[2],
          cohortID = i_cohort,
          height = as.vector(torch::as_array(cohortHeights)[i_site,,i_cohort,]),
          ba = as.vector(torch::as_array(ba)[i_site,,i_cohort]),
          comp = as.vector(torch::as_array(comp_out)[i_site,,i_cohort])
        )
      )
    }
    # }
  }
  out_dt = rbind(
    out_dt,
    comp_dt_temp
  )
  cat(i, "finished\n")
}
out_dt2 <- merge(test_cohorts_dt, out_dt, by = c("siteID"))


hist(out_dt$comp)
summary(cut(out_dt$height, breaks = 10))
boxplot(out_dt$comp ~ cut(out_dt$height, breaks = 100))

hist(out_dt$comp)
# growth parameters
parComp = torch_tensor(cbind(parComp[,1],parComp[,2]), requires_grad=TRUE, dtype=torch_float32())



# check output of full simulations

Ntimesteps = 200  # number of timesteps
Nsites = 1 # number of sites
patch_size = 0.1
# one species
Nsp = 5 # number of species


FINN.seed(1235)
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
  runif(Nsp, 0, 0.01) # the second growth parameter modulates the size dependent growth
),Nsp, 2)

parGrowthEnv = list(matrix(c(
  runif(Nsp, -1, 0.5), # intercept regulating the overall effect size
  runif(Nsp, -1, 0.5) # the second parameter modulates the effect of the environmental variable
),Nsp, 2))

# mortality parameters
parMort = matrix(c(
  shadeSP, # see above
  runif(Nsp, 1, 4) # the second growth parameter modulates the size dependent mortality
),Nsp, 2)
parMortEnv = list(matrix(c(
  runif(Nsp, 0, 0.5), # intercept regulating the overall effect size
  runif(Nsp, 0, 0) # the second parameter modulates the effect of the environmental variable
), Nsp, 2))

# allometric parameters for the calculation of tree height from a trees diameter
# parHeight = runif(Nsp, 0.3, 0.7)
# growth parameters
parComp = matrix(c(
  runif(Nsp, 0.3, 0.7), # parHeight
  runif(Nsp, 0.2, 0.2) # Competition strength
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
disturbance_intensity = rbinom(Ntimesteps*Nsites,1,0.2)*runif(Ntimesteps*Nsites, 0, 0.05)

# this will result in 0 to 20 % of the patches being disturbed at each timestep with a change of 1% that a disturbance occurs at the timestep
dist_dt$intensity = dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, disturbance_frequency)*disturbance_intensity

sim1 =
simulateForest(env = env_dt,
               sp = Nsp,
               patches = 100,
               patch_size = patch_size,
               competitionProcess = createProcess(func = competition, initSpecies = parComp),
               growthProcess = createProcess(~ 1 + env1, func = growth, initEnv = parGrowthEnv, initSpecies = parGrowth),
               mortalityProcess = createProcess(~ 1 + env1, func = mortality, initEnv = parMortEnv, initSpecies = parMort),
               regenerationProcess = createProcess(~ 1 + env1, func = regeneration, initEnv = parRegEnv, initSpecies = parReg),
               debug = T
)
summary(sim1$wide$site)

hist(sim1$wide$site$dbh)

summary(sim1$wide$cohort$dbh)
sim_dt <- sim1$wide$cohort

table(sim_dt[dbh > 500]$species)

sim_dt[,.(
  .N
  ), by = .(species,cut(dbh, c(seq(0,1000,200), Inf)))][order(cut, species)]


# combine all parameters in a data.table
species_pars = data.table()
for (i in 1:Nsp) {
  species_pars = rbind(species_pars, data.table(
    species = i,
    shade = shadeSP[i],
    reg = parReg[i],
    regEnv1 = parRegEnv[[1]][i,1],
    regEnv2 = parRegEnv[[1]][i,2],
    parGrowth1 = parGrowth[i,1],
    parGrowth2 = parGrowth[i,2],
    parGrowthEnv1 = parGrowthEnv[[1]][i,1],
    parGrowthEnv2 = parGrowthEnv[[1]][i,2],
    parMort1 = parMort[i,1],
    parMort2 = parMort[i,2],
    parMortEnv1 = parMortEnv[[1]][i,1],
    parMortEnv2 = parMortEnv[[1]][i,2],
    parComp1 = parComp[i,1],
    parComp2 = parComp[i,2]
  ))
}

species_pars

sim1$wide$site[AL > 1]
sim1$wide$site[AL <= 1]

sim1$wide$cohort
