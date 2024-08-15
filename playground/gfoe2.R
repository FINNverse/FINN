library(FINN)
library(data.table)

# one species, 100 patches, 50 sites, 1 env

Nsp = 1
Npatches = 100
Nsites = 200
Ntimesteps = 200L

site_dt <- data.table(
  expand.grid(
    list(
      siteID = 1:Nsites,
      year = 1:Ntimesteps
    )
  )
)

dist_dt <- site_dt
dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, 0.1)*rbeta(Ntimesteps*Nsites, 2, 2)

env_dt <- site_dt
env_dt$env1 = rep(seq(-2,2,length.out = Nsites), Ntimesteps)
sp = 1

system.time({
  predictions =
    simulateForest(env = env_dt,
                   #disturbance = dist_dt,
                   sp = 1L,
                   patches=10L,
                   growthProcess = createProcess(~1+env1, func = growth, initEnv = list(matrix(c(4, 0), sp, 2, byrow = TRUE)), initSpecies = matrix(c(0.8, 3.9), sp, 2, byrow = TRUE)),
                   mortalityProcess = createProcess(~1+env1, func = mortality, initEnv = list(matrix(c(1, -3), sp, 2, byrow = TRUE)), initSpecies = matrix(c(0.2, 3.9), sp, 2, byrow = TRUE)),
                   regenerationProcess = createProcess(~1+env1, func = regeneration, initEnv = list(matrix(c(3, 0), sp, 2, byrow = TRUE)) ),
                   device = "cpu")
})

predictions$Predictions
init = CohortMat$new(sp = 1L, dims = c(Nsites, Npatches, 5L))
pred = (predictions$model$nnGrowthEnv( torch::torch_tensor(cbind(1, unique(env_dt$env1)))))

light = predictions$model$competitionFunction(init$dbh, init$species, init$trees, predictions$model$get_parHeight(),
                                              h = NULL,
                                              minLight = predictions$model$minLight,
                                              patch_size_ha = predictions$model$patch_size_ha,
                                              ba = NULL,
                                              cohortHeights = NULL)

predictions$model$growthFunction(init$dbh+5.0, init$species,
                                 parGrowth = predictions$model$get_parGrowth(),
                                 pred = FINN::index_species(pred, init$species), light = light)

predictions$model$mortalityFunction(init$dbh, init$species,init$trees,
                                 parMort = predictions$model$get_parMort(),
                                 pred = FINN::index_species(pred, init$species), light = light)





predictions$model$nnGrowthEnv$parameters
predictions$model$nnMortEnv$parameters
predictions$model$get_parGrowth()
predictions$model$get_parMort()
predictions$model$get_parReg()

plot(unique(env_dt$env1), (predictions$model$nnGrowthEnv( torch::torch_tensor(cbind(1, unique(env_dt$env1)) )) |> as.matrix() |> exp())[,1])
plot(unique(env_dt$env1), predictions$model$nnMortEnv( torch::torch_tensor(cbind(1, unique(env_dt$env1)) )) |> as.matrix(), ylim = c(0,1))
plot(unique(env_dt$env1), (predictions$model$nnRegEnv( torch::torch_tensor(cbind(1, unique(env_dt$env1)) )) |> as.matrix() |> exp())[,1])

data <- pred2DF(predictions, format = "wide")$site

fit = finn(data = data,
           env = env_dt,
           #disturbance = dist_dt,
           growthProcess = createProcess(~1+env1,
                                         func = growth,
                                         initEnv = list(matrix(c(0, 0), sp, 2, byrow = TRUE)),
                                         initSpecies = matrix(c(0.2, 3.9), sp, 2, byrow = TRUE),
                                         optimizeSpecies = FALSE, optimizeEnv = TRUE),
           mortalityProcess = createProcess(~1+env1,
                                            func = mortality,
                                            initEnv = list(matrix(c(0, 0), sp, 2, byrow = TRUE)),
                                            initSpecies = matrix(c(0.2, 3.9), sp, 2, byrow = TRUE),
                                            optimizeSpecies = FALSE, optimizeEnv = TRUE),
           regenerationProcess = createProcess(~1+env1, func = regeneration,
                                               initEnv = list(matrix(c(4, 0), sp, 2, byrow = TRUE)),
                                               optimizeSpecies = FALSE, optimizeEnv = FALSE),
           device = "gpu", optimizeHeight= FALSE, lr = 0.3, epochs = 10L, batchsize = 100L,
           )
fit$model$parameters
continue_fit(fit, epochs = 20L, batchsize = 100L)
fit$model$parameters

fit$models_list[[3]]$model$check()
fit$models_list[[3]]$model$parMortEnv

rr=
  sapply(fit$models_list, function(m) {
    m = m$model
    m$check()
    return(as.matrix(m$parGrowthEnv[[1]][1,2]))
  })

mean(rr)
fit$model$parameters

predictions = predict(fit)

predictions = melt(predictions, id.vars = c("siteID", "species", "year"))

pred = pred2DF(predictions, format = "long")$site

par(mfrow = c(1,2))
plot(unique(env_dt$env1), (fit$model$nnGrowthEnv( torch::torch_tensor(cbind(1, unique(env_dt$env1)) , device = "cuda:0")) |> as.matrix())[,1], ylim = c(0,1))
points(unique(env_dt$env1), (predictions$model$nnGrowthEnv( torch::torch_tensor(cbind(1, unique(env_dt$env1)) )) |> as.matrix())[,1], col = "red")

plot(unique(env_dt$env1), (fit$model$nnMortEnv( torch::torch_tensor(cbind(1, unique(env_dt$env1)) , device = "cuda:0")) |> as.matrix())[,1], ylim = c(0,1))
points(unique(env_dt$env1), predictions$model$nnMortEnv( torch::torch_tensor(cbind(1, unique(env_dt$env1)) )) |> as.matrix(), col = "red")
as.matrix(fit$model$parameters[[1]])

pred = pred2DF(predictions, format = "long")$site

library(ggplot2)
# ggplot(pred[, .(value = mean(value) ), by = .(year, species, variable)], aes(x = year, y = value, color = factor(species))) +
ggplot(pred[siteID == 1], aes(x = year, y = value, group = siteID, color = siteID)) +
  geom_line() +
  labs(x = "Year",
       y = "value") +
  theme_minimal()+
  facet_wrap(~variable, scales = "free_y")

ggplot(pred[patch == 1, .(value = value ), by = .(year, species, variable)], aes(x = year, y = value, color = factor(species))) +
  geom_line() +
  labs(x = "Year",
       y = "value") +
  theme_minimal()+
  facet_wrap(~variable, scales = "free_y")


pred = pred2DF(predictions, format = "wide")$site
comp_dt <- merge(pred, env_dt, by = c("siteID", "year"))
plot(growth~env1, data = comp_dt)

summary(lm(growth~env1, data = comp_dt))




