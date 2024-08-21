library(FINN)
library(data.table)
library(ggplot2)

Nsp = 3
Npatches = 10
Nsites = 100
Ntimesteps = 200L

site_dt <- data.table(
  expand.grid(
    list(
      siteID = 1:Nsites,
      year = 1:Ntimesteps
    )
  )
)

env_dt <- site_dt
env_dt$env1 = rep(seq(-2,2,length.out = Nsites), Ntimesteps)

parGrowth = cbind(0.2, rep(1.0, 3))
shade_func = function(light_steepness = 10,  light = 0.5, inter = 0.0) {
((1 / (1 + exp(-light_steepness * (light - parGrowth[,1][species]) + inter)) - 1 / (1 + exp(light_steepness * parGrowth[,1][species] + inter))) /
           (1 / (1 + exp(-light_steepness * (1 - parGrowth[,1][species]) + inter)) - 1 / (1 + exp(light_steepness * parGrowth[,1][species] + inter))))
}

inp = seq(0.0, 1.0, length.out = 100)
plot(inp, sapply(inp, function(l) shade_func(light = l, )))
points(inp, sapply(seq(0.0, 1.0, length.out = 100), function(l) shade_func(light = l, light_steepness = 10, inter = -4)))
points(inp, sapply(seq(0.0, 1.0, length.out = 100), function(l) shade_func(light = l, light_steepness = 10, inter = + 4)))



growthCustom = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F){
  intercept = torch_tensor(c(0, -4.0, +4.0), device=self$device)

  shade = ((1 / (1 + torch::torch_exp(-light_steepness * (light - parGrowth[,1][species]) + intercept[species])) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]  + intercept[species]))) /
             (1 / (1 + torch::torch_exp(-light_steepness * (1 - parGrowth[,1][species])  + intercept[species])) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]  + intercept[species]))))

  environment = torch::torch_exp(pred)
  growth = shade * environment * (torch::torch_exp(-(dbh / (parGrowth[,2][species] * 100))))
  return(growth)
}

predictions =
  simulateForest(env = env_dt,
                 sp = Nsp,
                 patches=10L,
                 growthProcess = createProcess(~1+env1, func = growthCustom, initEnv = list(matrix(c(1.0, 1.0,
                                                                                               1.0, -1.0,
                                                                                               1.0, 0), Nsp, 2, byrow = TRUE)),
                                               initSpecies = matrix(c(0.8, 1.0), Nsp, 2, byrow = TRUE)),
                 mortalityProcess = createProcess(~1+env1, func = mortality, initEnv = list(matrix(c(0.0, -2.0,
                                                                                                     0.0, 2.0,
                                                                                                     2.0, 0.0), Nsp, 2, byrow = TRUE)),
                                                  initSpecies = matrix(c(0.8, 3.9), Nsp, 2, byrow = TRUE)),
                 regenerationProcess = createProcess(~1+env1, func = regeneration, initEnv = list(matrix(c(1.99, 0), Nsp, 2, byrow = TRUE)) ),
                 device = "cpu")


pred = predictions$long$site
ggplot(pred[siteID == 70], aes(x = year, y = value, group = siteID, color = factor(species))) +
  geom_line() +
  labs(x = "Year",
       y = "value") +
  theme_minimal()+
  facet_wrap(~variable, scales = "free_y")




# Hybrid modeling
growthCustom = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F){
  input =
    torch::torch_cat(
      list(dbh$unsqueeze(4L)/200., # rescale dbh
           pred$unsqueeze(2L)$unsqueeze(2L)$`repeat`(c(1L, dbh$shape[2], dbh$shape[3], 1L)),
           light$unsqueeze(4L),
           torch::nnf_one_hot(species, num_classes = 3L)),
      dim = 4L)
  growth = self$nnGrowthEnv(input)[,,,1]$exp()
  return(growth)
}

growthProcess = createProcess(~1+env1,
                              func = growthCustom,
                              hidden = c(50L, 50L),
                              initSpecies = matrix(c(0.8, 1.0), Nsp, 2, byrow = TRUE),
                              optimizeSpecies = FALSE,
                              optimizeEnv = TRUE,
                              inputNN = 7,
                              outputNN = 1)



data <- predictions$wide$site
predictions$model$check()
fit = finn(data = data,
           env = env_dt,
           patches = Npatches,
           height = as.numeric(predictions$model$parHeightT),
           growthProcess = growthProcess,
           mortalityProcess = createProcess(~1+env1,
                                            func = mortality,
                                            #initEnv = list(matrix(c(0.0, 0, 0.0, 0.0, 2.0, 0.0), sp, 2, byrow = TRUE)),
                                            hidden = c(50L, 50L),
                                            initSpecies = matrix(c(0.8, 3.9), Nsp, 2, byrow = TRUE),
                                            optimizeSpecies = FALSE, optimizeEnv = TRUE),
           regenerationProcess = createProcess(~1+env1, func = regeneration,
                                               initEnv = list(matrix(c(0.0, 0), Nsp, 2, byrow = TRUE)),
                                               initSpecies = (as.numeric(predictions$model$parRegT)),
                                               optimizeSpecies = FALSE, optimizeEnv = FALSE),
           weights = c(1/20, 1/2, 3,0.5,3,2),
           device = "cpu",
           optimizeHeight= FALSE,
           lr = 0.03,
           epochs = 5L,
           batchsize = 100L,
           speciesPars_ranges = predictions$model$speciesPars_ranges
)

continue_fit(fit, batchsize = 100L, lr = 0.01, epochs = 50L, weights = c(1/20, 1/2, 3,0.5,3,2, 1))



par(mfrow = c(1,1))

input = cbind(200/200 ,1, unique(env_dt$env1), 1.0, 0, 1, 0)
res =
  sapply(seq(1, 300, by = 50), function(i) {
    input = cbind(i/200 ,1, unique(env_dt$env1), 0.0, 1, 0, 0)
    ((fit$model$nnGrowthEnv( torch::torch_tensor( input , device = "cuda:0"))$relu()+0.001) |> as.matrix())[,1]

  })
matplot(x =unique(env_dt$env1) ,y = res, type = "l", col = viridis::viridis(6), lty = 1, lwd = 2.0, ylab = "growth", xlab = "Env")
legend("bottomright", bty = "n", col = viridis::viridis(6), lty = 1,legend = seq(1, 300, by = 50))
for(i in 1:3) {
  input = cbind(200/200 ,1, unique(env_dt$env1), 0.0, 0, 1, 0)
  plot(unique(env_dt$env1), (fit$model$nnGrowthEnv( torch::torch_tensor( input , device = "cuda:0"))$relu()+0.001 |> as.matrix())[,1])
  #points(unique(env_dt$env1), (predictions$model$nnGrowthEnv( torch::torch_tensor( input )) |> as.matrix())[,i] |> exp(), col = "red")

  plot(unique(env_dt$env1), (fit$model$nnMortEnv( torch::torch_tensor(cbind(1, unique(env_dt$env1)) , device = "cuda:0")) |> as.matrix())[,i] |> plogis(), ylim = c(0,1))
  points(unique(env_dt$env1), (predictions$model$nnMortEnv( torch::torch_tensor(cbind(1, unique(env_dt$env1)) )) |> as.matrix())[,i] |> plogis(), col = "red")
}


