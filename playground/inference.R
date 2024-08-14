library(FINN)
library(data.table)
library(ggplot2)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## env data input ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# Define the number of days in each month
days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Define the years and site IDs
years <- 1970:2019
site_ids <- 1:10

# Create an empty list to store data
climate_data_list <- list()

# Loop over each year, month, and site to generate data
for (site_id in site_ids) {
  for (year in years) {
    for (month in 1:12) {
      days <- 1:days_in_month[month]

      # Handle leap years
      if (month == 2 && (year %% 4 == 0 && (year %% 100 != 0 || year %% 400 == 0))) {
        days <- 1:29
      }

      # Create a data frame for the current site, year, and month
      climate_data <- data.frame(
        siteID = site_id,
        year = year,
        month = month,
        day = days,
        tmp = runif(length(days), -10, 30),  # Random temperature values
        pre = runif(length(days), 0, 200)    # Random precipitation values
      )

      # Add to the list
      climate_data_list <- append(climate_data_list, list(climate_data))
    }
  }
}

# Combine the list into a single data frame
climate_dt <- rbindlist(climate_data_list)

# Use the function to create the array
climate_dt_day <- climate_dt
resultDay <- climateDF2array(climate_dt_day, env_vars = c("tmp", "pre"))
str(resultDay)

climate_dt_month <- climate_dt[, .(tmp = mean(tmp), pre = sum(pre)), by = .(siteID,year,month)]
resultMonth <- climateDF2array(climate_dt = data.frame(climate_dt_month), env_vars = c("tmp", "pre"))
str(resultMonth)

climate_dt_year <- climate_dt[, .(tmp = mean(tmp), pre = sum(pre)), by = .(siteID,year)]
resultYear <- climateDF2array(climate_dt = climate_dt_year, env_vars = c("tmp", "pre"))
str(resultYear)


library(FINN)

env = data.table(siteID = rep(1:50, each = 20), year = rep(1:20, 50), tmp = runif(50*20), pre = runif(50*20))
system.time({

  predictions =
simulateForest(env,
               sp = 2L,
                patches=50L
               ,device = "gpu")
})


system.time({

  predictions =
    simulateForest(env,
                   disturbance = disturbance,
                   sp = 2L,
                   patches=50L
                   ,device = "gpu")
})


plot(predictions$Predictions$Site$dbh[1,,2])

data <- pred2DF(predictions, format = "wide")$site
head(data)
env
library(FINN)
result = finn(data = data, env = env, device = "gpu", batchsize = 50L, lr = 0.1, epochs = 5L, parallel=10L, NGPU = 4L, bootstrap = 10L)
result$model$get_parGrowth()
result$models_list[[1]]$model$check()
result$models_list[[1]]$model$parameters

result = finn(data = data,
              env = env,
              mortalityProcess = createProcess(~., optimizeAllometric = TRUE, optimizeEnv = FALSE, func = mortality),
              device = "gpu", batchsize = 50L, lr = 0.1, epochs = 5L)
result$model$parameters
preds = predict(result)


result$model$get_parHeight()
matplot(sapply(1:5, function(i) result$model$param_history[[i]]$parHeight) |> t(), type = "l")
preds = predict(result, env = env)
library(dplyr)
tmp = preds |> filter(species ==2) |> filter(siteID == 1)
plot(dbh~year, data = tmp)



rr=
sapply(result$models_list, function(m) {
  m = m$model
  m$check()
  return(as.matrix(m$parGrowthEnv[[1]][1,1]))
})
apply(rr, 1, mean)
apply(rr, 1, sd)
mean(rr)
sd(rr)
hist(rnorm(1000, mean(rr), sd(rr)))
result$models_list[[1]]$model

# with disturbance
disturbance = data.table(siteID = rep(1:50, each = 20), year = rep(1:20, 50), intensity = sample(c(0.0, 0.2), size = 20*50, replace = TRUE) )

result = finn(data = data,
              env = env,
              patches = 1L,
              disturbance = disturbance,
              mortalityProcess = createProcess(~., optimizeAllometric = TRUE, optimizeEnv = FALSE, func = mortality),
              device = "gpu", batchsize = 50L, lr = 0.1, epochs = 5L)
result$model$parameters
preds = predict(result)
preds

tmp = preds |> filter(species ==2) |> filter(siteID == 6)
plot(trees~year, data = tmp)







library(parallel)
library(FINN)
cl = makeCluster(3L)
sp = 3L
patches = 100L
sites = 500L
initCohort = CohortMat$new(dims = c(sites, patches, 10),
                           dbh = array(1, dim = c(sites, patches, 10)),
                           trees = array(1, dim = c(sites, patches, 10)),
                           sp = sp)

parallel::clusterEvalQ(cl, {library(torch);library(FINN)})
k = parallel::clusterExport(cl, "initCohort")
k = parallel::clusterEvalQ(cl, {initCohort$check()})

parallel::clusterEvalQ(cl, {rr = CohortMat$new(dims = c(10, 4, 10),
                                          dbh = array(1, dim = c(10, 4, 10)),
                                          trees = array(1, dim = c(10, 4, 10)),
                                          sp = 3)})

parallel::clusterEvalQ(cl, {initCohort$check()})
parLapply(cl, 1:3, fun = function(i) {
  p = initCohort$check()
  return(as_array(initCohort$dbh))
})


library(FINN)
#torch::torch_set_num_interop_threads(1L)
#torch::torch_set_num_threads(1L)
sp = 3L
patches = 10L
sites = 50L
initCohort = CohortMat$new(dims = c(sites, patches, 10),
                           dbh = array(1, dim = c(sites, patches, 10)),
                           trees = array(1, dim = c(sites, patches, 10)),
                           sp = sp)


finn = FINN$new(sp = sp, env = 2L, device = "cuda:0",
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
                    species = initCohort$species[,,], env = env[,,],
                    disturbance = torch::torch_rand(sites, 60L, 1L)$to(device = "cuda:0"),
                    patches = patches, debug = FALSE, verbose = TRUE)
})



finn$check()

library(FINN)
library(parallel)
cl = makeCluster(10L)

parallel::clusterEvalQ(cl, {library(torch);library(FINN)})
k = parallel::clusterExport(cl, c("initCohort", "finn"))

system.time({
parLapply(cl, 1:7, fun = function(i) {
  initCohort$check()
  finn$check()
  patches = 100L
  sites = 50L
  env = torch::torch_randn(size = c(sites, 60L, 2))
  pred = finn$predict(dbh = initCohort$dbh[1:50,,],
                      trees = initCohort$trees[1:50,,],
                      species = initCohort$species[1:50,,], env = env[1:50,,], patches = patches, debug = FALSE, verbose = TRUE)
  return(1)
})
})
parallel::stopCluster(cl)
parallel::clusterEvalQ(cl, {finn})

library(FINN)
finn$check()


# 1 -> dbh/ba, 2 -> counts, 3 -> AL, 4 -> growth rates, 5 -> mort rates, 6 -> reg rates
plot(pred[[1]][1,,1] , type = "l")

as_array(finn$nnMortEnv(env))

plot(as_array(env)[,30, 1], (as_array(pred[[1]]))[,30,1])
plot(as_array(env)[,30, 1], as_array(finn$nnGrowthEnv(env$to(device = "cuda:0")))[, 30, 3])


self = finn

system.time({

pred = FINN:::predict(dbh = initCohort$dbh,
                    trees = initCohort$trees,
                    species = initCohort$species,
                    response = "BA*T", env = env, patches = patches, debug = FALSE)

})
profvis::profvis({
pred = finn$predict(dbh = initCohort$dbh,
                    trees = initCohort$trees,
                    species = initCohort$species,
                    response = "BA*T", env = env, patches = patches, debug = FALSE)
})

system.time({
pred = finn$predict(dbh = initCohort$dbh,
                    trees = initCohort$trees,
                    species = initCohort$species,
                    response = "BA*T", env = env, patches = patches, debug = FALSE)

})

par(mfrow = c(1,1))
plot(torch::as_array(pred[[1]]$data())[1,,1], type = "l", col = "blue", ylim = c(0, 50))
lines(torch::as_array(pred[[1]]$data())[1,,2], type = "l", col = "green")
lines(torch::as_array(pred[[1]]$data())[1,,3], type = "l", col = "yellow")
lines(torch::as_array(pred[[1]]$data())[1,,4], type = "l", col = "red")

str(pred[[7]])
par(mfrow = c(2,2))
hist(torch::as_array(pred[[2]]$data())[1,,1], type = "l", col = "blue")
hist(torch::as_array(pred[[2]]$data())[1,,2], type = "l", col = "green")
hist(torch::as_array(pred[[2]]$data())[1,,3], type = "l", col = "yellow")
hist(torch::as_array(pred[[2]]$data())[1,,4], type = "l", col = "red")


plot(1:200, sin((1:200)/10), type = "l")

env1 <- torch::as_array(env)[1,,1]

plot(env1+(sin((1:200)/30))*10)
plot()

str(torch::as_array(env))

array(1,)

plot(torch::as_array(env)[1,,1], type = "l", lty = 1, col = "blue")
lines(torch::as_array(env)[1,,2], type = "l", lty = 2, col = "blue")

dim(torch::as_array(env))

dim(torch::as_array(pred[[1]]$data()))

library(torch)

Y = (torch_cat(list(pred[[1]]$unsqueeze(4),
                    pred[[2]]$unsqueeze(4),
                    pred[[3]]$unsqueeze(4),
                    pred[[4]]$unsqueeze(4),
                    pred[[5]]$unsqueeze(4),
                    pred[[6]]$unsqueeze(4),
                    pred[[7]]$unsqueeze(4)), 4))

#Y = Y$unsqueeze(3)
patches = 5L
initCohort2 = CohortMat$new(dims = c(sites, patches, 10),
                            dbh = array(1, dim = c(sites, patches, 10)),
                            trees = array(1, dim = c(sites, patches, 10)),
                            sp = sp)

finn2 = FINN$new(sp = sp, env = 2L, device = "cuda:0", which = "all" ,
                 parGrowth = matrix(c(0.5, 12), sp, 2, byrow = TRUE),
                 parMort = matrix(c(0.0, 2.5), sp, 2, byrow = TRUE),
                 parReg = runif(sp, 0.8, 0.9), # any value between 0 and 1. 0 = species needs no light for regeneration, 1 = species needs full light for regeneration
                 parHeight = runif(sp, 0.5, 0.55), # plausible range is 0 to 1. default should be between 0.3 and 0.7
                 parGrowthEnv = list(matrix(c(0.1, 0.1, 0.1), sp, 2)),
                 parMortEnv = list(matrix(c(0.1, 0.1), sp, 2)),
                 parRegEnv = list(matrix(1, sp, 2)),
                 patch_size_ha = 0.01)
start = lapply(finn2$parameters, as.matrix)
finn2$optimizer = NULL
system.time({
finn2$fit(initCohort = initCohort, X = (env)[,1:60,],Y = Y, patches = patches, batch_size = 50L, epochs = 200L, learning_rate = 0.01,update_step = 1L, weights = NULL)
})

 matplot(sapply(1:length(finn2$param_history), function(i) plogis(finn2$param_history[[i]]$H)) %>% t(), type = "l")

initCohort = initCohort
X = (env)[,1:20,]
Y = Y
patches = patches
batch_size = 100L
epochs = 1L
learning_rate = 1.01
response = "BA*T"
update_step = 5L

dbh = initCohort$dbh
trees = initCohort$trees
species = initCohort$species
x = X
start_time = 1L
response = "ba*T"
y = Y
c =  (Y[,,,2])$round()
update_step=1L
self = finn2
pred_growth = NULL
pred_morth = NULL
pred_reg = NULL
finn2$device = "cpu"
debug = FALSE
source("R/utils.R")
source("R/getter_setter.R")
source("R/distributions.R")
source("R/CohortMat.R")
source("R/FINN.R")
source("R/helpers.R")
finn2$predict( dbh, trees, species, x, start_time = start_time, response = response, y = y, c = c, update_step = update_step)




finn2$parameters$R$grad
finn2$parameters$
start$R

matplot(sapply(1:length(finn2$param_history), function(i) (finn2$param_history[[i]]$H)) %>% t(), type = "l")

profvis::profvis({
  finn2$fit(initCohort = initCohort, X = (env),Y = Y, patches = patches, batch_size = 100L, epochs = 2L, learning_rate = 0.001, response = "BA*T", update_step = 1L)
})


pred2 =
  finn2$predict(initCohort$dbh, initCohort$trees, initCohort$species,response = "ba*T",
                env = env,patches = patches, debug = TRUE)


pred2 = finn2$predict(dbh = initCohort$dbh,
                    trees = initCohort$trees,
                    species = initCohort$species,
                    response = "BA*T", env = env, patches = patches, debug = FALSE)

pred[[1]] - pred2[[1]]

par(mfrow = c(3, 2))
for(i in 1:3) {
plot(torch::as_array(pred2[[2]]$data())[1,,i], type = "l", main = "Fitted")
plot(torch::as_array(pred[[2]]$data())[1,,i], type = "l", main = "True")
}


par(mfrow = c(1, 1))
plot(as.matrix(finn$nnRegEnv(env)[,1,])[order(as.matrix(env[,1,]))])
plot(as.matrix(finn2$nnRegEnv(env)[,1,])[order(as.matrix(env[,1,]))])


plot(as.matrix(finn$nnGrowthEnv(env)[,1,2])[order(as.matrix(env[,1,]))])
plot(as.matrix(finn2$nnGrowthEnv(env)[,1,3])[order(as.matrix(env[,1,]))])


plot(as.matrix(finn$nnMortEnv(env)[,1,2])[order(as.matrix(env[,1,]))])
plot(as.matrix(finn2$nnMortEnv(env)[,1,3])[order(as.matrix(env[,1,]))], ylim = c(0.3, 0.7))


start



finn2$optimizer$zero_grad()

start_time = 1L
y = Y
c = y[,,,2]$round()
loss = torch_zeros(1L)
pred2 =
  finn2$predict(initCohort$dbh, initCohort$trees, initCohort$species,response = "dbh",
                env = env,patches = 10)
for(i in 1:6) {

  # 1 -> dbh/ba
  # 2 -> counts
  # 3 -> AL
  # 4 -> growth rates
  # 5 -> mort rates
  # 6 -> reg rates

  if(i != 2) loss = loss+torch::nnf_mse_loss(y[, start_time:dim(y)[2],,i], pred2[[i]][,start_time:dim(y)[2],])$sum()
  else loss = loss+torch::distr_poisson(pred[[2]][,start_time:pred2[[2]]$shape[2],]+0.001)$log_prob(c[, start_time:c$shape[2],])$sum()$negative()
}


g = torch::autograd_grad(loss, inputs = finn2$parameters[[6]], retain_graph = TRUE,create_graph = TRUE, allow_unused = TRUE)
g[[1]][1]

g2 = torch::autograd_grad(g[[1]][1,1], finn2$parameters[[6]], retain_graph = TRUE,create_graph = TRUE, allow_unused = TRUE)
g22 = torch::autograd_grad(g[[1]][1,2], finn2$parameters[[6]], retain_graph = TRUE,create_graph = TRUE, allow_unused = TRUE)

H = torch_cat(list(g2[[1]], g22[[1]]), dim = 1)
torch_sqrt(torch_diag(torch_inverse(H)))


finn2$parameters[[5]]$grad
finn2$parGrowth$grad



torch::torch_lgamma(0.9)


finn$parGrowth
finn2$parGrowth

finn$parMort
finn2$parMort

finn$parReg
finn2$parReg


finn$parHeight
finn2$parHeight
