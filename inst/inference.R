initCohort = CohortMat$new(dims = c(100, 10, 10),
                           dbh = array(1, dim = c(100, 10, 10)),
                           trees = array(1, dim = c(100, 10, 10)),
                           sp = 1)

finn = FINN$new(sp = 1L, env = 2L, device = "cpu", which = "all" ,
                parGrowth = matrix(c(-10, 10), 1, 2),
                parMort = cbind(0, 100, 0),
                parReg = -5,
                parHeight = 5,
                parGrowthEnv = list(matrix(5.0, 1, 2)),
                parMortEnv = list(matrix(-5, 1, 2)),
                parRegEnv = list(matrix(5, 1, 2)))
env = torch::torch_randn(size = c(100, 200, 2))
pred =
  finn$predict(initCohort$dbh, initCohort$trees, initCohort$species,response = "dbh",
               env = env,patches = 10)
plot(torch::as_array(pred[[2]]$data())[1,,1], type = "l")


Y = (torch_cat(list(pred[[1]], pred[[2]], pred[[3]], pred[[4]], pred[[5]], pred[[6]]), 3))

Y = Y$unsqueeze(3)
initCohort2 = CohortMat$new(dims = c(100, 20L, 10),
                            dbh = array(1, dim = c(100, 20L, 10)),
                            trees = array(1, dim = c(100, 20L, 10)),
                            sp = 1)
finn2 = FINN$new(sp = 1L, env = 2L, device = "cpu", which = "all", # all parameters, env, or species
                 # parGrowth = matrix(c(-10, 10), 1, 2),
                 # parMort = cbind(0, 100, 0),
                 # parReg = -5,
                 # parHeight = 5
                 # parGrowthEnv = list(matrix(5.0, 1, 2)),
                 # parMortEnv = list(matrix(-5, 1, 2)),
                 # parRegEnv = list(matrix(5, 1, 2))
)
start = lapply(finn2$parameters, as.matrix)
finn2$fit(initCohort = initCohort, X = (env),Y = Y, patches = 20L, batch_size = 100L, epochs = 200L, learning_rate = 1.1)


pred2 =
  finn2$predict(initCohort$dbh, initCohort$trees, initCohort$species,response = "dbh",
                env = env,patches = 10)
par(mfrow = c(1, 2))
plot(torch::as_array(pred2[[1]]$data())[1,,1], type = "l")
plot(torch::as_array(pred[[1]]$data())[1,,1], type = "l")



par(mfrow = c(3, 2))
plot(as.matrix(finn$nnRegEnv(env)[,1,])[order(as.matrix(env[,1,]))])
plot(as.matrix(finn2$nnRegEnv(env)[,1,])[order(as.matrix(env[,1,]))])


plot(as.matrix(finn$nnGrowthEnv(env)[,1,])[order(as.matrix(env[,1,]))])
plot(as.matrix(finn2$nnGrowthEnv(env)[,1,])[order(as.matrix(env[,1,]))])


plot(as.matrix(finn$nnMortEnv(env)[,1,])[order(as.matrix(env[,1,]))])
plot(as.matrix(finn2$nnMortEnv(env)[,1,])[order(as.matrix(env[,1,]))])



finn$parGrowth
finn2$parGrowth

finn$parMort
finn2$parMort

finn$parReg
finn2$parReg


finn$parHeight
finn2$parHeight
