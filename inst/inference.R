library(FINN)
initCohort = CohortMat$new(dims = c(100, 1, 10),
                           dbh = array(1, dim = c(100, 1, 10)),
                           trees = array(1, dim = c(100, 1, 10)),
                           sp = 2)

finn = FINN$new(sp = 2L, env = 2L, device = "cpu", which = "all" ,
                parGrowth = matrix(c(-10, 20), 2, 2, byrow = TRUE),
                parMort = matrix(c(-3, 100), 2, 2, byrow = TRUE),
                parReg = c(-5, -3),
                parHeight = c(5, 2),
                parGrowthEnv = list(matrix(5.0, 2, 2)),
                parMortEnv = list(matrix(c(0, 0), 2, 2)),
                parRegEnv = list(matrix(5, 2, 2)))
env = torch::torch_randn(size = c(100, 200, 2))
pred =
  finn$predict(initCohort$dbh, initCohort$trees, initCohort$species,response = "BA*T",
               env = env, patches = 1)
plot(torch::as_array(pred[[1]]$data())[1,,1], type = "l")

library(torch)
Y = (torch_cat(list(pred[[1]]$unsqueeze(4),
                    pred[[2]]$unsqueeze(4),
                    pred[[3]]$unsqueeze(4),
                    pred[[4]]$unsqueeze(4),
                    pred[[5]]$unsqueeze(4),
                    pred[[6]]$unsqueeze(4)), 4))

#Y = Y$unsqueeze(3)
initCohort2 = CohortMat$new(dims = c(100, 20L, 10),
                            dbh = array(1, dim = c(100, 20L, 10)),
                            trees = array(1, dim = c(100, 20L, 10)),
                            sp = 2)
finn2 = FINN$new(sp = 2L, env = 2L, device = "cpu", which = "all"
                 #hidden_growth = c(10L, 10L),
                 #hidden_mort = c(10L, 10L),
                 #hidden_reg = c(10L, 10L),  # all parameters, env, or species
                 # parGrowth = matrix(c(-10, 10), 1, 2),
                 # parMort = cbind(0, 100, 0),
                 # parReg = -5,
                 # parHeight = 5
                 # parGrowthEnv = list(matrix(5.0, 1, 2)),
                 # parMortEnv = list(matrix(-5, 1, 2)),
                 # parRegEnv = list(matrix(5, 1, 2))
)
start = lapply(finn2$parameters, as.matrix)
finn2$fit(initCohort = initCohort, X = (env),Y = Y, patches = 20L, batch_size = 100L, epochs = 200L, learning_rate = 1.0, response = "BA*T")


pred2 =
  finn2$predict(initCohort$dbh, initCohort$trees, initCohort$species,response = "ba*T",
                env = env,patches = 1L)
par(mfrow = c(1, 2))
plot(torch::as_array(pred2[[1]]$data())[1,,1], type = "l")
plot(torch::as_array(pred[[1]]$data())[1,,1], type = "l")



par(mfrow = c(3, 2))
plot(as.matrix(finn$nnRegEnv(env)[,1,])[order(as.matrix(env[,2,]))])
plot(as.matrix(finn2$nnRegEnv(env)[,1,])[order(as.matrix(env[,2,]))])


plot(as.matrix(finn$nnGrowthEnv(env)[,1,])[order(as.matrix(env[,1,]))])
plot(as.matrix(finn2$nnGrowthEnv(env)[,1,])[order(as.matrix(env[,1,]))])


plot(as.matrix(finn$nnMortEnv(env)[,1,])[order(as.matrix(env[,1,]))])
plot(as.matrix(finn2$nnMortEnv(env)[,1,])[order(as.matrix(env[,1,]))], ylim = c(0.3, 0.7))


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
