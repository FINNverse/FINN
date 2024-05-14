library(FINN)
FINN:::.onLoad()
FINN = FINN:::pkg.env$FINN
torch = FINN:::pkg.env$torch
as_ten = function(x, dtype=torch$float32) {torch$tensor(x, dtype=dtype)}
to_r = function(x) FINN:::force_r(x$cpu()$data$numpy())

env = matrix(runif(30*1,-2,2), 30, 1)
env = cbind(1, env)
env_a = abind::abind(lapply(1:100, function(i) env), along = 0)
dim(env_a)
env_a = aperm(env_a, perm = c(2, 1, 3))
dim(env_a)


model= FINN$FINN(sp = 2L, env = 2L,
                 parGlobal = rep(0.5, 2),
                 parGrowth = matrix(c(2, 5), 2, 2, byrow = TRUE),
                 parMort = matrix(c(1.0, 2.1), 2, 2, byrow = TRUE),
                 parReg = rep(-5.1, 2),
                 parGrowthEnv = list(matrix(c(0,2), 2, 2, byrow = TRUE)),
                 parMortEnv = list(matrix(c(-5,0), 2, 2, byrow = TRUE)),
                 parRegEnv = list(matrix(c(5.3, 0.0), 2, 2, byrow = TRUE))
)
plot(env_a[,1,2], to_r(model$nnGrowthEnv(as_ten(env_a[,1,])))[,1], ylim = c(0,1) )
plot(env_a[,1,2], to_r(model$nnRegEnv(as_ten(env_a[,1,])))[,1], ylim = c(0,1) )

pred = model$predict(env = as_ten(env_a), response = "dbh", patches = 15L)
dbh = to_r(pred[[0]])
nTree = to_r(pred[[1]])
plot(env_a[,1,2], dbh[,30,1]/nTree[,30,1])


plot(1:100, nTree[20,,2])

plot(1:50, dbh[which.max(env_a[,1,2]),,1]/nTree[which.max(env_a[,1,2]),,1], xlab = "time", ylab = "dbh/nTree")
points(1:50, dbh[which.min(env_a[,1,2]),,1]/nTree[which.min(env_a[,1,2]),,1], col = "red")

plot(env_a[,1,2], nTree[,30,1])

plot(1:50, nTree[which.max(env_a[,1,2]),,1], xlab = "time", ylab = "dbh/nTree")
points(1:50, nTree[which.min(env_a[,1,2]),,1], col = "red")



patches =
lapply(1:100, function(i) {
  pred = model$predict(env = as_ten(env_a), response = "ba*T", patches = 1L)
  dbh = to_r(pred[[0]])
  nTree = to_r(pred[[1]])
  return(list(dbh, nTree))
})
pred100 = model$predict(env = as_ten(env_a), response = "ba*T", patches = 100L)
dbh100 = to_r(pred100[[0]])


dbh1 = abind::abind(lapply(patches, function(i) i[[1]]), along = 3)
nTree1 = abind::abind(lapply(patches, function(i) i[[2]]), along = 3)
dbh1 = apply(dbh1, 1:2, mean)
nTree1 = apply(nTree1, 1:2, mean)
plot(env_a[,1,2], dbh1[,25])
plot(env_a[,1,2], dbh100[,25,1])
plot(env_a[,1,2], nTree1[,25])


pred = model$predict(env = as_ten(env_a), response = "ba*T", patches = 10L)
dbh = to_r(pred[[0]])
nTree = to_r(pred[[1]])
plot(env_a[,1,2], dbh[,30,1]/nTree[,30,1])

plot(env_a[,1,2], 1-to_r(model$nnMortEnv(as_ten(env_a[,1,]))) )
plot(env_a[,1,2], to_r(model$nnRegEnv(as_ten(env_a[,1,]))) )

growth(rep(30, 1), 0,rep(10, 1), rep(1., 1), matrix(c(1, 5), 1, 2), 0.0)
mort(rep(500000, 1), 0,rep(100, 1), rep(0.5, 1), matrix(c(0.1, 1), 1, 2), 1.0)



pred = model$predict(env = as_ten(env_a), response = "ba*T", patches = 1L)
dbh = to_r(pred[[0]])
nTree = to_r(pred[[1]])
plot(env_a[,1,2], dbh[,30,1]/nTree[,30,1])

plot(env_a[,1,2], nTree[,30,1])
plot(env_a[,1,2], apply(dbh, c(1, 3), mean)/apply(nTree, c(1, 3), mean))

nTree = to_r(pred[[1]])
plot(env_a[,1,2], dbh[,30,1])
plot(env_a[,1,2], nTree[,30,1])

plot(dbh[,,1])
plot(apply(nTree, c(1, 3), max))


ba_T_pred = to_r( pred[[0]] )
nTree_pred = to_r(pred[[1]] )
obs = abind::abind(ba_T_pred, nTree_pred, along = 3)
obs = to_r(as_ten(obs)$unsqueeze(2L))
fit= FINN$FINN(sp = 1L,
               which = "growth",
               env = 2L,
               parGlobal = 0.5,
               parGrowth = matrix(c(1, 2), 1, 2),
               parMort = matrix(c(1.0, 2.1), 1, 2),
               parReg = 0.1,
               parGrowthEnv = list(matrix(c(0,0), 1, 2)),
               parMortEnv = list(matrix(c(0,2), 1, 2)),
               parRegEnv = list(matrix(c(0.3, 0.0), 1, 2)))
fit$optimizer = reticulate::py_none()
fit$parameters
fit$fit(env_a, Y = obs, learning_rate = 0.01, epochs = 50L, batch_size = 50L, response = "ba*T", patches = 10L)

fit$`_parGlobal`
fit$RegEnv
fit$GrowthEnv
fit$MortEnv

growths = seq(-20, 20, length.out = 20)
growth_opt = sapply(growths, function(g) {
  g_m = list(matrix(c(0,g), 1, 2))
  fit$set_weights_nnGrowthEnv(g_m)
  fit$GrowthEnv
  p = fit$predict(env = as_ten(env_a), response = "dbh", patches = 10L)
  dbh_pred = to_r(p[[0]])
  nTree_pred = to_r(p[[1]])
  return(c(mean((dbh_pred-dbh)**2), mean((nTree_pred-nTree)**2)))

})
plot(growths, ((growth_opt[1,]/100/2)*pi)**2*(growth_opt[2,]))
plot(growths, log(growth_opt[2,]))


growths = expand.grid(seq(-5, 5, length.out = 10), seq(-5, 5, length.out = 10))
growth_opt = sapply(1:nrow(growths), function(g) {
  g_m = list(matrix(c(growths[g,1],growths[g,2]), 1, 2))
  fit$set_weights_nnGrowthEnv(g_m)
  fit$GrowthEnv
  p = fit$predict(env = as_ten(env_a**2), response = "dbh", patches = 10L)
  dbh_pred = to_r(p[[0]])
  nTree_pred = to_r(p[[1]])
  return(c(mean((dbh_pred-dbh)**2), mean((nTree_pred-nTree)**2)))

})

image(matrix(growth_opt[1,], 10, 10, byrow = TRUE))

plot(growths, (growth_opt[1,])*(growth_opt[2,]))

g_m = list(matrix(c(-2,1), 1, 2))
fit$set_weights_nnGrowthEnv(g_m)
mort = seq(-20, 20, length.out = 40)
mort_opt = sapply(mort, function(g) {
  m_m = list(matrix(c(-2,g), 1, 2))
  fit$set_weights_nnMortEnv(m_m)
  fit$MortEnv
  pred = fit$predict(env = as_ten(env_a**2), response = "dbh", patches = 20L)
  dbh_pred = to_r(pred[[0]])
  nTree_pred = to_r(pred[[1]])
  return(c(mean((dbh_pred-dbh)**2), mean((nTree_pred-nTree)**2)))

})
plot(mort, mort_opt[2,]*mort_opt[1,])


m_m = list(matrix(c(-2,1), 1, 2))
fit$set_weights_nnMortEnv(g_m)
r = seq(-5, 5, length.out = 20)
r_opt = sapply(r, function(g) {
  r_m = list(matrix(c(0.1,g), 1, 2))
  fit$set_weights_nnRegEnv(r_m)
  pred = fit$predict(env = as_ten(env_a**2), response = "dbh", patches = 10L)
  dbh_pred = to_r(pred[[0]])
  nTree_pred = to_r(pred[[1]])
  return(c(mean((dbh_pred-dbh)**2), mean((nTree_pred-nTree)**2)))

})
plot(r, log(r_opt[1,]*r_opt[2,]))


