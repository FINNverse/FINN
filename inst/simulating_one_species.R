library(FINN)
FINN = FINN:::pkg.env$FINN
torch = FINN:::pkg.env$torch
as_ten = function(x, dtype=torch$float32) {torch$tensor(x, dtype=dtype)}

to_r = function(x) FINN:::force_r(x$cpu()$data$numpy())

env = matrix(runif(100*1,-2,2), 100, 1)
env = cbind(1, env)
env_a = abind::abind(lapply(1:50, function(i) env), along = 0)
dim(env_a)
env_a = aperm(env_a, perm = c(2, 1, 3))
dim(env_a)


model= FINN$FINN(sp = 1L, env = 2L,
          parGlobal = 0.8,
          parGrowth = matrix(c(0.5, 0.5), 1, 2),
          parMort = matrix(c(1.0, 10.0), 1, 2),
          parReg = 0.2,
          parGrowthEnv = list(matrix(c(0,-10), 1, 2)),
          parMortEnv = list(matrix(c(0,-10), 1, 2)),
          parRegEnv = list(matrix(c(3,-10), 1, 2))
          )
plot(env_a[,1,2], to_r(model$nnMortEnv(as_ten(env_a[,1,]**2))) )
plot(env_a[,1,2], to_r(model$nnGrowthEnv(as_ten(env_a[,1,]**2))) )
plot(env_a[,1,2], to_r(model$nnRegEnv(as_ten(env_a[,1,]**2))) )

pred = model$predict(env = as_ten(env_a**2))
dbh = to_r(pred[[0]])
nTree = to_r(pred[[1]])
plot(dbh[,50,1])
plot(nTree[,50,1])

plot(dbh[10,,1])


