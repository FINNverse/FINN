library(dplyr)
library(FINN)
observed =read.csv("/Users/maximilianpichler/Downloads/evaluation2-obs.csv")
data = read.csv("/Users/maximilianpichler/Downloads/evaluation2-env.csv")
m = FINN(~0 + envAct1 + envAct2 ,
         data = data %>% filter(site < 14), epochs = 1L, batch_size = 5L,
         response = observed %>% filter(site < 14))
dim(m$data$env)
dim(m$data$observations)

# env: [sites, timestep, anzahl predictors] -> [100, 100, 4]
# observed: [sites, timestep, species, responses] -> 2x responses, dbh/ba/ba*T und nTRees

# dbh [50, 100, 150, 1] -> [5000, 150, 1]

m = FINN(X = m$data$env[,,1:2], Y = m$data$observations, learning_rate = 0.001,
         data = data %>% filter(site < 14), epochs = 20L, batch_size = 5L, hidden_growth = list(5L, 5L),
         response = observed %>% filter(site < 14), init_global = c(1, 2, 3, 0, 5), init_growth = matrix(1, 5, 2), patches = 2L)
coef(m)
torch = FINN:::pkg.env$torch
pred = m$model$nnGrowthEnv(torch$tensor(m$data$env[,,1:2], dtype=torch$float32))
pred = FINN:::force_r(pred$cpu()$data$numpy())
dim(pred)
plot(pred[,20,1])

pred = predict(m, newdata = NULL, patches = 200L)
dim(pred[[1]])
# learning rate = step size, batch size , und epochs --> stochastic gradient descent
# 100 sites, batch size = 5, -> epoch aus 20 schritten

en = m$data$env
obs = m$data$ba_observations

model = FINN:::pkg.env$FINN$FINN(sp = 5L, env = 2L, device = "cpu")
model$fit(X = en[,,1:2], Y = m$data$ba_observations[,,,], epochs = 10L, response = "ba*T",
          batch_size = 5L, start_time = 0.0, patches = 40L, learning_rate  =0.01)
model$nnGrowthEnv(torch$tensor(en[,,1:2], dtype = torch$float32))
model$parGrowth
model$continue_fit(X = en[,,1:2], Y = m$data$ba_observations[,,,], epochs = 1L, response = "ba*T",
          batch_size = 5L, start_time = 0.7, patches = 40L, learning_rate  =0.001)



model = FINN:::pkg.env$FINN$FINN(sp = 5L, env = 4L, device = "cpu")
model$fit(X = en[38:50,,], Y = m$data$ba_observations[38:50,,,], epochs = 1L, response = "ba*T",
          batch_size = 5L, start_time = 0.7, patches = 40L, learning_rate  =0.01)
model$continue_fit(X = en[38:50,,], Y = m$data$ba_observations[38:50,,,], epochs = 15L, response = "ba*T",
                   batch_size = 5L, start_time = 0.7, patches = 40L, learning_rate =0.01)

pred = model$predict(env = en[38:50,,], patches = 2L )
dbh = FINN:::force_r(pred[1]$cpu()$data$numpy())



FINN:::.onLoad()

torch = FINN:::pkg.env$torch
# dbh [sites, patches, cohorts, 1]

cohort = FINN:::pkg.env$FINN$CohortMat(dbh = array(1.0, dim = c(13, 40, 5, 1)),
                              nTree = array(1, dim = c(13, 40, 5, 1) ),
                              Species = array(matrix(0:4, 40, 5, byrow = TRUE), dim = c(13, 40, 5))
                              )


model = FINN:::pkg.env$FINN$FINN(sp = 5L, env = 2L, device = "cpu")
model$fit(X = m$data$env[,,1:2], Y = m$data$ba_observations[,,,], epochs = 10L, response = "ba*T", cohorts = cohort,
          batch_size = 5L, start_time = 0.0, patches = 40L, learning_rate  =0.01)
