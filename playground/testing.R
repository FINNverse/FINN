library(dplyr)
library(FINN)
library(data.table)
library(ggplot2)
observed = fread("data/testing-data/evaluation2-obs.csv")
data = fread("data/testing-data/evaluation2-env.csv")
m = FINN(~0 + envAct1 + envAct2 ,
         data = data %>% filter(site < 101), batch_size = 20, start_time=0.0, epochs=3, learning_rate=0.1,
         response = observed %>% filter(site < 101))
dim(m$data$env)
dim(m$data$observations)

pred = FINN:::predict.FINN(m)
my_array = pred[[2]]

# Generate the indices for each dimension
dim_indices <- expand.grid(site = 1:dim(my_array)[1],
                           timestep = 1:dim(my_array)[2],
                           Species = 1:dim(my_array)[3])

# Convert the 3D array to a vector and combine with indices
long_data <- cbind(dim_indices, ba = as.vector(my_array))

# Convert to data.table
dt_long <- as.data.table(long_data)


out_dt_all_pred <- merge(dt_long, data, by = c("site", "timestep"))
# out_dt_all_pred[,ba := (circle_area(dbh)*nTree)/0.1,]
species_dt = data.table(
  names = c("Fagus sylvatica", "Picea abies", "Acer pseudoplatanus", "Betula pendula", "Abies alba"),
  dbh2height = c(0.6,0.6,0.6,0.6,0.6),
  Light = c(0.1,0.15,0.5,0.6,0.07),
  Temp = c(0.6,0.3,0.7,0.4,0.3),
  Prec = c(0.6,0.7,0.4,0.4,0.5),
  maxD = c(350,300,200,100,350),
  maxG = c(3,5,7,8,4)
)

out_dt_all_pred[,SpeciesName := factor(Species, levels = 1:5, labels = paste(1:5, species_dt$names)),]

p2 <- ggplot(out_dt_all_pred[,.(
  ba = log(mean(ba))
), by = .(site,timestep,SpeciesName,envAvg1,envAvg2)], aes(x = timestep, y = ba, color = SpeciesName))+
  geom_line()+
  facet_grid(round(envAvg1,2)~round(envAvg2,2))
p2


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
