################################################################################
## simulate from model with predefined parameters
################################################################################
library(FINN)
library(data.table)

Nsites = 16
Ntimesteps = 500
Npatches = 10
NsamplesPerEnv = ceiling(sqrt(Nsites))
EnvSiteAvg = seq(0,1,length.out = NsamplesPerEnv)
EnvSiteAvgM = expand.grid(list(EnvSiteAvg,EnvSiteAvg))

cat("start at")
Sys.time()
envM_list = list()
env_out_all  <- data.table()
observed  <- data.table()
for(i_site in 1:Nsites){
  envAvg1 = EnvSiteAvgM[i_site,1]
  envAvg2 = EnvSiteAvgM[i_site,2]
  envM = matrix(c(rnorm(Ntimesteps, envAvg1, 0.1), rnorm(Ntimesteps, envAvg2, 0.1)), Ntimesteps, 2)
  # envM[envM > 1] = 1
  # envM[envM < 0] = 0
  envM_list[[i_site]] = envM
  env_out = data.table(
    timestep = 1:Ntimesteps,
    site = i_site,
    envAct1 = envM[,1],
    envAct2 = envM[,2],
    envAvg1 = envAvg1,
    envAvg2 = envAvg2
  )
  env_out_all = rbind(env_out_all, env_out)
  observed = rbind(observed, data.table(Species = rep(0:1, each = Ntimesteps), dbh = 0, nTree = 0, timestep = 1:Ntimesteps, site = i_site))
}

#array(cbind(c(0,1),c(1,1)), dim = c())
m2 =
FINN::FINN(
  formula = ~0+envAct1+envAct2,patches = Npatches,
  # X = m$data$env[,,1:2],
  # Y = m$data$ba_observations[,,,],
  data = env_out_all,
  response = observed,
  init_global = c(0.5,0.7),
  init_growth =
    rbind(c(1,3),
          c(5,5)),
  init_mort =
    rbind(c(1,40),
          c(.1,10)),
  init_reg = c(0.1,0.5),
  hidden_growth = list(),
  hidden_mort = list(),
  hidden_reg = list(),
  init_growth_env = list(array(rbind(c(0.1,.1),
                                     c(0.1,.1)),dim = c(2,2))),
  init_mort_env = list(array(.1,dim = c(2,2))),
  init_reg_env = list(array(rbind(c(0.01,0.01),
                                  c(0.04,0.04)),dim = c(2,2))),
  simulate = T
  )

# to_Tensor(parGlobal, dtype = dtype, device = dtype = "float32", FALSE, correct_zero = TRUE)
Sys.time()
pred = FINN:::predict.FINN(object = m2, newdata = m2$data$env[,,1:2])
Sys.time()
str(pred)
my_array = pred[[2]]

# Generate the indices for each dimension
dim_indices <- expand.grid(site = 1:dim(my_array)[1],
                           timestep = 1:dim(my_array)[2],
                           Species = 1:dim(my_array)[3])

# Convert the 3D array to a vector and combine with indices
long_data <- cbind(dim_indices, ba = as.vector(my_array))

# Convert to data.table
dt_long <- as.data.table(long_data)


out_dt_all_pred <- merge(dt_long, env_out_all, by = c("site", "timestep"))
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

# out_dt_all_pred[,SpeciesName := factor(Species, levels = 1:5, labels = paste(1:5, species_dt$names)),]

p2 <- ggplot(out_dt_all_pred[,.(
  ba = (mean(ba))
), by = .(site,timestep,Species,envAvg1,envAvg2)], aes(x = timestep, y = ba, color = factor(Species)))+
  geom_line()+
  facet_grid(round(envAvg1,2)~round(envAvg2,2))+
  coord_cartesian()
p2



fit1 = FINN::FINN(
  X = m2$data$env,
  Y = m2$data$observations, epochs = 4, batch_size = 16,
  init_global = c(0.5,0.7),
  init_growth =
    rbind(c(1,3),
          c(5,5)),
  init_mort =
    rbind(c(1,40),
          c(.1,10)),
  init_reg = c(0.1,0.5),
  hidden_growth = list(),
  hidden_mort = list(),
  hidden_reg = list(),
  learning_rate = 1, patches = 10, response = "ba*T")

Sys.time()
pred2 = FINN:::predict.FINN(object = fit1, newdata = fit1$data$env[,,1:2])
Sys.time()
my_array = pred2[[2]]

# Generate the indices for each dimension
dim_indices <- expand.grid(site = 1:dim(my_array)[1],
                           timestep = 1:dim(my_array)[2],
                           Species = 1:dim(my_array)[3])

# Convert the 3D array to a vector and combine with indices
long_data <- cbind(dim_indices, ba = as.vector(my_array))

# Convert to data.table
dt_long <- as.data.table(long_data)


out_dt_all_pred <- merge(dt_long, env_out_all, by = c("site", "timestep"))
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

# out_dt_all_pred[,SpeciesName := factor(Species, levels = 1:5, labels = paste(1:5, species_dt$names)),]

p2 <- ggplot(out_dt_all_pred[,.(
  ba = (mean(ba))
), by = .(site,timestep,Species,envAvg1,envAvg2)], aes(x = timestep, y = ba, color = factor(Species)))+
  geom_line()+
  facet_grid(round(envAvg1,2)~round(envAvg2,2))+
  coord_cartesian()
p2

