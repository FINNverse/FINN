library(data.table)
library(FINN)
library(torch)

df = fread("data/calibration-data/BCI-1h-patch/stand_dt.csv")
obs_df = fread("data/calibration-data/BCI-1h-patch/obs_df.csv")
env = fread("data/calibration-data/BCI-1h-patch/env_dt.csv")

# df = df[census >= 1985]
# obs_df = obs_df[year >= 1985]
# env = env[year >= 1985]

head(df)
df
df$AL = NA
colnames(df)[c(1,7, 9, 11)] = c("year", "growth", "reg", "mort")
head(df)
data = df

# fill sites
# Jedes jahr alle sites!
empty_res =
  lapply(unique(data$species), function(sp) {
    empty_sites =
      lapply(unique(data$year), function(yy) {
          expand.grid(siteID = setdiff(data$siteID |> unique(), data[species==sp & year == yy]$siteID),
                      year = yy,
                      species = sp
          )
      })
    rbindlist(empty_sites, fill = TRUE)
  })
data2 = rbindlist(list(data, rbindlist(empty_res, fill = TRUE)), fill = TRUE)
data = data2

dist_dt <- unique(data[,.(year,siteID)])

#
# dist_dt[,intensity := rbinom(.N, 1, runif(1,0.0043, 0.016)), by = year]
# Env prep.


# missing_years =
# lapply(1982:1984, function(year) {
#   tmp = env[2,]
#   tmp$year = year
#   return(tmp)
#   }) |> rbindlist()



# env = rbindlist(list(missing_years, env[-1,]))
env = env[!year %in% 2016:2019]

env =
  lapply(unique(data$siteID), function(site) {
    tmp = env
    tmp$siteID = site
    return(tmp)
  }) |> rbindlist()

env$Prec = scale(env$Prec)
env$T_mean = scale(env$T_mean)
env$T_max = scale(env$T_max)
env$T_min = scale(env$T_min)
env$SR_kW_m2 = scale(env$SR_kW_m2)
env$RH_prc = scale(env$RH_prc)
data$AL = as.numeric(data$AL)
str(data)

# Preparation of cohorts! BUG!!!!
cohort1 <- FINN::CohortMat$new(obs_df)
sp= cohort1$species_r
sp[sp==0] = 1L
sp[is.na(sp)] = 1L
cohort2 = FINN::CohortMat$new(dbh = cohort1$dbh_r, trees = cohort1$trees_r, species = sp)
# sp zum indexen, aber wenn 0 exisitiert -> seg fault / index fault


m = finn(env = env, data = data,
         patches = 1L,
         mortalityProcess = createProcess(~., func = mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
         growthProcess = createProcess(~., func = growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
         regenerationProcess = createProcess(~., func = regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
         device = "gpu",
         optimizeHeight = TRUE,
         init = cohort2,
         lr=0.01,
         epochs = 100L,
         patch_size = 1,
         batchsize = 350L,
         weights = c(1, 0.1, 3, 1.5, 3, 1.0))


m2 = finn(env = env, data = data,
          mortalityProcess = createProcess(~., func = mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE, initSpecies = m$model$parMortTR, initEnv = m$model$parMortEnv_r),
          growthProcess = createProcess(~., func = growth, optimizeSpecies = TRUE, optimizeEnv = TRUE, initSpecies = m$model$parGrowthTR, initEnv = m$model$parGrowthEnv_r),
          regenerationProcess = createProcess(~., func = regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE, initSpecies = m$model$parRegTR, initEnv = m$model$parRegEnv_r),
          device = "gpu",
          height = as.vector(m$model$parHeightTR),
          optimizeHeight = TRUE,
          init = cohort2,
          lr=0.01,
          epochs = 1000L,
          patch_size = 1,
          batchsize = 350L,
          weights = c(1, 0.1, 3, 1.5, 3, 1.0),
          file = "BCI_parameters_1ha.RDS")

m = readRDS(file = "gfoe2024/bci_model.RDS")
m$model$check(device = "cpu")
m$model$device = "cpu"
m$model$device_r = "cpu"

pred = predict(m)

library(ggplot2)
pred$species = as.factor(pred$species)
ggplot(pred[siteID == 25 & (species %in% 120:195)], aes(x = year, y = value,  color = (species))) +
  geom_line() +
  labs(x = "Year",
       y = "value") +
  theme_minimal()+
  facet_wrap(~variable, scales = "free_y")





