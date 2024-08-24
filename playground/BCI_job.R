library(data.table)
library(FINN)
library(torch)

df = fread("~/stand_dt.csv")
obs_df = fread("~/obs_df.csv")
env = fread("~/env_dt.csv")

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

# Env prep.


missing_years =
  lapply(1982:1984, function(year) {
    tmp = env[2,]
    tmp$year = year
    return(tmp)
  }) |> rbindlist()



env = rbindlist(list(missing_years, env[-1,]))
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




# Option A) Train model
m = finn(env = env, data = data,
         patches = 1L,
         mortalityProcess = createProcess(~., func = mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
         growthProcess = createProcess(~., func = growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
         regenerationProcess = createProcess(~., func = regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
         device = "gpu",
         optimizeHeight = TRUE,
         init = cohort2,
         lr=0.01,
         epochs = 13000L,
         patch_size = 0.1,
         batchsize = 350L,
         weights = c(1, 0.1, 3, 1.5, 3, 1.0))

# Check for convergence!!!
# matplot(sapply(1:5000, function(i) m$model$param_history[[i]]$parGrowth[1:10, 2]) |> t(), type = "l")


# Option B) Continue training!

saveRDS(m, file = "model.RDS")

# m = readRDS(file = "gfoe2024/bci_model.RDS")
# m$model$check()
#
# m2 = finn(env = env, data = data,
#           mortalityProcess = createProcess(~., func = mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE, initSpecies = m$model$parMortTR, initEnv = m$model$parMortEnv_r),
#           growthProcess = createProcess(~., func = growth, optimizeSpecies = TRUE, optimizeEnv = TRUE, initSpecies = m$model$parGrowthTR, initEnv = m$model$parGrowthEnv_r),
#           regenerationProcess = createProcess(~., func = regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE, initSpecies = m$model$parRegTR, initEnv = m$model$parRegEnv_r),
#           device = "gpu",
#           height = as.vector(m$model$parHeightTR),
#           optimizeHeight = TRUE,
#           init = cohort2,
#           lr=0.01,
#           epochs = 3000L,
#           patch_size = 0.1,
#           batchsize = 350L,
#           weights = c(1, 0.1, 3, 1.5, 3, 1.0),
#           file = "BCI_parameters.RDS")
#
# saveRDS(m2, file = "continue_model.RDS")
#
# #
# # pred = predict(m)
# # library(ggplot2)
# # pred = pred$long$site
# #
# # ggplot(pred[siteID == 25 & (species %in% 1:100)], aes(x = year, y = value,  color = factor(species))) +
# #   geom_line() +
# #   labs(x = "Year",
# #        y = "value") +
# #   theme_minimal()+
# #   facet_wrap(~variable, scales = "free_y")
# # #saveRDS(m, file = "gfoe2024/bci_model.RDS")
# #
# #
# #
# # fields::image.plot(m$model$parRegEnv_r[[1]] |> t())
# #
# # apply(m$model$parMortEnv_r[[1]], 2, mean)[-1]
# # #
#
#
# # Simulate with fitted parameters
# # library(FINN)
# # library(data.table)
# # getPars = function(internalPar, parRange) {
# #   internalPar = torch::torch_tensor(internalPar)
# #   if (is.vector(parRange)) {
# #     # Case where internalPar is a 1D tensor and parRange is a vector
# #     Npar <- length(parRange) / 2
# #     lower <- parRange[1:Npar]
# #     upper <- parRange[(Npar + 1):(2 * Npar)]
# #
# #     out <- internalPar$sigmoid() * (upper - lower) + lower
# #   } else {
# #     # Case where internalPar is a matrix and parRange is a matrix
# #     Npar <- ncol(internalPar)
# #     out <- list()
# #     for (i in 1:Npar) {
# #       lower <- parRange[i, 1, drop = FALSE]
# #       upper <- parRange[i, 2, drop = FALSE]
# #       out[[i]] <- internalPar[, i, drop = FALSE]$sigmoid() * (upper - lower) + lower
# #     }
# #     out <- torch::torch_cat(out, dim = 2L)
# #   }
# #   return(out |> as.matrix())
# # }
# #
# # w = readRDS("playground/BCI_parameters.RDS") # raw untransformierten Parameter
# # matplot(sapply(1:length(w), function(i) w[[i]]$parHeight[11:50]) |> t() |> plogis(), type = "l")
# #
# #
# # pars = w[[length(w)]]
# #
# # speciesPars_ranges = list(parGrowth = rbind(c(0.01, 0.99), c(0.01, 4)), parMort =
# #                             rbind(c(0.01, 0.99), c(0, 4)), parReg = c(0.01, 0.99), parHeight = c(0.3, 0.7),
# #                           parGrowthEnv = rbind(c(-1, 1), c(-1, 1)), parMortEnv = rbind(c(-2, 2), c(-2, 2)),
# #                           parRegEnv = rbind(c(-2, 2), c(-2, 2)))
# #
# # # laeuft nicht...
# # simulation = simulateForest(env = env,
# #                     init = cohort2,
# #                     sp = 195L,
# #                     mortalityProcess = createProcess(~., func = mortality, optimizeSpecies = FALSE, optimizeEnv = FALSE, initSpecies = getPars(pars$parMort, speciesPars_ranges$parMort), initEnv = list(pars$nnMort.0.weight)),
# #                     growthProcess = createProcess(~., func = growth, optimizeSpecies = FALSE, optimizeEnv = FALSE, initSpecies = getPars(pars$parGrowth, speciesPars_ranges$parGrowth), initEnv = list(pars$nnGrowth.0.weight)),
# #                     regenerationProcess = createProcess(~., func = regeneration, optimizeSpecies = FALSE, optimizeEnv = FALSE, initSpecies = getPars(pars$parReg, speciesPars_ranges$parReg)[,1], initEnv = list(pars$nnReg.0.weight)),
# #                     device = "cpu",
# #                     height = getPars(pars$parHeight, speciesPars_ranges$parHeight)[,1])
# #
#
#
#
#
#
