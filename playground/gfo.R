library(FINN)
library(data.table)

# one species, 100 patches, 50 sites, 1 env

Nsp = 5


speciesPars_ranges = list(
  parGrowth = rbind(
    c(0.01, 0.99),
    c(1, 4)
  ),
  parMort = rbind(
    c(0.01, 0.99),
    c(1, 4)
  ),
  parReg = c(0.01, 0.99),
  parHeight = c(0.3, 0.7),
  parGrowthEnv = rbind(
    c(0, 2),
    c(-2, 2)
  ),
  parMortEnv = rbind(
    c(0, 2),
    c(-2, 2)
  ),
  parRegEnv = rbind(
    c(0, 2),
    c(-2, 2)
  ))

speciesPars = list(
  speciesID = 1:Nsp,
  parGrowth = matrix(c(
    runif(Nsp, 0.01, 0.99),
    runif(Nsp, 1, 4)
    ),Nsp, 2),
  parMort = matrix(c(
    runif(Nsp, 0.01, 0.99),
    runif(Nsp, 1, 4)
    ),Nsp, 2),
  parReg = runif(Nsp, 0.01, 0.99),
  parHeight = runif(Nsp, 0.3, 0.7),
  parGrowthEnv = list(matrix(c(
    runif(Nsp, 0, 2),
    runif(Nsp, -2, 2)
    ),Nsp, 2)),
  parMortEnv = list(matrix(c(
    runif(Nsp, 0, 2),
    runif(Nsp, -2, 2)
    ), Nsp, 2)),
  parRegEnv = list(matrix(c(
    runif(Nsp, 0, 2),
    runif(Nsp, -2, 2)
    ),Nsp, 2))
  )

# # speciesPars_ranges$parMort[,1,drop=FALSE]$sigmoid()]
#
# parRange = speciesPars_ranges$parMort
# parRange = speciesPars_ranges$parReg
# inputPar = speciesPars$parReg
# # inputPar = speciesPars$parMortEnv
# # checkPars = function(inputPar, parRange){
# #   if(is.list(inputPar)) inputPar = inputPar[[1]]
# #   if(is.vector(inputPar)) {
# #     Npar = 1
# #     Nsp = length(inputPar)
# #     inputPar = matrix(inputPar, nrow = Nsp, ncol = Npar)
# #     parRange = matrix(parRange, nrow = Npar, ncol = 2)
# #   }else if(is.matrix(inputPar)){
# #     Npar = ncol(inputPar)
# #     Nsp = nrow(inputPar)
# #   }else{
# #     stop("speciesPars and speciesPars_ranges must contain vectors or a matrices")
# #   }
# #   if(Npar != nrow(parRange)){
# #     stop("speciesPars and speciesPars_ranges must have the same number of parameters")
# #   }
# #   checkedPar = matrix(nrow = Nsp, ncol = Npar)
# #   j=1
# #   for(j in 1:Nsp){
# #     for(i in 1:Npar){
# #       lower = parRange[i,1]
# #       upper = parRange[i,2]
# #       checkedPar[j,i] = inputPar[j,i] < upper & inputPar[j,i] > lower
# #     }
# #   }
# #   return(list(invalid = any(!checkedPar), checkedPar = checkedPar, inputPar = inputPar, parRange = parRange))
# # }
# #
# # checkParInput = function(speciesPars, speciesPars_ranges){
# #   checked = list()
# #   valid_pars = T
# #   for(i in names(speciesPars_ranges)){
# #     checked[[i]] = checkPars(speciesPars[[i]], speciesPars_ranges[[i]])
# #     if(checked[[i]]$invalid) valid_pars = F
# #   }
# #   if(!valid_pars){
# #     stop_message = paste("speciesPars must be within the range of speciesPars_ranges", sep = "\n")
# #     stop_message <- paste(stop_message, "The following parameters are out of range:", sep = "\n")
# #     for(i in names(speciesPars_ranges)){
# #       if(any(!checked[[i]]$checkedPar)){
# #         false_idx = which(!checked[[i]]$checkedPar, arr.ind = TRUE)
# #         for(ipar in 1:nrow(false_idx)){
# #           stop_message <- paste(stop_message, paste0("speciesPars$", i, "[", false_idx[ipar,1], ",", false_idx[ipar,2], "] = ", checked[[i]]$inputPar[false_idx[ipar,1], false_idx[ipar,2]]), sep = "\n")
# #         }
# #       }
# #     }
# #     stop_message <- paste(stop_message, "Adjust the upper and lower range in speciesPars_ranges or adjust the parameters in ", i, sep = "\n")
# #     stop(stop_message)
# #   }
# # }
#
# checkParInput(speciesPars, speciesPars_ranges)
#
inputPar = torch::torch_tensor(speciesPars$parMort)
inputPar = speciesPars$parMort
parRange = speciesPars_ranges$parMort
parRange = speciesPars_ranges$parReg
inputPar = speciesPars$parReg
# setPars <- function(inputPar, parRange){
#   # inputPar = torch::torch_tensor(inputPar, requires_grad = TRUE)
#   if(is.vector(inputPar)) {
#     Npar = 1
#     NPsp = length(inputPar)
#     inputPar = matrix(inputPar, nrow = NPsp, ncol = Npar)
#     parRange = matrix(parRange, nrow = Npar, ncol = 2)
#   }else if(is.matrix(inputPar)){
#     Npar = ncol(inputPar)
#     Nsp = nrow(inputPar)
#   }else{
#     stop("speciesPars and speciesPars_ranges must contain vectors or a matrices")
#   }
#   j=1
#   out = matrix(nrow = Nsp, ncol = 0)
#   # for(j in 1:Nsp){
#     for(i in 1:Npar){
#       lower = parRange[i,1,drop=FALSE]
#       upper = parRange[i,2,drop=FALSE]
#       out <- cbind(out, plogis(inputPar[,i, drop=FALSE])*as.numeric((upper - lower)+lower))
#     }
#   if(Npar == 1) out = as.vector(out) # TODO replace vectors with matrices as species input everywhere
#   # }
#   # out <- torch::torch_tensor(out, requires_grad = TRUE, device = "cpu", dtype = "float32")
#   out <- torch::torch_tensor(out, requires_grad = TRUE, device = self$device, dtype = self$dtype)
#   return(out)
# }
# setPars(inputPar, parRange)
#
# for(j in 1:Nsp){
#   for(i in 1:NPar){
#     inputPar[i, x, drop = FALSE]$sigmoid()*speciesPars_ranges$parMort[,x,drop=FALSE]
#   }
# }
#
#
# checked_parMort = lapply(1:nrow(self$parMort), function(x) self$parMort[, x, drop = FALSE]$sigmoid()*speciesPars_ranges$parMort[,x,drop=FALSE])
#



Npatches = 50
Nsites = 1
Ntimesteps = 500

site_dt <- data.table(
  expand.grid(
    list(
      siteID = 1:Nsites,
      year = 1:Ntimesteps
    )
  )
)
hist(rbeta(Ntimesteps*Nsites, 1, 10), xlim = c(0,1))
dist_dt <- site_dt
dist_dt$intensity = rbinom(Ntimesteps*Nsites, 1, 0)*rbeta(Ntimesteps*Nsites, 1, 10)

env_dt <- site_dt
env_dt$env1 = rep(1, Ntimesteps)
# env_dt$env1 = rep(seq(-2,2,length.out = Nsites), Ntimesteps)

system.time({
  predictions =
    simulateForest(env = env_dt,
                   disturbance = dist_dt,
                   sp = Nsp,
                   patches=Npatches,
                   height = speciesPars$parHeight,
                   # growthProcess = createProcess(~1+env1, func = growth, initEnv = list(matrix(c(runif(),3,4), Nsp, 2)), initSpecies = matrix(c(0.2, 3.9), Nsp, 2)),
                   # mortalityProcess = createProcess(~1+env1, func = mortality, initEnv = list(matrix(c(1, -3), Nsp, 2)), initSpecies = matrix(c(0.2, 3.9), Nsp, 2)),
                   # regenerationProcess = createProcess(~1+env1, func = regeneration, initEnv = list(matrix(c(4, 0), 1, 2)), initSpecies = matrix(c(0.2), Nsp)),
                   growthProcess = createProcess(~1+env1, func = growth, initEnv = speciesPars$parGrowthEnv, initSpecies = speciesPars$parGrowth),
                   mortalityProcess = createProcess(~1+env1, func = mortality, initEnv = speciesPars$parMortEnv, initSpecies = speciesPars$parMort),
                   regenerationProcess = createProcess(~1+env1, func = regeneration, initEnv = speciesPars$parRegEnv, initSpecies = speciesPars$parReg),
                   device = "cpu")
})

pred = pred2DF(predictions, format = "long")$site

library(ggplot2)
ggplot(pred[, .(value = mean(value) ), by = .(year, species, variable)], aes(x = year, y = value, color = factor(species))) +
# ggplot(pred[, value = ,], aes(x = year, y = value, group = siteID, color = siteID)) +
  geom_line() +
  labs(x = "Year",
       y = "value") +
  theme_minimal()+
  facet_wrap(~variable, scales = "free_y")

pred = pred2DF(predictions, format = "wide")$site
comp_dt <- merge(pred, env_dt, by = c("siteID", "year"))
comp_dt <- merge(comp_dt, dist_dt, by = c("siteID", "year"))
breaks = cut(comp_dt$year, breaks = 50)
plot(growth ~ env1, data = comp_dt, col = colorRampPalette(c("red","white", "blue"))(50)[as.integer(breaks)])
summary(lm(growth ~ env1, data = comp_dt))


plot(growth ~ siteID, data = comp_dt)

plot(intensity ~ env1, data = comp_dt)
plot(env1 ~ env1, data = comp_dt)
plot(intensity ~ siteID, data = comp_dt)

predictions$Predictions$Site

FINN::simulateForest()



finn = f(sp = sp, env = 2L, device = "cuda:0", which = "all" ,
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
                      species = initCohort$species[,,], env = env[,,], patches = patches, debug = FALSE, verbose = TRUE)
})







