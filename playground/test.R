#### TEST + Debugging
library(torch)
self = list()
self$sp = 5
self$env = 2
self$device = "cpu"
self$parHeight = NULL
self$parGrowth = NULL
self$parMort = NULL
self$parReg = NULL
self$parGrowthEnv = NULL
self$parMortEnv = NULL
self$parRegEnv = NULL
self$hidden_growth = list()
self$hidden_mort = list()
self$hidden_reg = list()
self$build_NN = build_NN

sp = self$sp
env = self$env
device = self$device
parHeight = self$parHeight
parGrowth = self$parGrowth
parMort = self$parMort
parReg = self$parReg
parGrowthEnv = self$parGrowthEnv
parMortEnv = self$parMortEnv
parRegEnv = self$parRegEnv
hidden_growth = self$hidden_growth
hidden_mort = self$hidden_mort
hidden_reg = self$hidden_reg
which = "all"



self = init_FINN(sp = 5, device = "cpu" )
self$nnMortEnv

patches = 50
env = torch::torch_randn(size = c(100, 50, 2))
pred_growth = NULL
pred_morth = NULL
pred_reg = NULL
debug = TRUE
record_time = 0

cohort = CohortMat$new(dims = c(100, 30, 10), sp = 5)
env = env
patches = 30L
dbh = cohort$dbh
nTree = cohort$nTree
Species = cohort$Species
response = "dbh"
predict(self,
        env = env,
        patches = 30L,
        dbh = cohort$dbh,
        nTree = cohort$nTree,
        Species = cohort$Species)
parHeight = self$parHeight
compF_P(dbh, Species, nTree, self$parHeight)

X = torch::torch_randn(size = c(100, 50, 2))
Y = torch::torch_ones(size = c(100, 50, 5, 2))
initCohort = cohort
initCohort = NULL
epochs = 2L
batch_size = 20L
learning_rate = 0.1
start_time = 0.5
patches = 30L
response = "dbh"

mortFP(dbh, Species, nTree+0.00001, self$parMort, pred_morth[,i,], AL) #.unsqueeze(3)

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


r = torch_tensor(1.0)
r$names = c()

wregeneration2 = function(species, parReg, pred, light) {
  regP = torch_sigmoid((light + (1-parReg) - 1)/1e-3) # TODO masking? better https://pytorch.org/docs/stable/generated/torch.masked_select.html
  environment = pred
  regeneration = sample_poisson_relaxed((regP*(environment[,NULL])$`repeat`(c(1, species$shape[2], 1))+0.2 )) # TODO, check if exp or not?! lambda should be always positive!
  regeneration = regeneration + regeneration$round()$detach() - regeneration$detach()
  return(regeneration)
}




res = replicate(50, {

  for(i in 1:100) {



    cat("Regeneration...\n")
    with_detect_anomaly({

      # species = initCohort2$species
      # parReg = finn2$parReg
      # pred = finn2$nnRegEnv(env)[,i,]
      # light = torch_zeros_like(initCohort2$dbh)
      #
      # regP <<- torch_sigmoid((light + (1-parReg) - 1)/1e-7) # TODO masking? better https://pytorch.org/docs/stable/generated/torch.masked_select.html
      # environment2 <<- pred
      # regeneration2 <<- sample_poisson_relaxed((regP*(environment2[,NULL])$`repeat`(c(1, species$shape[2], 1))+0.2 )) # TODO, check if exp or not?! lambda should be always positive!
      # #regeneration2 = 5*(regP*(environment2[,NULL])$`repeat`(c(1, species$shape[2], 1))+0.2 )
      # r = regeneration2 + regeneration2$round()$detach() - regeneration2$detach()


      r = regeneration(initCohort2$species, finn2$parReg, finn2$nnRegEnv(env)[,i,], torch_zeros_like(initCohort2$dbh))
      r<<-r
      cat("R: ", as.numeric(r$sum()))
      r$sum()$backward()
      grad_par = as.matrix(finn2$parReg$grad)
      grad_pred = as.matrix(finn2$nnRegEnv$parameters[[1]]$grad)
      print(grad_pred)
      print(grad_par)
    })

    # cat("Mortality...\n")
    #
    # with_detect_anomaly({
    #   mortality(initCohort2$dbh, initCohort2$species, initCohort2$trees, finn2$parMort, finn2$nnMortEnv(env)[,i,], torch_ones_like(initCohort2$dbh))$sum()$backward()
    #   grad_par = as.matrix(finn2$parMort$grad)
    #   grad_pred = as.matrix(finn2$nnMortEnv$parameters[[1]]$grad)
    # })
    # #
    # cat("Growth....\n")
    # with_detect_anomaly({
    #   growth(initCohort2$dbh, initCohort2$species, finn2$parGrowth,finn2$parMort, finn2$nnGrowthEnv(env)[,i,], torch_zeros_like(initCohort2$dbh))$sum()$backward()
    #   grad_par = as.matrix(finn2$parMort$grad)
    #   grad_par = as.matrix(finn2$parGrowth$grad)
    #   grad_pred = as.matrix(finn2$nnGrowthEnv$parameters[[1]]$grad)
    # })

  }
})


res = replicate(10000, {
  with_detect_anomaly({
    pars = torch_tensor(0.2, requires_grad = TRUE)
    lmbd = (pars+5-5)*torch_ones(c(200, 200))

    #samples=50
    temperature = 1e-2
    num_samples = 1

    samples <<- torch::torch_rand(c(num_samples,lmbd$shape)) #+ 0.001
    te <<- (samples$log()$negative()/lmbd)$cumsum(dim=1L)
    print(as.numeric(te$max()))
    relaxed_indicator <<- torch_sigmoid((1.0 - te)  )#/ temperature)
    r = relaxed_indicator$sum(1)
    #r = t$sum(1)

    #r = regeneration2 + regeneration2$round()$detach() - regeneration2$detach()
    r$sum()$backward()
    grad_par = as.matrix(pars)
  })

})



library(torch)
growth = function(dbh, species, parGrowth, parMort, pred, light){

  shade = torch_sigmoid((light + (1-parGrowth[,1][species]) - 1)/1e-1)
  print(shade)
  environment = index_species(pred, species)
  pred = (shade*environment)
  # growth = (1.- torch.pow(1.- pred,4.0)) * parGrowth[species,1]
  growth = pred/2 * parGrowth[,2][species] * ((parMort[,2][species]-dbh/100) / parMort[,2][species])^(2)
  # growth = parGrowth[species,1]
  # return torch.nn.functional.softplus(growth)
  return(torch_clamp(growth, min = 0.0))
}

growth(dbh = initCohort$dbh+1000, species = initCohort$species, parMort = cbind(0, 1000, 0)  |> torch_tensor() ,parGrowth = matrix(c(-100, 100), 1, 2) |> torch_tensor(), finn$nnGrowthEnv(env)[,1,], light = torch_ones(100, 1, 1) )$max()

(torch_ones_like(initCohort$dbh)+5 - mortality(dbh = initCohort$dbh+1000, species = initCohort$species, trees = torch_ones_like(initCohort$dbh)+5, parMort = cbind(0, 1000, 0)  |> torch_tensor(), finn$nnMortEnv(env)[,1,], light = torch_ones(100, 1, 1)))

regeneration(initCohort$species, parReg = -torch_ones(1), pred = finn$nnRegEnv(env)[,1,], light = torch_ones(100, 1, 1))





growth


pred =
  finn$predict(initCohort$dbh, initCohort$trees, initCohort$species,
               env = torch::torch_randn(size = c(100, 31, 2)),patches = 30)


finn$fit(initCohort = initCohort,epochs = 20,
         X =  torch::torch_randn(size = c(100, 31, 2)),
         Y = torch::torch_ones(size = c(100, 31, 15, 2)),batch_size = 100L,
         patches = 30L, learning_rate = 0.1)

initCohort = CohortMat$new(dims = c(100, 30, 2),
                           dbh = array(100, dim = c(100, 30, 2)),
                           trees = array(50, dim = c(100, 30, 2)),
                           sp = 3)
finn = FINN$new(sp = 3L, env = 2L, device = "cpu", which = "all", parMort = matrix(0, 3, 3), parGrowth = matrix(0.1, 3, 2))
finn$parMort
pred =
  finn$predict(initCohort$dbh, initCohort$trees, initCohort$species,
               env = torch::torch_randn(size = c(100, 31, 2)),patches = 30)
plot(torch::as_array(pred[[1]]$data())[1,,1])



