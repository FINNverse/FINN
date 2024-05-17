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

initCohort = CohortMat$new(dims = c(100, 1, 10),
                           dbh = array(1, dim = c(100, 1, 10)),
                           trees = array(1, dim = c(100, 1, 10)),
                           sp = 1)

finn = FINN$new(sp = 1L, env = 2L, device = "cpu", which = "all" ,parGrowth = matrix(555, 1, 2), parGrowthEnv = list(matrix(5.0, 1, 2)), parMortEnv = list(matrix(0, 1, 2)), parRegEnv = list(matrix(10, 1, 2)))
pred =
  finn$predict(initCohort$dbh, initCohort$trees, initCohort$species,
               env = torch::torch_randn(size = c(100, 200, 2)),patches = 1)
plot(torch::as_array(pred[[1]]$data())[1,,1])






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



