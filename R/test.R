#### TEST + Debugging
self = list()
self$sp = 5
self$env = 2
self$device = "cpu"
self$parGlobal = NULL
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
parGlobal = self$parGlobal
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



