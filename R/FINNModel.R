#' @title FINN: Forest Informed Neural Network
#'
#' @description
#' \loadmathjax
#'
#' The `FINNModel` class provides tools to initialize, train, and predict using a Forest Informed Neural Network. This model is designed for predicting tree growth, mortality, and regeneration across multiple species. The class supports various configurations, including the use of different devices (CPU or CUDA) and hidden layers in the neural network models.
#'
#' @include FINNBase.R
#' @export
FINNModel = R6::R6Class(
  classname = 'FINNModel',
  inherit = FINNBase,
  lock_objects = FALSE,
  lock_class = FALSE,
  public = list(

    #' @field sp (`integer(1)`) \cr
    #' Number of species
    sp = NULL,

    #' @field env (`integer(1)` or `list(3)` of `integer(1)`) \cr
    #' Number of environmental predictors
    env = NULL,

    #' @field output (`integer(1)` or `list(3)` of `integer(1)`) \cr
    #' Output dimensions of the NNs
    output = NULL,

    #' @field dtype (`character(1)`) \cr
    #' Device, "cpu" or cuda devices (e.g. "cuda:0")
    dtype = NULL,

    #' @field parameters (`list(n)` with n objects of class `torch_tensor`) \cr
    #' List of parameters that should be optimized
    parameters = NULL,

    #' @field optimizer (`torch_optimizer`) \cr
    #' First-order gradient descent optimizer, automatically set by the `fit()` method
    optimizer = NULL,

    #' @field history (`numeric(n)` for n epochs) \cr
    #' Vector of epoch losses
    history = NULL,

    #' @field param_history (`list(numeric(n))`) \cr
    #' Weights after each optimization step
    param_history = NULL,

    #' @field pred (`list(torch_tensor)`) \cr
    #' List of predictions
    pred = NULL,

    #' @field device (`character(1)` or `torch_device`) \cr
    #' Device, "cpu" or "cuda:0"
    device = 'cpu',

    #' @field device_r (`character(1)`) \cr
    #' Device, "cpu" or "cuda:0" (backup R object)
    device_r = "cpu",

    #' @field parComp (`torch_tensor`) \cr
    #' Competition species parameters, Dimensions: \mjseqn{[n_{species}]}
    parComp = NULL,

    #' @field parGrowth (`torch_tensor`) \cr
    #' Growth species parameters, Dimensions: \mjseqn{[n_{species}, m_{growth}]}
    parGrowth = NULL,

    #' @field parMort (`torch_tensor`) \cr
    #' Mortality species parameters, Dimensions: \mjseqn{[n_{species}, m_{growth}]}
    parMort = NULL,

    #' @field parReg (`torch_tensor`) \cr
    #' Regeneration species parameters, Dimensions: \mjseqn{[n_{species}]}
    parReg = NULL,

    #' @field parGrowthEnv (`list(torch_tensor, torch_tensor, ...)`) \cr
    #' Parameters for environmental growth model
    parGrowthEnv = NULL,

    #' @field parMortEnv (`list(torch_tensor, torch_tensor, ...)`) \cr
    #' Parameters for environmental mortality model
    parMortEnv = NULL,

    #' @field parRegEnv (`list(torch_tensor, torch_tensor, ...)`) \cr
    #' Parameters for environmental regeneration model
    parRegEnv = NULL,

    #' @field nnRegEnv (`nn_sequential`) \cr
    #' Neural Network for environmental regeneration model
    nnRegEnv = NULL,

    #' @field nnGrowthEnv (`nn_sequential`) \cr
    #' Neural Network for environmental growth model
    nnGrowthEnv = NULL,

    #' @field nnMortEnv (`nn_sequential`) \cr
    #' Neural Network for environmental mortality model
    nnMortEnv = NULL,

    #' @field hidden_growth (`list(integer(1), integer(1))`) \cr
    #' Each integer in the list is one hidden layer
    hidden_growth = list(),

    #' @field hidden_mort (`list(integer(1), integer(1))`) \cr
    #' Each integer in the list is one hidden layer
    hidden_mort = list(),

    #' @field hidden_reg (`list(integer(1), integer(1))`) \cr
    #' Each integer in the list is one hidden layer
    hidden_reg = list(),

    #' @field speciesPars_ranges (`list()`) \cr
    #' Ranges for species parameters
    speciesPars_ranges = NULL,

    #' @field runEnvGrowth (`logical(1)`) \cr
    #' Make predictions with GrowthNN on env or not
    runEnvGrowth = TRUE,

    #' @field runEnvMort (`logical(1)`) \cr
    #' Make predictions with MortNN on env or not
    runEnvMort = TRUE,

    #' @field runEnvReg (`logical(1)`) \cr
    #' Make predictions with RegNN on env or not
    runEnvReg = TRUE,

    #' @field bias (`logical(1)`) \cr
    #' Use bias in neural networks or not
    bias = FALSE,

    #' @field patch_size_ha (`numeric(1)`) \cr
    #' Patch size in ha
    patch_size_ha = 0.1,

    #' @field minLight (`numeric(1)`) \cr
    #' Minimum Light
    minLight = 50,

    #' @field growthFunction (`function()`) \cr
    #' Growth process function
    growthFunction = NULL,

    #' @field mortalityFunction (`function()`) \cr
    #' Mortality process function
    mortalityFunction = NULL,

    #' @field regenerationFunction (`function()`) \cr
    #' Regeneration process function
    regenerationFunction = NULL,

    #' @field competitionFunction (`function()`) \cr
    #' Competition function
    competitionFunction = NULL,

    #' @field nnMortConfig (`list()`) \cr
    #' Neural Network Configuration, internal useage
    nnMortConfig = NULL,

    #' @field nnGrowthConfig (`list()`) \cr
    #' Neural Network Configuration, internal useage
    nnGrowthConfig = NULL,

    #' @field nnRegConfig (`list()`) \cr
    #' Neural Network Configuration, internal useage
    nnRegConfig = NULL,

    #' @field year_sequence (`numeric(n)`) \cr
    #' Time points for loss calculation
    year_sequence = NULL,

    #' @description
    #' Initializes the FINNModel model with the specified parameters.
    #' @param sp integer. Number of species.
    #' @param env vector of input dimensions for the neural networks
    #' @param output vector of output dimensions for the neural networks
    #' @param dtype dtype, should be either `torch::torch_float32()` or `torch::torch_float64()`
    #' @param device character. Device to use ('cpu' or 'cuda').
    #' @param parComp torch.Tensor. Global parameters for competition
    #' @param parGrowth torch.Tensor. Growth parameters.
    #' @param parMort torch.Tensor. Mortality parameters.
    #' @param parReg torch.Tensor. Regeneration parameters.
    #' @param parGrowthEnv torch.Tensor. Growth environment parameters.
    #' @param parMortEnv torch.Tensor. Mortality environment parameters.
    #' @param parRegEnv torch.Tensor. Regeneration environment parameters.
    #' @param nnRegEnv list of matrices. Initial parameters for the NN.
    #' @param nnGrowthEnv list of matrices. Initial parameters for the NN.
    #' @param nnMortEnv list of matrices. Initial parameters for the NN.
    #' @param hidden_growth list. Hidden layers for the growth model.
    #' @param hidden_mort list. Hidden layers for the mortality model.
    #' @param hidden_reg list. Hidden layers for the regeneration model.
    #' @param runEnvGrowth logical. Should the runEnvGrowth NN be run outside of the process functions on the environment.
    #' @param runEnvMort logical. Should the runEnvMort NN be run outside of the process functions on the environment.
    #' @param runEnvReg logical. Should the runEnvReg NN be run outside of the process functions on the environment.
    #' @param speciesPars_ranges list. List of species parameter ranges.
    #' @param bias logical. Whether to include a bias term in the neural networks.
    #' @param patch_size_ha numeric. Patch size in hectares.
    #' @param minLight numeric. Minimum light level.
    #' @param growthFunction function or NULL. Which process function should be used. If NULL, package process functions will be used, see `?growth`
    #' @param mortalityFunction function or NULL. Which process function should be used. If NULL, package process functions will be used, see `?mortality`
    #' @param regenerationFunction function or NULL. Which process function should be used. If NULL, package process functions will be used, see `?regeneration`
    #' @param competitionFunction function or NULL. Which process function should be used. If NULL, package process functions will be used, see `?competition`

    #' @return Invisible self.
    initialize = function(
                  sp = NULL,
                  env = NULL,
                  output = NULL,
                  dtype = NULL,
                  device = 'cpu',
                  # parHeight = NULL, # must be dim [species]
                  parComp = NULL, # must be dim [species, 2], first for height alometry, second for species competition strength
                  parGrowth = NULL, # must be dim [species, 2], first for shade tolerance
                  parMort = NULL, # must be dim [species, 2], first for shade tolerance,
                  parReg = NULL, # must be dim [species]
                  parGrowthEnv = NULL, # must be dim [species, 2], first for shade tolerance
                  parMortEnv = NULL, # must be dim [species, 2], first for shade tolerance,
                  parRegEnv = NULL, # must be dim [species]
                  nnRegEnv = NULL,
                  nnGrowthEnv = NULL,
                  nnMortEnv = NULL,
                  hidden_growth = list(),
                  hidden_mort = list(),
                  hidden_reg = list(),
                  runEnvGrowth = TRUE,
                  runEnvMort = TRUE,
                  runEnvReg = TRUE,
                  speciesPars_ranges = NULL,
                  bias = FALSE,
                  patch_size_ha = 0.1,
                  minLight = 50,
                  growthFunction = NULL,
                  mortalityFunction = NULL,
                  regenerationFunction = NULL,
                  competitionFunction = NULL
      ){

      # check input parameters
      # if(any(parHeight < 0 | parHeight > 1)) stop("parHeight cannot be <0 or >1")
      # if(any(parReg < 0 | parReg > 1)) stop("parReg cannot be <0 or >1")

      self$sp = sp
      self$device = device


      # ###### Defaults #####
      if(is.null(speciesPars_ranges)) speciesPars_ranges = default_speciesPars_ranges
      # if(is.null(parHeight)) parHeight = runif(sp, min = 0.45, max = 0.55)
      if(is.null(parComp)) parComp = cbind(
        runif(sp, min = 0.45, 0.55),
        runif(sp, min = 0.2, 0.2) # 0.2 corresponds to 0 light at 50m2/ha basal area
      )
      if(is.null(parGrowth)) parGrowth = cbind(
        runif(sp, min = 0.5, 0.55),
        runif(sp, min = 1.90, 2.0)
      )
      if(is.null(parMort)) parMort = cbind(
        runif(sp, min = 0.5, 0.55),
        runif(sp, min = 1.90, 2.0)
      )
      if(is.null(parReg)) parReg = runif(sp, min = 0.45, max = 0.55)
      self$speciesPars_ranges = speciesPars_ranges

      # if(!is.null(speciesPars_ranges)){
      checkParInput(
        speciesPars = list(
          parGrowth = parGrowth,
          parMort = parMort,
          parReg = parReg,
          parComp = parComp
        ),
        speciesPars_ranges = list(
          parGrowth = speciesPars_ranges$parGrowth,
          parMort = speciesPars_ranges$parMort,
          parReg = speciesPars_ranges$parReg,
          parComp = speciesPars_ranges$parComp
        )
      )
      # }

      # self$set_parHeight(parHeight)
      # self$set_parGrowth(parGrowth)
      # self$set_parMort(parMort)
      # self$set_parReg(parReg)

      self$parGrowth = self$setPars(parGrowth, speciesPars_ranges$parGrowth)
      self$parMort = self$setPars(parMort, speciesPars_ranges$parMort)
      self$parReg = self$setPars(parReg, speciesPars_ranges$parReg)
      self$parComp = self$setPars(parComp, speciesPars_ranges$parComp)

      self$parGrowthEnv = parGrowthEnv
      self$parMortEnv = parMortEnv
      self$parRegEnv = parRegEnv
      self$hidden_growth = hidden_growth
      self$hidden_mort = hidden_mort
      self$hidden_reg = hidden_reg
      self$patch_size_ha = patch_size_ha
      self$minLight = minLight
      self$device_r = device
      self$device = torch_device(device)
      self$sp = sp
      self$pred = NULL
      if(length(env) == 1) env = rep(env, 3)
      self$env = env
      self$optimizer = NULL
      self$dtype = torch_float32()
      self$runEnvGrowth = runEnvGrowth
      self$runEnvMort = runEnvMort
      self$runEnvReg = runEnvReg
      if(is.null(output)) output = c(sp, sp, sp)
      self$output = output

      self$nnMortConfig = list(input_shape=env[1], output_shape=output[1], hidden=hidden_mort, activation="selu", bias=bias, dropout=-99, last_activation = "linear")
      self$nnMortEnv = do.call(self$build_NN, self$nnMortConfig)

      self$nnGrowthConfig = list(input_shape=env[2], output_shape=output[2], hidden=hidden_growth, activation="selu", bias=bias, dropout=-99, last_activation = "linear")
      self$nnGrowthEnv = do.call(self$build_NN, self$nnGrowthConfig)

      self$nnRegConfig = list(input_shape=env[3], output_shape=output[3], hidden=hidden_reg, activation="selu", bias=bias, dropout=-99, last_activation = "linear")
      self$nnRegEnv = do.call(self$build_NN, self$nnRegConfig)

      if(!is.null(parGrowthEnv)) self$set_weights_nnGrowthEnv(parGrowthEnv)
      if(!is.null(parMortEnv)) self$set_weights_nnMortEnv(parMortEnv)
      if(!is.null(parRegEnv)) self$set_weights_nnRegEnv(parRegEnv)

      if(is.null(growthFunction)) self$growthFunction = growth |> self$set_environment()
      else self$growthFunction = growthFunction |> self$set_environment()

      if(is.null(mortalityFunction)) self$mortalityFunction = mortality |> self$set_environment()
      else self$mortalityFunction = mortalityFunction |> self$set_environment()

      if(is.null(regenerationFunction)) self$regenerationFunction = regeneration |> self$set_environment()
      else self$regenerationFunction = regenerationFunction |> self$set_environment()

      if(is.null(competitionFunction)) self$competitionFunction = competition |> self$set_environment()
      else self$competitionFunction = competitionFunction |> self$set_environment()

      self$nnMortEnv$to(device = self$device)
      self$nnGrowthEnv$to(device = self$device)
      self$nnRegEnv$to(device = self$device)

      self$parGrowthEnv = self$nnGrowthEnv$parameters
      self$parMortEnv = self$nnMortEnv$parameters
      self$parRegEnv = self$nnRegEnv$parameters

      # if(is.null(parHeight)){
      #   parHeight = np_runif(0.3, 0.7, size = self$sp)
      # }
      if(is.null(parGrowth)){
        parHeight = np_runif(0.3, 0.7, size = c(self$sp,1))
        parCompStr = np_runif(0, 2, size = c(self$sp,1))
        parComp = cbind(parHeight, parCompStr)
      }
      if(is.null(parGrowth)){
        first = np_runif(0, 6, size = c(self$sp,1))
        second = np_runif(0, 6, size = c(self$sp,1))
        parGrowth = cbind(first, second)
      }

      if(is.null(parMort)){
        first = np_runif(0, 2, size = c(self$sp,1))
        second = np_runif(0.1, 5, size = c(self$sp,1))
        parMort = cbind(first, second)
      }

      if(is.null(parReg)){
        self$parReg = np_runif(0, 1, size = self$sp)
      }

      # self$set_parHeight(parHeight)
      # self$set_parGrowth(parGrowth)
      # self$set_parMort(parMort)
      # self$set_parReg(parReg)

      # browser()

      self$parGrowth = self$setPars(parGrowth, speciesPars_ranges$parGrowth)
      self$parMort = self$setPars(parMort, speciesPars_ranges$parMort)
      self$parReg = self$setPars(parReg, speciesPars_ranges$parReg)
      # self$parHeight = self$setPars(parHeight, speciesPars_ranges$parHeight)
      self$parComp = self$setPars(parComp, speciesPars_ranges$parComp)

      self$scale_2 = torch::torch_tensor(1.0, requires_grad = TRUE, dtype = self$dtype, device = self$device)
      self$scale_4 = torch::torch_tensor(1.0, requires_grad = TRUE, dtype = self$dtype, device = self$device)
      self$scale_5 = torch::torch_tensor(1.0, requires_grad = TRUE, dtype = self$dtype, device = self$device)
      self$scale_6 = torch::torch_tensor(1.0, requires_grad = TRUE, dtype = self$dtype, device = self$device)

      # self$set_parHeight(parHeight)
      # self$set_parGrowth(parGrowth)
      # self$set_parMort(parMort)
      # self$set_parReg(parReg)

      self$parameters = c(
        self$parComp, self$parGrowth, self$parMort,self$parReg, self$nnRegEnv$parameters,
        self$nnGrowthEnv$parameters, self$nnMortEnv$parameters,self$scale_2, self$scale_4,
        self$scale_5, self$scale_6)
      names(self$parameters) = c(
        "parComp" , "parGrowth", "parMort", "parReg", "nnReg",
        "nnGrowth", "nnMort", "scale_2",
        "scale_4", "scale_5", "scale_6")
      self$parameter_to_r()

      return(invisible(self)) # Only for testing now
    },


    #' @description
    #' Predicts the growth, mortality, and regeneration of trees based on the given inputs.
    #'
    #' The `predict` method generates predictions for tree growth, mortality, and regeneration for the specified species across different environmental conditions. It uses the initialized model parameters and can handle optional input tensors like diameter at breast height (dbh), number of trees, and species. If these are not provided, they will be initialized internally.
    #'
    #' @param dbh torch.Tensor (Optional). Diameter at breast height of the trees.
    #' @param trees torch.Tensor (Optional). Number of trees.
    #' @param species torch.Tensor (Optional). Species of the trees.
    #' @param env torch.Tensor. Environmental data.
    #' @param disturbance torch.Tensor. Disturbance rates.
    #' @param start_time integer. Time at which to start recording the results.
    #' @param pred_growth torch.Tensor (Optional). Predicted growth values.
    #' @param pred_mort torch.Tensor (Optional). Predicted mortality values.
    #' @param pred_reg torch.Tensor (Optional). Predicted regeneration values.
    #' @param patches numeric. Number of patches.
    #' @param debug logical. Run in debug mode if TRUE.
    #' @param y torch.Tensor. Response tensor for target data.
    #' @param c torch.Tensor. Number of tree tensor.
    #' @param update_step integer. Backpropagation step length.
    #' @param verbose logical. Print progress if TRUE.
    #' @param weights weights for reweighting the loss
    #' @param year_sequence at which year indices should the predictions compared with the observed values
    #' @return list. A list of predicted values for dbh, number of trees, and other recorded time points. If `debug` is TRUE, raw results and cohorts are also returned.
    predict = function(dbh = NULL,
                       trees = NULL,
                       species = NULL,
                       env = NULL,
                       disturbance = NULL,
                       start_time = 1L,
                       pred_growth = NULL,
                       pred_mort = NULL,
                       pred_reg = NULL,
                       patches = 100L,
                       debug = FALSE,
                       y = NULL,
                       c = NULL,
                       update_step = 1L,
                       verbose = TRUE,
                       weights = NULL,
                       year_sequence = NULL){


      # if no cohorts exist initialize empty cohort array
      if(is.null(dbh)){
        cohorts = CohortMat$new(dims = c(env$shape[1], patches, self$sp), sp = self$sp, device=self$device)
        trees = cohorts$trees
        species = cohorts$species
        dbh = cohorts$dbh
      }

      # repeat env if same env for the three processes
      if(is.list(env)) {
        env = lapply(env, function(e) torch_tensor(e, dtype=self$dtype, device=self$device))
      } else {
        # processes 1, 2, 3 are  mortality, growth, and regeneration
        env = lapply(1:3, function(i) torch_tensor(env, dtype=self$dtype, device=self$device))
      }
      names(env) <- c("mort", "growth", "reg")

      # get dimensions
      sites = env[[1]]$shape[1]
      time =  env[[1]]$shape[2]
      patches = dbh$shape[2]
      sp = self$sp

      # check dtype and device of disturbance
      if(!is.null(disturbance)) {
        disturbance = disturbance$to(dtype=self$dtype, device=self$device)
        disturbances_tens = torch::distr_bernoulli(probs = disturbance$squeeze(3L))$sample(patches)$permute(c(2, 3, 1))
        disturbances_tens = 1*(disturbances_tens==0)
      }

      # if y (response) is null -> simulation mode, gradients are not required (not neccessary to set them to zero but it is cleaner)
      if(is.null(y)) {
        lapply(self$parameters, function(p) p$requires_grad_(FALSE) )
      }

      # loss weights / weighting of the different losses, must be cast to a tensor
      if(!is.null(y)) {
        # if weights are zero, set all weights to 1
        if(is.null(weights)) weights = torch_ones(dim(y)[4], dtype = self$dtype, device = self$device)
      }

      # send data to correct devices (and dtypes)
      dbh = torch_tensor(dbh, dtype=self$dtype, device=self$device)
      trees = torch_tensor(trees, dtype=self$dtype, device=self$device)
      species = torch_tensor(species, dtype=torch_int64(), device=self$device)

      # create cohort ids
      cohort_ids = torch_tensor(array(
        1:(prod(species$shape)+1),
        dim = species$shape), dtype=torch_int32(), device = self$device
      )

      # init Result tensors
      Result = lapply(1:8,function(tmp) torch_zeros(list(sites, time, self$sp), device=self$device))
      names(Result) =  c("dbh","ba", "trees", "AL", "growth", "mort", "reg", "r_mean_ha")

      # if debug modus, create empty lists
      if(debug) {
        Raw_cohort_results = list()
        Raw_cohort_ids = list()
        Raw_patch_results = list()
      }

      # if simulation mode, environmental predictions can be made a priori (it would interrupt the gradients in inference mode)
      if(is.null(y)) {
        if(self$runEnvGrowth) predGrowthGlobal = self$nnGrowthEnv(env[["growth"]])
        if(self$runEnvMort) predMortGlobal = self$nnMortEnv(env[["mort"]])
        if(self$runEnvReg) predRegGlobal = self$nnRegEnv(env[["reg"]])
      }

      # create process bar
      if(verbose) cli::cli_progress_bar(format = "Year: {cli::pb_current}/{cli::pb_total} {cli::pb_bar} ETA: {cli::pb_eta} ", total = time, clear = FALSE)

      for(i in 1:time){

        # In inference mode, make env predictions in each time step (to get the gradients)
        # otherwise, just take the i-th prediction
        if(!is.null(y)) {
          if(self$runEnvGrowth) pred_growth = self$nnGrowthEnv(env[["growth"]][,i,])
          if(self$runEnvMort) pred_mort = self$nnMortEnv(env[["mort"]][,i,])
          if(self$runEnvReg) pred_reg = self$nnRegEnv(env[["reg"]][,i,])
        } else {
          if(self$runEnvGrowth) pred_growth = predGrowthGlobal[,i,]
          if(self$runEnvMort) pred_mort = predMortGlobal[,i,]
          if(self$runEnvReg) pred_reg = predRegGlobal[,i,]
        }


        # empty rate objects/tensors
        light = torch_zeros(list(sites, time,  dbh$shape[3]), device=self$device)
        g = torch_zeros(list(sites, time, dbh$shape[3]), device=self$device)
        m = torch_zeros(list(sites, time, dbh$shape[3]), device=self$device)
        r = torch_zeros(list(sites, time, dbh$shape[3]), device=self$device)
        if(debug) trees_before = torch::torch_zeros_like(g)

        # detach previous cohort objects (to interrupt the gradients)
        dbh=dbh$detach()
        trees=trees$detach()
        species=species$detach()
        cohort_ids=cohort_ids$detach()

        # Apply disturbance
        if(!is.null(disturbance)) {
          trees = trees*disturbances_tens[,i,]$unsqueeze(3L)
        }

        # Model - get Parameters parameter constrains
        # TODO change to active bindings see https://collinerickson.github.io/2018/01/10/using-active-bindings-in-r6-classes/
        parMort = self$getPars(self$parMort, self$speciesPars_ranges$parMort)
        parGrowth = self$getPars(self$parGrowth, self$speciesPars_ranges$parGrowth)
        parReg = self$getPars(self$parReg, self$speciesPars_ranges$parReg)
        # parHeight = self$getPars(self$parHeight, self$speciesPars_ranges$parHeight)
        parComp = self$getPars(self$parComp, self$speciesPars_ranges$parComp)


        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## Demographic processes ####
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

        # check if there is at least one cohort alive (otherwise just jump directly to regeneration?)
        if(dbh$shape[3] > 0.5){
          # calculate available light for each cohort
          light = self$competitionFunction(
            dbh = dbh,
            species = species,
            trees = trees,
            # parHeight = parHeight,
            parComp = parComp,
            h = NULL,
            # minLight = self$minLight,
            patch_size_ha = self$patch_size_ha,
            ba = NULL,
            cohortHeights = NULL
          )
          # Convert tensor comparison result to a logical scalar
          if (light$gt(1.0)$sum()$item() > 0) {
            stop("Light values > 1")
          }

          if(self$runEnvGrowth) pred = index_species(pred_growth, species)
          else pred = env[["growth"]][,i,]

          g = self$growthFunction(
            dbh = dbh,
            species = species,
            parGrowth = parGrowth,
            pred = pred,
            light = light
          )
          # print("Mort:")
          # print(g)
          # print("++++")
          dbh_growth = dbh*g
          dbh = dbh + dbh_growth

          light = self$competitionFunction(
            dbh = dbh,
            species = species,
            trees = trees,
            parComp = parComp,
            h = NULL,
            # minLight = self$minLight,
            patch_size_ha = self$patch_size_ha,
            ba = NULL,
            cohortHeights = NULL
          )
          # cat("Second section  B2\n")

          if(self$runEnvMort) pred = index_species(pred_mort, species)
          else pred = env[[2]][,i,]

          m = self$mortalityFunction(
            dbh = dbh,
            species = species,
            trees = trees + 0.001,
            parMort = parMort,
            pred = pred,
            light = light
          )
          # print("Mort:")
          # print(m)
          # print("++++")
          trees_dead = binomial_from_gamma(trees+trees$le(0.5)$float(), m+0.001)*trees$ge(0.5)$float() #rbinom(n, trees, m) --> m*trees
          trees_dead = trees_dead + trees_dead$round()$detach() - trees_dead$detach()
          trees_before = trees

          #.unsqueeze(3) # TODO check!
          #trees$sub_(m)$clamp_(min = 0.0)
          trees = torch_clamp(trees - trees_dead, min = 0) #### TODO if trees = 0 then NA...prevent!
        }

        # start reg
        AL_reg = self$competitionFunction( # must have dimension = n species in last dim
          dbh = dbh,
          species = species,
          trees = trees,
          parComp = parComp,
          h = 1,
          # minLight = self$minLight,
          patch_size_ha = self$patch_size_ha,
          ba = NULL,
          cohortHeights = NULL
        )

        if(self$runEnvReg) pred = pred_reg
        else pred = env[["reg"]][,i,]

        r_mean = self$regenerationFunction(species = species,
                                           parReg = parReg,
                                           pred = pred,
                                           light = AL_reg)


        r = sample_poisson_gaussian(r_mean*self$patch_size_ha)
        # ein nummerischer Trick um den Gradienten für die Zahlen beim Runden zu behalten
        r = r + r$round()$detach() - r$detach()

        # print("Mort:")
        # print(r)
        # print("++++")

        # New recruits
        dbh_new = ((r-1+0.1)/1e-3)$sigmoid() # TODO: check!!! --> when r 0 dann dbh = 0, ansonsten dbh = 1 dbh[r==0] = 0
        trees_new = r
        species_new = torch_arange(1, sp, dtype=torch_int64(), device = self$device)$unsqueeze(1)$`repeat`(c(r$shape[1], r$shape[2], 1))

        # assign cohortIDs
        max_id = max(c(1,as_array(cohort_ids), na.rm = TRUE))
        new_cohort_id = torch_tensor(array(
          (max_id+1):(max_id+prod(r$shape)+1),
          dim = r$shape), dtype=torch_int32(), device = self$device
        ) #TODO check for performance

        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## Aggregation of rates####
        # from previous cohort
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        if(dbh$shape[3] != 0){
          samples = vector("list", 4)
          names(samples) = c(
            "light*trees_before",
            "g*trees_before",
            "m*trees_before",
            "trees_before"
          )
          samples[[1]] = light*trees_before
          samples[[2]] = g*trees_before
          samples[[3]] = m*trees_before
          samples[[4]] = trees_before # original number of trees

          Results_tmp = replicate(length(samples), torch_zeros_like(Result[[1]][,i,]))
          # calculate the sum of all elements in samples over all patches for each species, patch and site
          tmp_res = aggregate_results(species, samples, Results_tmp, aggregation = "sum")

          # Aggregation [sites, patches, cohorts] --> [sites, species]:
          if(as.numeric(trees$sum() > 0)) {
            # Light
            mask = tmp_res[[4]]$gt(0.5)

            Result[[4]][,i,][mask] =  Result[[4]][,i,][mask] + tmp_res[[1]][mask]/tmp_res[[4]][mask]

            # Growth
            Result[[5]][,i,][mask] = Result[[5]][,i,][mask]+tmp_res[[2]][mask]/tmp_res[[4]][mask] # summe dbh_growth / summe trees  #/cohort_counts[alive_species])

            # Mort
            Result[[6]][,i,][mask] = Result[[6]][,i,][mask]+tmp_res[[3]][mask]/tmp_res[[4]][mask] # summe dbh_growth / summe trees  #/cohort_counts[alive_species])
            # Result[[6]][,i,] = Result[[6]][,i,]$add(tmp_res[[3]]/patches) # summe dbh_growth / summe trees  #/cohort_counts[alive_species])
            #Result[[6]][,i,][Result[[6]][,i,]$isinf()] = 0.0
          }
        }
        # reg extra
        ## Regeneration count
        tmp_res = aggregate_results(species_new, list(r), list(torch::torch_zeros(Result[[1]][,i,]$shape[1], sp, device = self$device )))
        Result[[7]][,i,] = Result[[7]][,i,]$add(tmp_res[[1]])/patches

        ## Regeneration rate mean
        tmp_res = aggregate_results(species_new, list(r_mean), list(torch::torch_zeros(Result[[1]][,i,]$shape[1], sp, device = self$device )))
        r_mean = tmp_res[[1]]/patches
        Result[[8]][,i,] = r_mean
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## Update arrays ####
        # 1. Combine old cohorts and recruit cohorts
        # 2. Find dead cohorts and remove them, if possible (by finding the minimal required cohort dimension),
        #    a few dead cohorts are preserved for padding.
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

        dbh = torch::torch_cat(list(dbh, dbh_new), 3)
        trees = torch::torch_cat(list(trees, trees_new), 3)
        species = torch::torch_cat(list(species, species_new), 3)
        cohort_ids = torch::torch_cat(list(cohort_ids, new_cohort_id), 3)

        if (debug) {
          tmp = torch::torch_zeros_like(dbh_new)
          tmp[] = NaN
          g = torch::torch_cat(list(g, tmp), 3)
          m = torch::torch_cat(list(m, tmp), 3)
          trees_before = torch::torch_cat(list(trees_before, tmp), 3)
        }
        rm(dbh_new,trees_new,species_new,new_cohort_id )

        # Pad tensors, expensive, currently each timestep
        if(i %% 1 == 0){

          # Gradient shouldn't be required, also expensive for backpropagation because of reshape/view operations!
          #torch::with_no_grad({
          # Masks to find alive cohorts
          mask = (trees > 0.5)$flatten(start_dim = 1, end_dim = 2)
          org_dim = species$shape[1:2]
          org_dim_t = torch::torch_tensor(org_dim, dtype = torch_long(), device = "cpu")

          # Minimal number of cohort dimension
          non_zero_counts = mask$sum(dim=2)
          max_non_zeros = non_zero_counts$max()$item()
          sorted_tensor = torch::torch_sort(mask$float(), dim=2, descending=TRUE)[[2]]

          dbh = dbh$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
          trees = trees$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
          cohort_ids = cohort_ids$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
          species = species$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)

          if(debug) {
            g = g$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
            m = m$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
            trees_before = trees_before$flatten(start_dim = 1, end_dim = 2)$gather(2, sorted_tensor)[, 1:max_non_zeros]$unflatten(1, org_dim)
          }

          # })
        }

        # aggregate results
        # TODO: Position of the first block?
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## Aggregation of stand variables ####
        # from updated cohorts
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        if(i > 0){
          if(dbh$shape[3] != 0){

              #dead_trees_mask = trees == 0
              dbh = dbh*trees$gt(0.5)$float()

              BA_stem_values = BA_stem(dbh = dbh)*trees
              species = species
              samples = vector("list", 3)
              samples[[1]] = trees * dbh #* mask
              samples[[2]] = BA_stem_values #* mask# torch_sigmoid((trees - 0.5)/1e-3) # better to just use greater? (Masking!) Gradients shouldn't be needed! (I think?)
              samples[[3]] = trees # torch_sigmoid((trees - 0.5)/1e-3)
              Results_tmp = replicate(length(samples), torch_zeros_like(Result[[1]][,i,]))

              tmp_res = aggregate_results(species, samples, Results_tmp)
              # BA and number of trees Result[[1]] and Result[[2]]
              if(as.numeric(trees$sum() > 0)) {
                # mask only alive trees
                mask = tmp_res[[3]]$gt(0.5)

                Result[[1]][,i,][mask] = Result[[1]][,i,][mask]+ (tmp_res[[1]][mask]/tmp_res[[3]][mask]) # /cohort_counts[alive_species]

                # BA
                Result[[2]][,i,]= Result[[2]][,i,]$add(tmp_res[[2]]$div_(patches))

                # Trees
                Result[[3]][,i,]= Result[[3]][,i,]$add(tmp_res[[3]]$div_(patches))
              }
              rm(BA_stem_values)
          }
        }
        # end aggregation


        # if (debug) {
        #   Raw_cohort_results[[i]]$species = torch::as_array(species$cpu())
        #   Raw_cohort_results[[i]]$trees = torch::as_array(trees$cpu())
        #   Raw_cohort_results[[i]]$dbh = torch::as_array(dbh$cpu())
        #   Raw_cohort_ids[[i]] = torch::as_array(cohort_ids$cpu())
        # }

          if (debug) {
            Raw_cohort_results[[i]] = list(
              "species" = torch::as_array(species$cpu()),
              "trees" = torch::as_array(trees$cpu()),
              "dbh" = torch::as_array(dbh$cpu()),
              "m" = torch::as_array(m$cpu()),
              "g" = torch::as_array(g$cpu()),
              "trees_before" = torch::as_array(trees_before$cpu()) # noob
            )
            Raw_patch_results[[i]] = list(
              "r" = torch::as_array(r$cpu())
              # "r_mean" = torch::as_array(r_mean$cpu())
            )
            Raw_cohort_ids[[i]] = torch::as_array(cohort_ids$cpu())
          }



        loss = torch_zeros(6L, device = self$device)
        if(i > 0 && dbh$shape[3] != 0 && !is.null(y) && (i %% update_step == 0)) {
          if(i %in% year_sequence) {
            tmp_index = which(year_sequence %in% i, arr.ind = TRUE)
            for(j in 2:7) {
              # 1 -> dbh, 2 -> ba, 3 -> counts, 4 -> AL, 5 -> growth rates, 6 -> mort rates, 7 -> reg rates
              if(j == 3) {
                mask = c[, tmp_index,]$isnan()$bitwise_not()
                if(as.logical(mask$max()$data()))  loss[j-1] = loss[j-1]+
                    torch::distr_poisson(Result[[2]][,i,][mask]+0.01)$log_prob(c[, tmp_index,][mask])$mean()$negative()*(weights[j-1]+0.0001)
              } else if(j == 7) {
                mask = y[, tmp_index,,j]$isnan()$bitwise_not()
                if(as.logical(mask$max()$data()))  loss[j-1] = loss[j-1]+
                    torch::distr_poisson((r_mean+0.001)[mask])$log_prob(y[, tmp_index,,j][mask])$mean()$negative()*(weights[j-1]+0.0001)
              } else {
                mask = y[, tmp_index,,j]$isnan()$bitwise_not()
                #if(as.logical(mask$max()$data())) loss[j-1] = loss[j-1]+torch::nnf_mse_loss(y[, tmp_index,,j][mask], Result[[j]][,i,][mask])$mean()*(weights[j-1]+0.0001)
                if(as.logical(mask$max()$data()))
                loss[j-1] = loss[j-1]+ torch::distr_normal(Result[[j]][,i,][mask],  self$parameters[[paste0("scale_",j)]]$relu()+0.0001)$log_prob(y[, tmp_index,,j][mask])$mean()$negative()*(weights[j-1]+0.0001) + 0.01*(self$parameters[[paste0("scale_",j)]]$relu()+0.0001)**2
                #loss[j-1] = loss[j-1]+ torch::nnf_mse_loss(Result[[j]][,i,][mask],  y[, tmp_index,,j][mask])

              }
            }
            loss$sum()$backward()
          }

          # cat("Backprop....\n")
          for(j in 1:8) Result[[j]] = Result[[j]]$detach()
        }
        loss$detach_()

        if(verbose) cli::cli_progress_update()

      }

      #browser()

      names(Result) =  c("dbh","ba", "trees", "AL", "growth", "mort", "reg", "r_mean_ha")
      if(debug){
        Result_out = list(
          Predictions = list(
            Site = lapply(Result, function(x) as_array(x)),
            Patch = Raw_patch_results,
            Cohort = list(
              cohortStates = Raw_cohort_results,
              cohortID = Raw_cohort_ids
              # cohortID = lapply(Raw_cohort_ids, function(x) torch::as_array(x))
              )
            ),
          loss = loss)
      }else if(!debug){
        Result_out = list(
          Predictions = list(
            Site = lapply(Result, function(x) torch::as_array(x))),
          loss = loss)
      }

      if(is.null(y)) {
        lapply(self$parameters, function(p) p$requires_grad_(TRUE) )
      }

      return(Result_out)
    },
    #' @description
    #' Fits the FINN model to the provided training data.
    #'
    #' The `fit` method trains the FINN model using the input (X) and target (Y) data. It iteratively adjusts the model parameters to minimize the loss over the specified number of epochs. This method also supports customizing batch size, learning rate, and the starting time for predictions.
    #'
    #' @param X torch.Tensor (Optional). Input data, typically environmental variables over time.
    #' @param Y torch.Tensor (Optional). Target data, representing observed tree metrics.
    #' @param disturbance torch.Tensor. Disturbance rates.
    #' @param initCohort list (Optional). Initial cohort data, including initial dbh, trees, and species information.
    #' @param epochs integer. Number of epochs to train the model. Default is 2.
    #' @param batch_size integer. Batch size for training. Default is 20.
    #' @param learning_rate float. Learning rate for the optimizer. Default is 0.1.
    #' @param start_time float. The starting time for predictions. Default is 1.
    #' @param patches integer. Number of patches in the dataset. Default is 50.
    #' @param response character. Response variable to predict ("dbh", "BA_stem", or "BA_stem*nT"). Default is "dbh".
    #' @param update_step integer. Backpropagation step length.
    #' @param weights reweight likelihood
    #' @param year_sequence at which year indices should the predictions compared with the observed values
    #' @param file path if weights should be saved after each epoch (for monitoring).
    #' @param thin thin saving of weights
    #' @param optimizer torch optimizer
    #' @return None. The trained model is stored within the class instance.
    fit = function(
          X = NULL,
          Y = NULL,
          disturbance = NULL,
          initCohort = NULL,
          epochs = 2L,
          batch_size = 20L,
          learning_rate = 0.1,
          start_time = 1,
          patches = 100,
          update_step = 1L,
          weights = NULL,
          year_sequence = NULL,
          file = NULL,
          thin = 100L,
          optimizer = NULL){

        if(is.null(optimizer)) optimizer = torch::optim_adagrad

        if(is.null(self$optimizer)){
          self$optimizer = optimizer(params = self$parameters, lr = learning_rate) # AdaBound was also good
          # TODO scheduler implementieren
          # self.scheduler = torch.optim.lr_scheduler.ExponentialLR(self.optimizer, gamma=0.9)
        }

        if(!is.list(X)) X = list(X, X, X)
        stepSize = floor(X[[1]]$shape[1] / batch_size)
        time = X[[1]]$shape[2]

        if(self$device$type == 'cuda'){
          #torch.cuda.set_device(self$device)
          pin_memory = FALSE
        }else{
          pin_memory = TRUE
        }
        Counts = (Y[,,,3])$round() # Trees!!!
        if(is.null(disturbance)) {
          data = tensor_dataset(torch_tensor(X[[1]], dtype=self$dtype, device=torch_device('cpu')),
                                torch_tensor(X[[2]], dtype=self$dtype, device=torch_device('cpu')),
                                torch_tensor(X[[3]], dtype=self$dtype, device=torch_device('cpu')),
                                torch_tensor(Y[,,,], dtype=self$dtype, device=torch_device('cpu')),
                                torch_tensor(Counts, dtype=self$dtype, device=torch_device('cpu')),
                                torch_arange(1, X[[1]]$shape[1])$to(dtype = torch_int64() , device=torch_device('cpu'))
          )
        } else {
          data = tensor_dataset(torch_tensor(X[[1]], dtype=self$dtype, device=torch_device('cpu')),
                                torch_tensor(X[[2]], dtype=self$dtype, device=torch_device('cpu')),
                                torch_tensor(X[[3]], dtype=self$dtype, device=torch_device('cpu')),
                                torch_tensor(Y[,,,], dtype=self$dtype, device=torch_device('cpu')),
                                torch_tensor(Counts, dtype=self$dtype, device=torch_device('cpu')),
                                torch_arange(1, X[[1]]$shape[1])$to(dtype = torch_int64() , device=torch_device('cpu')),
                                torch_tensor(disturbance, dtype=self$dtype, device=torch_device('cpu'))
          )
        }
        DataLoader = torch::dataloader(data, batch_size=batch_size, shuffle=TRUE, num_workers=0, pin_memory=pin_memory, drop_last=TRUE)

        self$history = list()


        if(!is.null(weights)) {
          weights = torch::torch_tensor(weights, dtype = self$dtype, device = self$device)
        } else {
          weights = torch::torch_ones(dim(Y)[1], dtype = self$dtype, device = self$device)
        }

        cli::cli_progress_bar(format = "Epoch: {cli::pb_current}/{cli::pb_total} {cli::pb_bar} ETA: {cli::pb_eta} Loss: {bl}", total = epochs, clear = FALSE)
        if(is.null(year_sequence)) year_sequence = 1:X[[1]]$shape[2]

        for(epoch in 1:epochs){
          counter = 1
          coro::loop(for (b in DataLoader) {
            batch_loss = matrix(NA, nrow = 10000, ncol = 6L)
            x_mort = b[[1]]$to(device = self$device, non_blocking=TRUE)
            x_growth = b[[2]]$to(device = self$device, non_blocking=TRUE)
            x_reg = b[[3]]$to(device = self$device, non_blocking=TRUE)
            y = b[[4]]$to(device = self$device, non_blocking=TRUE)
            c = b[[5]]$to(device = self$device, non_blocking=TRUE)
            ind = b[[6]]$to(device = self$device, non_blocking=TRUE)
            dist = NULL
            if(!is.null(disturbance)) dist = b[[7]]$to(device = self$device, non_blocking=TRUE)
            self$optimizer$zero_grad()


            if (is.null(initCohort)) {
              cohorts = CohortMat$new(dims = c(x$shape[1], patches, self$sp),
                                      sp = self$sp,
                                      device=self$device)
              trees = cohorts$trees
              species = cohorts$species
              dbh = cohorts$dbh
            } else {
              trees = initCohort$trees$to(device = self$device, non_blocking=TRUE)[ind,]
              species = initCohort$species$to(device = self$device, non_blocking=TRUE)[ind,]
              dbh = initCohort$dbh$to(device = self$device, non_blocking=TRUE)[ind,]
            }

            pred_tmp = self$predict(dbh = dbh,
                                    trees = trees,
                                    species = species,
                                    env = list(x_mort, x_growth, x_reg),
                                    disturbance = dist,
                                    start_time = start_time,
                                    y = y,
                                    c = c,
                                    update_step = update_step,
                                    verbose = FALSE,
                                    weights = weights,
                                    year_sequence = year_sequence)
            pred = pred_tmp[[1]]
            loss = pred_tmp[[2]]

            self$optimizer$step()
            batch_loss[counter, ] =  as.numeric(loss$data()$cpu())
            counter <<- counter + 1

            self$param_history = c(self$param_history, list(lapply(self$parameters, function(p) as.matrix(p$cpu()))))
            #sf.scheduler.step()
          })
          bl = mean(rowSums(batch_loss), na.rm = TRUE)
          bl = round(bl, 5)
          #cat("Epoch: ", epoch, "Loss: ", bl, "\n")
          self$history[[epoch]] = colMeans(batch_loss, na.rm = TRUE)
          cli::cli_progress_update()

          if(!is.null( file )) {
            if(epoch %% thin == 0) saveRDS( self$param_history, file = file )
          }

        }
        cli::cli_progress_done()

        self$parameter_to_r()

        if(torch::cuda_is_available()) torch::cuda_empty_cache()

        # ignore debugging method
        pred$Patch = NULL
        pred$Cohort = NULL
        self$pred = list(long = pred2DF(list(Predictions = pred), "long"), wide = pred2DF(list(Predictions = pred), "wide"))
        self$year_sequence = year_sequence
        self$optimizer = NULL # TODO: best solution?
      }

    ))


