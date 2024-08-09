#' @include getter_setter.R
NULL

#' Calculate the BA_stemsal area of a tree given the diameter at breast height (dbh).
#'
#' This function calculates the BA_stemsal area of a tree given the diameter at breast height (dbh).
#'
#' @param dbh torch.Tensor The diameter at breast height of the tree.
#'
#' @return torch.Tensor The BA_stemsal area of the tree.
#'
#' @examples
#' dbh = torch::torch_tensor(50)
#' BA_stemsal_area = BA_stem(dbh)
#' print(BA_stemsal_area)
#'
#' @import torch
#' @export
BA_stem = function(dbh) {
  return(pi*(dbh/100./2.)^2.0)
}


#' BA_stemsal aera multiplied with number of species
#'
#' @param dbh dbh
#' @param trees number of Trees
#'
#' @export
BA_stand = function(dbh, trees, patch_size_ha) {
  return((pi*(dbh/100./2.)$pow(2.0)*trees)/patch_size_ha)
}



#' Calculate the height of a tree BA_stemsed on its diameter at breast height and a global parameter.
#'
#' This function calculates the height of a tree BA_stemsed on the diameter at breast height (dbh) and a parameter.
#'
#' @param dbh A numeric value representing the diameter at breast height of the tree.
#' @param parHeight A numeric value representing the global parameter used in the height calculation.
#'
#' @return A numeric value representing the calculated height of the tree.
#'
#' @examples
#' height(30, 0.5)
#'
#' @export
height = function(dbh, parHeight) {
  height = (exp((((dbh * parHeight) / (dbh+100))))-1)*100 + 0.001
  return(height)
}



#' Compute the fraction of available light (light) for each cohort BA_stemsed on the given parameters.
#'
#' This function calculates the fraction of available light for each cohort of trees BA_stemsed on their diameter at breast height (dbh), species, number of trees, and global parameters.
#'
#' @param dbh torch.Tensor Diameter at breast height for each cohort.
#' @param species torch.Tensor species index for each cohort.
#' @param trees Number of trees.
#' @param parHeight torch.Tensor global parameters for all species.
#' @param h torch.Tensor (Optional) Height of each cohort. Defaults to NULL.
#' @param minLight float (Optional) Minimum light requirement. Defaults to 50.
#'
#' @return torch.Tensor Fraction of available light (light) for each cohort.
#' @import torch
#' @examples
#' competition(dbh = torch::torch_tensor(c(10, 15, 20)), species = torch::torch_tensor(c(1, 2, 1)),
#'         trees = 100, parHeight = torch::torch_tensor(c(0.3, 0.5)), h = torch::torch_tensor(c(5, 7, 6)), minLight = 40)
#' @export
competition = function(dbh, species, trees, parHeight, h = NULL, minLight = 50., patch_size_ha, ba = NULL, cohortHeights = NULL){


  if(is.null(ba)) ba = BA_stand(dbh = dbh, trees = trees, patch_size_ha = patch_size_ha)
  if(is.null(cohortHeights)) cohortHeights = height(dbh, parHeight[species])$unsqueeze(4)
  if(is.null(h)) {
    h = cohortHeights

    # TODO:
    # Rework the height comparison, the sigmoid function preserves the gradients but they are not great (either -1 or 1)

    ba_height = (ba$unsqueeze_(4)$multiply(((cohortHeights - h$permute(c(1,2, 4, 3)) - 0.1)/1e-2)$sigmoid_() ))$sum(-2) # AUFPASSEN
  }else{
    ba_height = (ba$unsqueeze_(4)$multiply_(((cohortHeights - 0.1)/1e-2)$sigmoid_() ))$sum(-2)
  }
  light = 1.-ba_height/minLight
  light = torch_clamp(light, min = 0)
  return(light)
}



#' Mortality
#'
#' @param dbh dbh
#' @param species species
#' @param trees trees
#' @param parMort parMort
#' @param pred predictions
#' @param light available light
#'
#' @export
mortality = function(dbh, species, trees, parMort, pred, light, debug = F) {
  # TODO remove constant part
  shade = 1-torch_sigmoid((light + (1-parMort[,1][species]) - 1)/(1/10^(1.5 + torch_abs(light-0.5))))
  environment = index_species(pred, species)
  gPSize = torch_clamp(0.1*(dbh/torch_clamp((parMort[,2][species]*100), min = 0.00001))$pow(2.3), max = 1.0)
  # gPSize = torch_sigmoid(gPSize)
  # TODO
  # clamp can lead to vanishing gradients, sigmoid is not perfect but probably better here!
  # shade -> [0,1]
  # gPsize -> [0, 1]
  # environment -> [0, 1]
  # -> raw pred -> [0, 5]
  predM = torch_clamp((environment*(shade+gPSize) + shade*gPSize + shade + gPSize), min = 0.0001, max = 0.9999)
  # predM = torch_sigmoid((environment*(shade+gPSize) + shade*gPSize + shade + gPSize -1.5)*2 )
  mort1 = binomial_from_gamma(trees+trees$le(0.5)$float(), predM)*trees$ge(0.5)$float()
  mort2 = mort1 + mort1$round()$detach() - mort1$detach()
  if(debug == TRUE) out = list(shade = shade, light = light, environment = environment, gPSize = gPSize, predM = predM, mort1 = mort1, mort2 = mort2) else out = mort2
  return(out)
}



#' Calculate growth
#'
#' This function calculates growth BA_stemsed on specified parameters.
#'
#' @param dbh torch.Tensor Diameter at breast height.
#' @param species torch.Tensor species of tree.
#' @param parGrowth torch.Tensor Growth parameters.
#' @param parMort torch.Tensor Mortality parameters.
#' @param pred torch.Tensor Predicted values.
#' @param light torch.Tensor Accumulated Light.
#'
#' @return torch.Tensor A tensor representing the forest plot growth.
#'
#' @import torch
#'
#' @export
growth = function(dbh, species, parGrowth, parMort, pred, light){

  shade = torch_sigmoid((light + (1-parGrowth[,1][species]) - 1)/1e-1)
  environment = index_species(pred, species)
  pred = (shade*environment)
  # growth = (1.- torch.pow(1.- pred,4.0)) * parGrowth[species,1]
  growth = pred/2 * parGrowth[,2][species] * ((parMort[,2][species]-dbh/100) / parMort[,2][species])^(2)
  # growth = parGrowth[species,1]
  # return torch.nn.functional.softplus(growth)
  return(torch_clamp(growth, min = 0.0))
}

#' Calculate the regeneration of forest patches BA_stemsed on the input parameters.
#'
#' This function calculates the regeneration of forest patches BA_stemsed on species information, regeneration parameters, prediction values, and available light.
#'
#' @param species torch.Tensor species information.
#' @param parReg torch.Tensor Regeneration parameters. 0 <= parReg <= 1
#' This parameter denotes the fraction of light needed for a species to regenerate.
#' In general low values for high regeneration and high values for low regeneration.
#' @param pred torch.Tensor Prediction values.
#' @param light torch.Tensor Available light variable for calculation.
#'
#' @return torch.Tensor Regeneration values for forest patches.
#'
#' @import torch
#' @importFrom torch torch_sigmoid
#' @export
regeneration = function(species, parReg, pred, light, patch_size_ha, debug = F) {
  if("matrix" %in% class(pred)) pred = torch::torch_tensor(pred)
  environment = pred
  regP = (1 / (1 + torch_exp(-10 * (light - parReg))) - 1 / (1 + torch_exp(10 * parReg))) / (1 - 1 / (1 + torch_exp(10 * (1 - parReg))))
  #regP = torch_sigmoid((light + (1-parReg) - 1)/1e-3) # TODO masking? better https://pytorch.org/docs/stable/generated/torch.masked_select.html
  mean = (regP*(environment[,NULL])$`repeat`(c(1, species$shape[2], 1))+0.2)
  regeneration1 = sample_poisson_gaussian(mean*patch_size_ha) # TODO, check if exp or not?! lambda should be always positive!
  regeneration2 = regeneration1 + regeneration1$round()$detach() - regeneration1$detach()
  if(debug == T) out = list(regP = regP, mean = mean, regeneration1 = regeneration1, regeneration2 = regeneration2) else out = regeneration2
  return(out)
}



#' Initialize the FINN model with specified parameters.
#'
#' This function initializes the FINN model with specified species, device, global parameters, growth parameters, mortality parameters, regeneration parameters, and neural network settings.
#'
#' @param sp integer Number of species.
#' @param device character Device to use ('cpu' or 'cuda').
#' @param parHeight torch.Tensor global parameters.
#' @param parGrowth torch.Tensor Growth parameters.
#' @param parMort torch.Tensor Mortality parameters.
#' @param parReg torch.Tensor Regeneration parameters.
#' @param parGrowthEnv torch.Tensor Growth environment parameters.
#' @param parMortEnv torch.Tensor Mortality environment parameters.
#' @param parRegEnv torch.Tensor Regeneration environment parameters.
#' @param hidden_growth list Hidden layers for growth model.
#' @param hidden_mort list Hidden layers for mortality model.
#' @param hidden_reg list Hidden layers for regeneration model.
init_FINN = function(
    sp = self$sp,
    env = self$env,
    device = self$device,
    parHeight = self$parHeight,
    parGrowth = self$parGrowth,
    parMort = self$parMort,
    parReg = self$parReg,
    parGrowthEnv = self$parGrowthEnv,
    parMortEnv = self$parMortEnv,
    parRegEnv = self$parRegEnv,
    hidden_growth = self$hidden_growth,
    hidden_mort = self$hidden_mort,
    hidden_reg = self$hidden_reg,
    bias = self$bias,
    patch_size_ha = self$patch_size_ha,
    minLight = self$minLight,
    disturbance = self$disturbance,
    which = c("all", "env", "species")
    ){

  # check input parameters
  if(any(parHeight < 0 | parHeight > 1)) stop("parHeight cannot be <0 or >1")
  if(any(parReg < 0 | parReg > 1)) stop("parReg cannot be <0 or >1")

  which = match.arg(which)

  self$sp = sp
  self$device = device
  self$parHeight = parHeight
  self$parGrowth = parGrowth
  self$parMort = parMort
  self$parReg = parReg
  self$parGrowthEnv = parGrowthEnv
  self$parMortEnv = parMortEnv
  self$parRegEnv = parRegEnv
  self$hidden_growth = hidden_growth
  self$hidden_mort = hidden_mort
  self$hidden_reg = hidden_reg
  self$patch_size_ha = patch_size_ha
  self$minLight = minLight

  self$device = torch_device(device)
  self$sp = sp
  self$pred = NULL
  self$env = env
  self$optimizer = NULL
  self$dtype = torch_float32()
  self$nnRegEnv = self$build_NN(input_shape=env, output_shape=sp, hidden=hidden_reg, activation="selu", bias=bias, dropout=-99, last_activation = "relu")
  self$nnGrowthEnv = self$build_NN(input_shape=env, output_shape=sp, hidden=hidden_growth, activation="selu", bias=bias, dropout=-99)
  self$nnMortEnv = self$build_NN(input_shape=env, output_shape=sp, hidden=hidden_mort, activation="selu", bias=bias, dropout=-99)
  if(!is.null(parGrowthEnv)) self$set_weights_nnGrowthEnv(parGrowthEnv)
  if(!is.null(parMortEnv)) self$set_weights_nnMortEnv(parMortEnv)
  if(!is.null(parRegEnv)) self$set_weights_nnRegEnv(parRegEnv)


  self$nnMortEnv$to(device = self$device)
  self$nnGrowthEnv$to(device = self$device)
  self$nnRegEnv$to(device = self$device)

  if(is.null(parHeight)){
    parHeight = np_runif(0.3, 0.7, size = self$sp)
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

  self$set_parHeight(parHeight)
  self$set_parGrowth(parGrowth)
  self$set_parMort(parMort)
  self$set_parReg(parReg)

  if(which == "env"){
    self$parameters = c(self$nnRegEnv$parameters, self$nnGrowthEnv$parameters, self$nnMortEnv$parameters)
    self$parHeight$requires_grad_(FALSE)
    self$parGrowth$requires_grad_(FALSE)
    self$parMort$requires_grad_(FALSE)
    self$parReg$requires_grad_(FALSE)
    names(self$parameters) = c("R_E", "G_E", "M_E")
  }else if(which == "species"){
    self$parameters = c(self$parHeight, self$parGrowth, self$parMort, self$parReg)
    names(self$parameters) = c("H", "G", "M", "R")
    .n = lapply(self$nnRegEnv$parameters, function(p) p$requires_grad_(FALSE))
    .n = lapply(self$nnGrowtEnv$parameters, function(p) p$requires_grad_(FALSE))
    .n = lapply(self$nnMortEnv$parameters, function(p) p$requires_grad_(FALSE))

  }else{
    self$parameters = c(self$parHeight, self$parGrowth, self$parMort,self$parReg, self$nnRegEnv$parameters, self$nnGrowthEnv$parameters, self$nnMortEnv$parameters)
    names(self$parameters) = c("H", "G", "M", "R","R_E", "G_E", "M_E" )
  }
  return(self) # Only for testing now
}

#' Build a neural network
#'
#' This function builds a neural network with specified input shape, output shape, hidden layers, bias, activation functions, dropout rate, and last activation function.
#'
#' @param input_shape integer Number of predictors.
#' @param output_shape integer Number of species.
#' @param hidden vector List of hidden layers.
#' @param bias vector Boolean values indicating whether to use bias in hidden layers.
#' @param activation vector List of activation functions.
#' @param dropout float Dropout rate.
#' @param last_activation character Last activation function.
#'
#' @return torch.nn.modules.container.Sequential Sequential neural network object.
#'
#' @examples
#' build_NN(input_shape = 2, output_shape = 3, hidden = c(1, 3, 2), bias = TRUE, activation = "relu", dropout = 1, last_activation = "sigmoid")
build_NN <- function(self,
              input_shape,
              output_shape,
              hidden, # vector
              bias, # vector
              activation, # vector
              dropout,
              last_activation ="sigmoid"){
    # Build neural network
    #
    # Args:
    # input_shape (int): Number of predictors
    # output_shape (int): Number of species
    # hidden (vector): List of hidden layers
    # bias (vector[bool]): Biases in hidden layers
    # activation (vector[str]): List of activation functions
    # dropout (float): Dropout rate
    #
    # Returns:
    # torch.nn.modules.container.Sequential: Sequential neural network object

  model_list = list()
  if(length(hidden) != length(activation)){
    activation = rep(activation[1], length(hidden))
  }

  if(length(bias) == 1){
    bias = rep(bias[1],length(hidden))
  }
  bias = c(FALSE, bias)

  if(length(hidden) > 0){
    for(i in 1:length(hidden)){
      if(i == 1){
        model_list = c(model_list,torch::nn_linear(input_shape, hidden[i], bias=bias[i]))
      }else{
        model_list = c(model_list,torch::nn_linear(hidden[i-1], hidden[i], bias=bias[i]))
      }
      if(activation[i] == "relu") model_list = c(model_list, torch::nn_relu())
      if(activation[i] == "selu") model_list = c(model_list, torch::nn_selu())
      if(activation[i] == "leakyrelu") model_list = c(model_list, torch::nn_leaky_relu())
      if(activation[i] == "tanh") model_list = c(model_list, torch::nn_tanh())
      if(activation[i] == "sigmoid") model_list = c(model_list, torch::nn_sigmoid())
      if(dropout > 0.0) model_list = c(model_list, torch::nn_dropout(p=dropout))
    }
  }

  if(length(hidden) > 0){
    model_list = c(model_list, torch::nn_linear(hidden[length(hidden)], output_shape, bias=bias[length(hidden)]))
  }else{
    model_list = c(model_list, torch::nn_linear(input_shape, output_shape, bias=FALSE))
  }
  if(last_activation == "sigmoid") model_list = c(model_list, torch::nn_sigmoid())
  if(last_activation == "relu") model_list = c(model_list, torch::nn_relu())
  return(do.call(torch::nn_sequential,model_list))
}


#' Predict the growth and mortality of trees BA_stemsed on the given inputs.
#'
#' This function predicts the growth and mortality of trees BA_stemsed on the given inputs, including diameter at breast height (dbh), number of trees, species, environmental data, and other parameters.
#'
#' @param dbh torch.Tensor (Optional) The diameter at breast height of the trees. If NULL, it will be initialized using CohortMat.
#' @param trees torch.Tensor (Optional) The number of trees. If NULL, it will be initialized using CohortMat.
#' @param species torch.Tensor (Optional) The species of the trees. If NULL, it will be initialized using CohortMat.
#' @param env torch.Tensor The environmental data.
#' @param start_time integer The time at which to start recording the results.
#' @param response character The response variable to use for aggregating results. Can be "dbh", "BA_stem", or "BA_stem".
#' @param pred_growth torch.Tensor (Optional) The predicted growth values.
#' @param pred_morth torch.Tensor (Optional) The predicted mortality values.
#' @param pred_reg torch.Tensor (Optional) The predicted regeneration values.
#' @param patches float (Optional) The number of patches.
#' @param debug logical Whether to run in debug mode.
#' @param y response tensor
#' @param c number of tree tensor
#' @param update_step backpropagation step length
#'
#' @return list A list of predicted values for dbh and number of trees for the recorded time points.
predict = function(
            dbh = NULL,
            trees = NULL,
            species = NULL,
            env = NULL,
            start_time = 1L,
            response ="dbh",
            pred_growth = NULL,
            pred_morth = NULL,
            pred_reg = NULL,
            patches = 50.,
            debug = TRUE,
            y = NULL,
            c = NULL,
            update_step = 1L){


  if(is.null(dbh)){
    cohorts = CohortMat$new(dims = c(env$shape[1], patches, self$sp), sp = self$sp, device=self$device)
    trees = cohorts$trees
    species = cohorts$species
    dbh = cohorts$dbh
  }

  env = torch_tensor(env, dtype=self$dtype, device=self$device)


  if(is.null(y)) {
    lapply(self$parameters, function(p) p$requires_grad_(FALSE) )
  }


  # Predict env niches for all sites and timesteps
  if(is.null(pred_growth)) pred_growth = self$nnGrowthEnv(env)
  if(is.null(pred_morth)) pred_morth = self$nnMortEnv(env)
  if(is.null(pred_reg)) pred_reg = self$nnRegEnv(env)


  dbh = torch_tensor(dbh, dtype=self$dtype, device=self$device)
  trees = torch_tensor(trees, dtype=self$dtype, device=self$device)
  species = torch_tensor(species, dtype=torch_int64(), device=self$device)
  cohort_ids = torch_randint(0, 50000, size=species$shape, device=self$device)
  light = torch_zeros(list(env$shape[1], env$shape[2],  dbh$shape[3]), device=self$device)
  g = torch_zeros(list(env$shape[1], env$shape[2], dbh$shape[3]), device=self$device)
  m = torch_zeros(list(env$shape[1], env$shape[2], dbh$shape[3]), device=self$device)
  r = torch_zeros(list(env$shape[1], env$shape[2], dbh$shape[3]), device=self$device)


  Result = lapply(1:6,function(tmp) torch_zeros(list(env$shape[1], env$shape[2], self$sp), device=self$device))

  # Can get memory intensive...
  if(debug) {
    Raw_results = list()
    Raw_cohorts = list()
  }

  time = env$shape[2]
  patches = dbh$shape[2]
  sp = self$sp
  sites = env$shape[1]

  predGrowthGlobal = self$nnGrowthEnv(env)
  predMortGlobal = self$nnMortEnv(env)
  predRegGlobal = self$nnRegEnv(env)
  predGrowthGlobal = torch::torch_split(predGrowthGlobal,split_size = 1, dim = 2L)
  predMortGlobal = torch::torch_split(predMortGlobal,split_size = 1, dim = 2L)
  predRegGlobal = torch::torch_split(predRegGlobal,split_size = 1, dim = 2L)


  # Run over timesteps
  for(i in 1:time){
    # cat("Time: ", i, "\n")

     pred_growth = self$nnGrowthEnv(env[,i,])
    # pred_growth = predGrowthGlobal[,i,]
     pred_morth = self$nnMortEnv(env[,i,])
    # pred_morth = predMortGlobal[,i,]
     pred_reg = self$nnRegEnv(env[,i,])
    # pred_reg = predRegGlobal[,i,]
    # pred_growth = predGrowthGlobal[[i]]$squeeze(dim = 2L)
    # pred_morth = predMortGlobal[[i]]$squeeze(dim = 2L)
    # pred_reg = predRegGlobal[[i]]$squeeze(dim = 2L)

    light = torch_zeros(list(env$shape[1], env$shape[2],  dbh$shape[3]), device=self$device)
    g = torch_zeros(list(env$shape[1], env$shape[2], dbh$shape[3]), device=self$device)
    m = torch_zeros(list(env$shape[1], env$shape[2], dbh$shape[3]), device=self$device)
    r = torch_zeros(list(env$shape[1], env$shape[2], dbh$shape[3]), device=self$device)

    dbh=dbh$detach()
    trees=trees$detach()
    species=species$detach()
    cohort_ids=cohort_ids$detach()

    # Model
    parHeight = self$get_parHeight()
    parMort = self$get_parMort()
    parGrowth = self$get_parGrowth()
    parReg = self$get_parReg()
    if(dbh$shape[3] > 0.5){

      # get Parameters parameter constrains


      light = competition(
        dbh = dbh,
        species = species,
        trees = trees,
        parHeight = parHeight,
        h = NULL,
        minLight = self$minLight,
        patch_size_ha = self$patch_size_ha,
        ba = NULL,
        cohortHeights = NULL
        )

      g = growth(
        dbh = dbh,
        species = species,
        parGrowth = parGrowth,
        parMort = parMort,
        pred = pred_growth,
        light = light
      )

      dbh = dbh + g

      light = competition(
        dbh = dbh,
        species = species,
        trees = trees,
        parHeight = parHeight,
        h = NULL,
        minLight = self$minLight,
        patch_size_ha = self$patch_size_ha,
        ba = NULL,
        cohortHeights = NULL
        )
      # cat("Second section  B2\n")
      m = mortality(
        dbh = dbh,
        species = species,
        trees = trees + 0.001,
        parMort = parMort,
        pred = pred_morth,
        light = light
      ) #.unsqueeze(3) # TODO check!
      #trees$sub_(m)$clamp_(min = 0.0)
      trees = torch_clamp(trees - m, min = 0) #### TODO if trees = 0 then NA...prevent!
    }


    AL_reg = competition( # must have dimension = n species in last dim
      dbh = dbh,
      species = species,
      trees = trees,
      parHeight = parHeight,
      h = 1,
      minLight = self$minLight,
      patch_size_ha = self$patch_size_ha,
      ba = NULL,
      cohortHeights = NULL
    )
    r = regeneration(species = species,
                     parReg = parReg,
                     pred = pred_reg,
                     light = AL_reg,
                     patch_size_ha = self$patch_size_ha)

    # New recruits
    new_dbh = ((r-1+0.1)/1e-3)$sigmoid() # TODO: check!!! --> when r 0 dann dbh = 0, ansonsten dbh = 1 dbh[r==0] = 0
    new_trees = r
    new_species = torch_arange(1, sp, dtype=torch_int64(), device = self$device)$unsqueeze(1)$`repeat`(c(r$shape[1], r$shape[2], 1))
    new_cohort_id = torch_randint(0, 50000, size = list(sp), device = self$device)$unsqueeze(1)$`repeat`(c(r$shape[1], r$shape[2], 1))

    if(dbh$shape[3] != 0){
      labels = species
      samples = vector("list", 3)
      samples[[1]] = light
      samples[[2]] = g
      samples[[3]] = m

      # count number of cohorts
      Sp_tmp = species$to(dtype = g$dtype)
      cohort_counts = aggregate_results(labels, list((Sp_tmp)/(Sp_tmp)), list(torch_zeros_like(Result[[1]][,i,], dtype=g$dtype, device = self$device)))

      Results_tmp = replicate(length(samples), torch_zeros_like(Result[[1]][,i,]))
      tmp_res = aggregate_results(labels, samples, Results_tmp)
      for(v in c(3,4,5)){
        Result[[v]][,i,] = Result[[v]][,i,]$add(tmp_res[[v-2]]/torch::torch_clamp(cohort_counts[[1]], min = 1.0)/patches) # TODO
      }

      # cohort ids
      if(debug) Raw_cohorts = c(Raw_cohorts, cohort_ids)

      # reg extra
      tmp_res = aggregate_results(new_species, list(r), list(torch::torch_zeros(Result[[1]][,i,]$shape[1], sp, device = self$device )))
      Result[[6]][,i,] = Result[[6]][,i,]$add(tmp_res[[1]]/torch::torch_clamp(cohort_counts[[1]], min = 1.0)/patches)


      if(debug) {
        Raw_results = c(Raw_results,list(list(list(torch::as_array(species$cpu())),
                                              list(torch::as_array(trees$cpu())),
                                              list(torch::as_array(dbh$cpu())),
                                              list(torch::as_array(m$cpu())),
                                              list(torch::as_array(g$cpu())),
                                              list(torch::as_array(r$cpu())))))
      }
    }

    # Combine
    dbh = torch::torch_cat(list(dbh, new_dbh), 3)
    trees = torch::torch_cat(list(trees, new_trees), 3)
    species = torch::torch_cat(list(species, new_species), 3)
    cohort_ids = torch::torch_cat(list(cohort_ids, new_cohort_id), 3)
    rm(new_dbh,new_trees,new_species,new_cohort_id )

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
     # })
    }


    # aggregate results
    # TODO: Position of the first block?
    if(i > 0){
      if(dbh$shape[3] != 0){
        if(response == "dbh"){
          BA_stem = dbh
        }else if(response == "BA_stem"){
          BA_stem = BA_stem(dbh)
        }else{
          BA_stem = BA_stand(dbh = dbh, trees = trees, patch_size_ha = self$patch_size_ha)
          # BA_stem * trees
        }
        labels = species
        samples = vector("list", 2)
        samples[[1]] = BA_stem * torch_sigmoid((trees - 0.5)/1e-3)
        samples[[2]] = trees * torch_sigmoid((trees - 0.5)/1e-3)
        Results_tmp = replicate(2, torch_zeros_like(Result[[1]][,i,]))

        tmp_res = aggregate_results(labels, samples, Results_tmp)
        # BA and number of trees Result[[1]] and Result[[2]]
        for(v in 1:2){
          Result[[v]][,i,] = Result[[v]][,i,]$add(tmp_res[[v]]$div_(patches))
        }
        rm(BA_stem)
      }
    }




    loss = torch_zeros(1L, device = self$device)
    if(i > 0 && dbh$shape[3] != 0 && !is.null(y) && (i %% update_step == 0)) {
      for(j in 1:6) {
        # 1 -> dbh/ba, 2 -> counts, 3 -> AL, 4 -> growth rates, 5 -> mort rates, 6 -> reg rates
        if(j != 2) {
          loss = loss+torch::nnf_mse_loss(y[, i,,j], Result[[j]][,i,])$mean()
        } else {
          loss = loss+torch::distr_poisson(Result[[2]][,i,]+0.001)$log_prob(c[, i,])$mean()$negative()
        }
      }
      loss$backward()
      # cat("Backprop....\n")
      for(j in 1:6) Result[[j]] = Result[[j]]$detach()
    }
    loss$detach_()
    # predGrowthGlobal[,i,] = predGrowthGlobal[,i,]$detach()
    # predMortGlobal[,i,] = predMortGlobal[,i,]$detach()
    # predRegGlobal[,i,] = predRegGlobal[,i,]$detach()

    # predGrowthGlobal[[i]] = predGrowthGlobal[[i]]$detach()
    # predMortGlobal[[i]] = predMortGlobal[[i]]$detach()
    # predRegGlobal[[i]] = predRegGlobal[[i]]$detach()

  }

  if(debug){
    Result = list(Predictions = c(Result, list(Raw_results, Raw_cohorts)), loss = loss)
  }


  if(is.null(y)) {
    lapply(self$parameters, function(p) p$requires_grad_(TRUE) )
  }

  return(Result)
}




#' Fit the model to the given data.
#'
#' This function fits the model to the given data using specified epochs, BA_stemtch size, learning rate, start time, and response variable.
#'
#' @param X torch.Tensor (Optional) The input data. Default is NULL.
#' @param Y torch.Tensor (Optional) The target data. Default is NULL.
#' @param initCohort list (Optional) Initial cohort data. Default is NULL.
#' @param epochs integer The number of epochs to train the model. Default is 2.
#' @param batch_size integer The BA_stemtch size for training. Default is 20.
#' @param learning_rate float The learning rate for the optimizer. Default is 0.1.
#' @param start_time float The start time for prediction. Default is 0.5.
#' @param patches integer The number of patches. Default is 50.
#' @param response character The response variable. Default is "dbh", other options are 'BA_stem' or 'BA_stem*nT'.
#' @param update_step backpropagation step length
#'
#' @return None
#'
#' @export
fit = function(
        X = NULL,
        Y = NULL,
        initCohort = NULL,
        epochs = 2L,
        batch_size = 20L,
        learning_rate = 0.1,
        start_time = 1,
        patches = 50L,
        response = "dbh",
        update_step = 1L){

  if(is.null(self$optimizer)){
    self$optimizer = optim_adagrad(params = self$parameters, lr = learning_rate) # AdaBound was also good
    # TODO scheduler implementieren
    # self.scheduler = torch.optim.lr_scheduler.ExponentialLR(self.optimizer, gamma=0.9)
  }
  time = X$shape[2]

  stepSize = floor(X$shape[1] / batch_size) # type: ignore

  if(self$device$type == 'cuda'){
    #torch.cuda.set_device(self$device)
    pin_memory = FALSE
  }else{
    pin_memory = TRUE
  }
  Counts = (Y[,,,2])$round()
  data = tensor_dataset(torch_tensor(X, dtype=self$dtype, device=torch_device('cpu')),
                        torch_tensor(Y[,,,], dtype=self$dtype, device=torch_device('cpu')),
                        torch_tensor(Counts, dtype=torch_int64(), device=torch_device('cpu')),
                        torch_arange(1, X$shape[1])$to(dtype = torch_int64() , device=torch_device('cpu'))
  )
  DataLoader = torch::dataloader(data, batch_size=batch_size, shuffle=TRUE, num_workers=0, pin_memory=pin_memory, drop_last=TRUE)

  self$history = torch::torch_zeros(epochs)

  cli::cli_progress_bar(format = "Epoch: {cli::pb_current}/{cli::pb_total} {cli::pb_bar} ETA: {cli::pb_eta} Loss: {bl}", total = epochs)

  for(epoch in 1:epochs){

    coro::loop(for (b in DataLoader) {
      batch_loss = c()
      x = b[[1]]$to(device = self$device, non_blocking=TRUE)
      y = b[[2]]$to(device = self$device, non_blocking=TRUE)
      c = b[[3]]$to(device = self$device, non_blocking=TRUE)
      ind = b[[4]]$to(device = self$device, non_blocking=TRUE)
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


      pred_tmp = self$predict( dbh, trees, species, x, start_time = start_time, response = response, y = y, c = c, update_step = update_step)
      pred = pred_tmp[[1]]
      loss = pred_tmp[[2]]

      self$optimizer$step()
      batch_loss =  c(batch_loss, loss$item())

      self$param_history = c(self$param_history, list(lapply(self$parameters, function(p) as.matrix(p$cpu()))))
      #sf.scheduler.step()
    })
    bl = mean(batch_loss)
    bl = round(bl, 3)
    #cat("Epoch: ", epoch, "Loss: ", bl, "\n")
    self$history[epoch] = bl
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  # torch::cuda_empty_cache()
  self$pred = pred
}





#' FINN: Forest Inventory Neural Network
#'
#' A class to initialize and train a Forest Inventory Neural Network (FINN) for predicting tree growth, mortality, and regeneration.
#'
#' @field np_runif Function to generate random numbers from a uniform distribution similar to np.random.uniform in Python.
#' @field sp integer Number of species.
#' @field device character Device to use ('cpu' or 'cuda').
#' @field parHeight torch.Tensor global parameters.
#' @field parGrowth torch.Tensor Growth parameters.
#' @field parMort torch.Tensor Mortality parameters.
#' @field parReg torch.Tensor Regeneration parameters.
#' @field parGrowthEnv torch.Tensor Growth environment parameters.
#' @field parMortEnv torch.Tensor Mortality environment parameters.
#' @field parRegEnv torch.Tensor Regeneration environment parameters.
#' @field hidden_growth list Hidden layers for growth model.
#' @field hidden_mort list Hidden layers for mortality model.
#' @field hidden_reg list Hidden layers for regeneration model.
#'
#' @export
FINN = R6::R6Class(
  classname = 'FINN',
  public = list(
    # helper functions
    np_runif = np_runif, # a random uniform function similar to np.random.uniform
    # default model parameters
    sp = NULL,
    env = NULL,
    dtype = NULL,
    parameters = NULL,
    optimizer = NULL,
    history = NULL,
    param_history = NULL,
    pred = NULL,
    device = 'cpu',
    parHeight = NULL, # must be dim [species]
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
    bias = FALSE,
    patch_size_ha = 0.1,
    minLight = 50,
    disturbance = 0.0,
    #initialisation function
    initialize = init_FINN,
    build_NN = build_NN,
    predict = predict,
    fit = fit,
    set_weights_nnGrowthEnv = set_weights_nnGrowthEnv,
    set_weights_nnMortEnv = set_weights_nnMortEnv,
    set_weights_nnRegEnv = set_weights_nnRegEnv,
    get_parMort = get_parMort,
    get_parGrowth = get_parGrowth,
    get_parHeight = get_parHeight,
    get_parReg = get_parReg,
    set_parMort = set_parMort,
    set_parGrowth = set_parGrowth,
    set_parHeight = set_parHeight,
    set_parReg = set_parReg
    ))


