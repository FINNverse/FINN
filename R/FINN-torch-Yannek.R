
#' Calculate the basal area of a tree given the diameter at breast height (dbh).
#'
#' Calculate the basal area of a tree given the diameter at breast height (dbh).
#'
#' @param dbh (torch.Tensor) The diameter at breast height of the tree.
#'
#' @return (torch.Tensor) The basal area of the tree.
#'
#' @examples
#' dbh = torch.tensor(50)
#' basal_area = BA_P(dbh)
#' print(basal_area)
#'
#' @import torch
#' @export
BA_P = function(dbh) {
  return(pi*(dbh/100./2.)^2.0)
}


#' Calculate the height of a tree based on its diameter at breast height and a global parameter.
#'
#' This function calculates the height of a tree based on the diameter at breast height (dbh) and a parameter.
#'
#' @param dbh A numeric value representing the diameter at breast height of the tree.
#' @param parGlobal A numeric value representing the global parameter used in the height calculation.
#'
#' @return A numeric value representing the calculated height of the tree.
#'
#' @examples
#' height_P(30, 0.5)
#'
#' @export
height_P = function(dbh, parGlobal) {
  height = (exp((((dbh * parGlobal) / (dbh+100))))-1)*100 + 0.001
  return(height)
}


#' Compute the fraction of available light (AL) for each cohort based on the given parameters.
#'
#' @param dbh \code{torch.Tensor} Diameter at breast height for each cohort.
#' @param Species \code{torch.Tensor} Species index for each cohort.
#' @param nTree Number of trees.
#' @param parGlobal \code{torch.Tensor} Global parameters for all species.
#' @param h \code{torch.Tensor} (Optional) Height of each cohort. Defaults to \code{NULL}.
#' @param minLight \code{float} (Optional) Minimum light requirement. Defaults to 50.
#' @return \code{torch.Tensor} Fraction of available light (AL) for each cohort.
#' @import torch
#' @examples
#' compF_P(dbh = torch_tensor(c(10, 15, 20)), Species = torch_tensor(c(1, 2, 1)),
#'         nTree = 100, parGlobal = torch_tensor(c(0.3, 0.5)), h = torch_tensor(c(5, 7, 6)), minLight = 40)
compF_P = function(dbh, Species, nTree, parGlobal, h = NULL, minLight = 50.){

  ba = (BA_P(dbh)*nTree)/0.1
  cohortHeights = height_P(dbh, parGlobal[Species])$unsqueeze(3)
  if(is.null(h)) {
    h = cohortHeights
    BA_height = (ba$unsqueeze(3)*torch_sigmoid((cohortHeights - h$permute(c(1,2, 4, 3)) - 0.1)/1e-3) )$sum(-2) # AUFPASSEN
  }else{
    BA_height = (ba$unsqueeze(3)*torch.sigmoid((cohortHeights - 0.1)/1e-3))$sum(-2)
  }
  AL = 1.-BA_height/minLight
  AL = torch_clamp(AL, min = 0)
  return(AL)
}

# YANNEK
#' Calculate forest plot growth
#'
#' This function calculates forest plot growth based on specified parameters.
#'
#' @param dbh Diameter at breast height
#' @param Species Species of tree
#' @param parGrowth Growth parameters
#' @param parMort Mortality parameters
#' @param pred Predicted values
#' @param AL Accumulated Light
#'
#' @return A numeric value representing the forest plot growth
#'
#' @import torch
#' @importFrom torch nn functional
#' @importFrom torch nn functional torch_sigmoid
#' @importFrom torch nn functional torch_clamp
#'
#' @export
growthFP = function(dbh, Species, parGrowth, parMort, pred, AL){

  shade = torch_sigmoid((AL + (1-parGrowth[,1][Species]) - 1)/1e-1)
  environment = index_species(pred, Species)
  pred = (shade*environment)
  # growth = (1.- torch.pow(1.- pred,4.0)) * parGrowth[Species,1]
  growth = pred/2 * parGrowth[,2][Species] * ((parMort[,2][Species]-dbh/100) / parMort[,2][Species])^(2)
  # growth = parGrowth[Species,1]
  # return torch.nn.functional.softplus(growth)
  return(torch_clamp(growth, min = 0.0))
}

#' Calculate the regeneration of forest patches based on the input parameters.
#'
#' @param Species Species information.
#' @param parReg Regeneration parameters.
#' @param pred Prediction values.
#' @param AL Available light variable for calculation.
#'
#' @return Regeneration values for forest patches.
#'
#' @import torch
#' @importFrom torch torch_sigmoid
#'
regFP = function(Species, parReg, pred, AL) {
  regP = torch_sigmoid((AL + (1-parReg) - 1)/1e-3)
  environment = pred
  regeneration = sample_poisson_relaxed((regP*(environment[,NULL])$`repeat`(c(1, Species$shape[2], 1)) + 1e-20)) # TODO, check if exp or not?! lambda should be always positive!
  regeneration = regeneration + regeneration$round()$detach() - regeneration$detach()
  return(regeneration)
}

# This function generates random numbers from a uniform distribution
# with specified low and high values and size similar to np.random.uniform in Python

# Function definition
np_runif = function(low, high, size) {
  N = prod(size)  # Calculate total number of values to generate
  array(runif(N, low, high), dim = size)  # Generate random numbers and reshape into desired size
}
init_FINN = function(
    sp = self$sp,
    device = self$device,
    parGlobal = self$parGlobal,
    parGrowth = self$parGrowth,
    parMort = self$parMort,
    parReg = self$parReg,
    parGrowthEnv = self$parGrowthEnv,
    parMortEnv = self$parMortEnv,
    parRegEnv = self$parRegEnv,
    hidden_growth = self$hidden_growth,
    hidden_mort = self$hidden_mort,
    hidden_reg = self$hidden_reg
    ){
  self$sp = sp
  self$device = device
  self$parGlobal = parGlobal
  self$parGrowth = parGrowth
  self$parMort = parMort
  self$parReg = parReg
  self$parGrowthEnv = parGrowthEnv
  self$parMortEnv = parMortEnv
  self$parRegEnv = parRegEnv
  self$hidden_growth = hidden_growth
  self$hidden_mort = hidden_mort
  self$hidden_reg = hidden_reg

  self$device = torch_device(device)
  self$sp = sp
  self$pred = NULL
  self$env = env
  self$optimizer = NULL
  self$dtype = torch_float32()
  self$nnRegEnv = self$build_NN(input_shape=env, output_shape=sp, hidden=hidden_reg, activation="selu", bias=list(FALSE), dropout=-99, last_activation = "relu")
  self$nnRegEnv$to(self$device)
  self$nnGrowthEnv = self$build_NN(input_shape=env, output_shape=sp, hidden=hidden_growth, activation="selu", bias=list(FALSE), dropout=-99)
  self$nnGrowthEnv$to(self$device)


  self$nnMortEnv = self$build_NN(input_shape=env, output_shape=sp, hidden=hidden_mort, activation="selu", bias=list(FALSE), dropout=-99)
  self$nnMortEnv$to(self$device)

  if(!is.null(parGrowthEnv)) self$set_weights_nnGrowthEnv(parGrowthEnv)
  if(!is.null(parMortEnv)) self$set_weights_nnMortEnv(parMortEnv)
  if(!is.null(parRegEnv)) self$set_weights_nnRegEnv(parRegEnv)

  if(is.null(parGlobal)){
    self$parGlobal = torch_tensor(uniform(0.3, 0.7, size = self$sp), requires_grad=TRUE, dtype=torch_float32(), device=self$device)
  }else{
    parGlobal = parGlobal$reshape(-1)
    self$parGlobal = torch_tensor(parGlobal, requires_grad=TRUE, dtype=torch_float32(), device=self$device)
  }

  if(is.null(parGrowth)){
    first = uniform(0, 6, size = c(self$sp,1))
    second = uniform(0, 6, size = c(self$sp,1))
    self$parGrowth = torch_tensor(cbind(first, second, 1), requires_grad=TRUE, dtype=torch.float32(), device=self$device)
  }else{
    self$parGrowth = torch_tensor(parGrowth, requires_grad=TRUE, dtype=torch_float32(), device=self$device)
  }

  if(is.null(parMort)){
    first = uniform(0, 2, size = c(self$sp,1))
    second = uniform(0.1, 5, size = c(self$sp,1))
    self$parMort = torch.tensor(cbind(first, second, 1), requires_grad=TRUE, dtype=torch_float32, device=self$device)
    # self._parMort = torch.tensor(np.random.uniform(0, 500, size = [self.sp,2]), requires_grad=True, dtype=torch.float32, device=self.device)
  }else{
    self$parMort = torch.tensor(parMort, requires_grad=True, dtype=torch_float32(), device=self$device)
  }

  if(is.null(parReg)){
    self$parReg = torch_tensor(uniform(0, 1, size = self$sp), requires_grad=True, dtype=torch_float32(), device=self$device)
  }else{
    parReg = parReg$reshape(-1)
    self$parReg = torch_tensor(parReg, requires_grad=TRUE, dtype=torch_float32(), device=self$device)
  }

  if(which == "env"){
    self$parameters = list(self$nnRegEnv$parameters(), self$nnGrowthEnv$parameters(), self$nnMortEnv$parameters())
  }else if(which == "species"){
    self$parameters = list(list(self$parGlobal), list(self$parGrowth), list(self$parMort), list(self$parReg))
  }else{
    self$parameters = list(list(self$parGlobal), list(self$parGrowth), list(self$parMort), list(self$parReg), list(self$nnRegEnv$parameters()), list(self$nnGrowthEnv$parameters()), list(self$nnMortEnv$parameters()))
  }

}
# input_shape = 2
# output_shape = 3
# hidden = c(1,3,2)
# bias = TRUE
# activation = "relu"
# dropout = 1
# last_activation ="sigmoid"

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
    for(i in seq_along(hidden)){}
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

  if(length(hidden) > 0){
    model_list = c(model_list, torch::nn_linear(hidden[length(hidden)], output_shape, bias=bias[length(hidden)]))
  }else{
    model_list = c(model_list, torch.nn.Linear(input_shape, output_shape, bias=FALSE))
  }
  if(last_activation == "sigmoid") model_list = c(model_list, torch::nn_sigmoid())
  if(last_activation == "relu") model_list = c(model_list, torch::nn_relu())
  return(do.call(torch::nn_sequential,model_list))
}

predict = function(self,
            dbh = NULL,
            nTree = NULL,
            Species = NULL,
            env = NULL,
            record_time = 0L,
            response ="dbh",
            pred_growth = NULL,
            pred_morth = NULL,
            pred_reg = NULL,
            patches = 50.,
            debug = TRUE){

    # Predicts the growth and mortality of trees based on the given inputs.
    #       Args:
    #           dbh (Optional[torch.Tensor]): The diameter at breast height of the trees. If None, it will be initialized using CohortMat.
    #           nTree (Optional[torch.Tensor]): The number of trees. If None, it will be initialized using CohortMat.
    #           Species (Optional[torch.Tensor]): The species of the trees. If None, it will be initialized using CohortMat.
    #           env (Optional[torch.Tensor]): The environmental data.
    #           record_time (int): The time at which to start recording the results.
    #           response (str): The response variable to use for aggregating results. Can be "dbh", "ba", or "ba_p".
    #           pred_growth (Optional[torch.Tensor]): The predicted growth values.
    #           pred_morth (Optional[torch.Tensor]): The predicted mortality values.
    #           pred_reg (Optional[torch.Tensor]): The predicted regeneration values.
    #           patches (Optional[float]): The number of patches.
    #
    #       Returns:
    #           Tuple[torch.Tensor, torch.Tensor]: The predicted dbh and number of trees for the recorded times points.

  if(is.null(dbh)){
    cohorts = CohortMat(dims = list(env.shape[1], patches, self$sp), sp = self$sp, device=self$device)
    nTree = cohorts$nTree
    Species = cohorts$Species
    dbh = cohorts$dbh
  }

  env = torch_tensor(env, dtype=self$dtype, device=self$device)

  # Predict env niches for all sites and timesteps
  if(is.na(pred_growth)) pred_growth = self$nnGrowthEnv(env)
  if(is.na(pred_morth)) pred_morth = self$nnMortEnv(env)
  if(is.na(pred_reg)) pred_reg = self$nnRegEnv(env)

  dbh = torch_tensor(dbh, dtype=self$dtype, device=self$device)

  nTree = torch_tensor(nTree, dtype=self$dtype, device=self$device)

  Species = torch_tensor(Species, dtype=torch$int64, device=self$device)

  cohort_ids = torch_randint(0, 50000, size=Species$shape)

  # Result arrays
  ## dbh / ba / ba*nTRee
  Result = torch_zeros(list(env$shape[1],env$shape[2],  dbh$shape[3]), device=self$device) # sites, time, species
  ## number of trees
  Result_T = torch_zeros(list(env$shape[1],env$shape[2],  dbh$shape[3]), device=self$device)

  AL = torch_zeros(list(env$shape[1], env$shape[2],  dbh$shape[3]), device=self$device)
  g = torch_zeros(list(env$shape[1], env$shape[2], dbh$shape[3]), device=self$device)
  m = torch_zeros(list(env$shape[1], env$shape[2], dbh$shape[3]), device=self$device)
  r = torch_zeros(list(env$shape[1], env$shape[2], dbh$shape[3]), device=self$device)


  if(debug){
    Result = lapply(1:6,torch_zeros(list(env$shape[1],env$shape[2],  dbh$shape[3]), device=self$device))
    Raw_results = list()
    Raw_cohorts = list()
  }else{
    Result = lapply(1:2,torch_zeros(list(env$shape[1],env$shape[2],  dbh$shape[3]), device=self$device))
    }

  time = env$shape[2]
  patches = dbh$shape[2]
  sp = self$sp

  # Run over timesteps
  for(i in 1:time){
      # aggregate results
      sites = env$shape[1]
    if(i > record_time){
      if(dbh$shape[3] != 0){
        if(response == "dbh"){
          ba = dbh
        }else if(response == "ba"){
          ba = BA_P(dbh)
        }else{
          ba = BA_T_P(dbh, nTree)
          # ba * nTree
        }

        labels = Species
        samples = list()
        samples = c(samples,(ba * torch_sigmoid((nTree - 0.5)/1e-3)))
        samples = c(samples, (nTree * torch_sigmoid((nTree - 0.5)/1e-3)))

        torch_zeros_like()
        Results_tmp = replicate(length(samples),torch_zeros_like(Result[[0]][,i,]))
        tmp_res = aggregate_results(labels, samples, Results_tmp)
        for(v in seq_along(2)){
          Result[[v]][,i,] = Result[[v]][,i,] + tmp_res[[v]]/patches
        }
      }
    }
  }
    # Model

  torch::with_no_grad({
    nTree_clone =  nTree$clone() #torch.zeros_like(nTree).to(self.device)+0
  })

  if(dbh$shape[3] != 0){
    AL = compF_P(dbh, Species, nTree_clone, self$parGlobal)
    g = growthFP(dbh, Species, self$parGlobal, self$parGrowth, self$parMort, pred_growth[,i,], AL)
    dbh = dbh+g
    AL = compF_P(dbh, Species, nTree_clone, self$parGlobal)

    m = mortFP(dbh, g, Species, nTree, self$parGlobal, self$parMort, pred_morth[,i,], AL) #.unsqueeze(3)
    nTree = torch.clamp(nTree - m, min = 0)
  }


  torch::with_no_grad({
    nTree_clone =  nTree$clone() #torch.zeros_like(nTree).to(self.device)+0
  })

  AL_reg = compF_P(dbh, Species, nTree_clone, self$parGlobal, h = torch_zeros(list(1, 1)))

  r = regFP(dbh, Species, self$parGlobal, self$parReg, pred_reg[,i,], AL_reg)
  # New recruits
  new_dbh = ((r-1+0.1)/1e-3)$sigmoid() # TODO: check!!! --> when r 0 dann dbh = 0, ansonsten dbh = 1 dbh[r==0] = 0
  new_nTree = r
  new_Species = torch_arange(0, sp, dtype=torch_int64(), device = self$device)$unsqueeze(1)$`repeat`(c(r$shape[1], r$shape[2], 1))
  new_cohort_id = torch_randint(0, 50000, size = list(sp, 1))$unsqueeze(1)$`repeat`(c(r$shape[1], r$shape[2], 1))


  if(debug){
    if(dbh$shape[3] != 0){
      labels = Species
      }
    }

  samples = list()
  samples = c(samples, AL)
  samples = c(samples, g)
  samples = c(samples, m)
  # samples = c(samples, r)

  # count number of cohorts
  Sp_tmp = Species$type(g$dtype)
  cohort_counts = aggregate_results(labels, list((Sp_tmp+1)/(Sp_tmp+1)), list(torch_zeros_like(Result[[1]][,i,], dtype=g$dtype)))

  Results_tmp = replicate(length(samples), torch_zeros_like(Result[[1]][,i,]))
  tmp_res = aggregate_results(labels, samples, Results_tmp)
  for(v in c(3,4,5)){
    Result[v][,i,] = Result[[v]][,i,] + tmp_res[[v-2]]/cohort_counts[[1]]/patches
  }

  # cohort ids
  Raw_cohorts = c(Raw_cohorts, cohort_ids)

  # reg extra
  tmp_res = aggregate_results(new_Species, list(r), list(torch.zeros(Result[[1]][,i,]$shape[1], sp )))
  Result[[6]][,i,] = Result[[6]][,i,] + tmp_res[[1]]/cohort_counts[1]/patches

  Raw_results = c(Raw_results, as.matrix(Species$cpu()$data), as.matrix(dbh$cpu()$data), as.matrix(m$cpu()$data), as.matrix(g$cpu()$data), as.matrix(r$cpu()$data))

  # Combine
  dbh = torch_cat(list(dbh, new_dbh), 3)
  nTree = torch_cat(list(nTree, new_nTree), 3)
  Species = torch_cat(list(Species, new_Species), 3)
  cohort_ids = torch_cat(list(cohort_ids, new_cohort_id), 3)

  # Pad tensors, expensive
  if(i %% 10 == 0){
    indices = (nTree > 0.5)$flatten(1, 2)
  }

  org_dim = Species$shape[1:3]
  org_dim_t = torch_tensor(org_dim, dtype = torch_long(), device = "cpu")
  dbh = pad_tensors_speed_up(dbh, indices, org_dim_t)$unflatten(1, org_dim)#$unsqueeze(3)
  nTree = pad_tensors_speed_up(nTree, indices, org_dim_t)$unflatten(1, org_dim)#$unsqueeze(3)
  cohort_ids = pad_tensors_speed_up(cohort_ids, indices, org_dim_t)$unflatten(1, org_dim)
  Species = pad_tensors_speed_up(Species, indices, org_dim_t)$unflatten(1, org_dim)


  if(debug){
    Result = c(Results, Raw_results)
    Result = c(Results, Raw_cohorts)
  }

  return(Result)
}


fit = function(self,
        X = NULL,
        Y = NULL,
        initCohort = NULL,
        epochs = 2L,
        batch_size = 20L,
        learning_rate = 0.1,
        start_time = 0.5,
        patches = 50L,
        response = "dbh"){

    # Fits the model to the given data.
    #
    #       Args:
    #           X (Optional[torch.Tensor]): The input data. Default is None.
    #           Y (Optional[torch.Tensor]): The target data. Default is None.
    #           epochs (int): The number of epochs to train the model. Default is 2.
    #           batch_size (int): The batch size for training. Default is 20.
    #           learning_rate (float): The learning rate for the optimizer. Default is 0.1.
    #           start_time (float): The start time for prediction. Default is 0.5.
    #           patches (int): The number of patches. Default is 50.
    #           response (str): The response variable. Default is "dbh", other options are 'ba' or 'ba*nT'.
    #
    #       Returns:
    #           None

  torch::optim_adagrad()
  if(is.null(self.optimizer)){
    self$optimizer = optim_adagrad(params = self.parameters, lr = learning_rate) # AdaBound was also good
    # TODO scheduler implementieren
    # self.scheduler = torch.optim.lr_scheduler.ExponentialLR(self.optimizer, gamma=0.9)
  }
  time = X$shape[2]
  start_time = round(start_time*time)

  stepSize = floor(X$shape[1] / batch_size) # type: ignore

  if(self$device$type == 'cuda'){
    torch.cuda.set_device(self$device)
    pin_memory = FALSE
  }else{
    pin_memory = True
  }
  Counts = round(Y[,,,2])
  indices = torch_arange(1, X.shape[1])
  data = tensor_dataset(torch_tensor(X, dtype=self$dtype, device=torch$device('cpu')),
                        torch_tensor(Y[,,,], dtype=self$dtype, device=torch$device('cpu')),
                        torch_tensor(Counts, dtype=torch_int64(), device=torch$device('cpu')),
                        torch_tensor(indices, dtype=torch_int64(), device=torch$device('cpu'))
  )
  DataLoader = torch::dataloader(data, batch_size=batch_size, shuffle=TRUE, num_workers=0, pin_memory=pin_memory, drop_last=TRUE)

  batch_loss = torch_zeros(stepSize)
  self$history = torch_zeros(epochs)


  #### bis hier habe ich gemacht
  ep_bar = tqdm(range(epochs),bar_format= "Iter: {n_fmt}/{total_fmt} {l_bar}{bar}| [{elapsed}, {rate_fmt}{postfix}]", file=sys.stdout)
  for(epoch in ep_bar){

    for step, (x, y, c, ind) in enumerate(DataLoader):
    self.optimizer.zero_grad()

      if initCohort is None:
        cohorts = CohortMat(dims = [x.shape[0], patches, self.sp], sp = self.sp, device=self.device)
      nTree = cohorts.nTree
      Species = cohorts.Species
      dbh = cohorts.dbh
      else:
        nTree = (initCohort.nTree[ind,...]).to(self.device, non_blocking=True)
      Species = (initCohort.Species[ind,...]).to(self.device, non_blocking=True)
      dbh = (initCohort.dbh[ind,...]).to(self.device, non_blocking=True)

      x = x.to(self.device, non_blocking=True)
      y = y.to(self.device, non_blocking=True)
      c = c.to(self.device, non_blocking=True)
      pred = self.predict(dbh, nTree, Species, x, start_time, response)
      loss1 = torch.nn.functional.mse_loss(y[:, start_time:,:,0], pred[0][:,start_time:,:]).mean() # dbh / ba
      #loss2 = torch.nn.functional.mse_loss(y[:, start_time:,:,1], pred[1][:,start_time:,:]).mean() # nTree
      loss2 = torch.distributions.Poisson(pred[1][:,start_time:,:]+0.001).log_prob(c[:, start_time:,:]).mean().negative() # nTree
      loss = (loss1 + loss2)
      #loss = torch.nn.functional.mse_loss((y[:, start_time:,:,0]+1).log(), (pred[0][:,start_time:,:]*pred[1][:,start_time:,:] + 1).log()).mean() # dbh / ba

      loss.backward()
      self.optimizer.step()
      batch_loss[step] = loss.item()
      #sf.scheduler.step()
      bl = np.mean(batch_loss)
      bl = np.round(bl, 3)
      ep_bar.set_postfix(loss=f'{bl}')
      self.history[epoch] = bl
  }
  torch.cuda.empty_cache()
  self.pred = pred
}

library(R6)

FINN = R6Class(
  classname = 'FINN',
  public = list(
    # helper functions
    np_runif = np_runif, # a random uniform function similar to np.random.uniform
    # default model parameters
    sp = NULL,
    device = 'cpu',
    parGlobal = NULL, # must be dim [species]
    parGrowth = NULL, # must be dim [species, 2], first for shade tolerance
    parMort = NULL, # must be dim [species, 2], first for shade tolerance,
    parReg = NULL, # must be dim [species]
    parGrowthEnv = NULL, # must be dim [species, 2], first for shade tolerance
    parMortEnv = NULL, # must be dim [species, 2], first for shade tolerance,
    parRegEnv = NULL, # must be dim [species]
    hidden_growth = list(),
    hidden_mort = list(),
    hidden_reg = list(),
    #initialisation function
    initialize = init_FINN,
    build_NN = build_NN,
    predict = predict

    ))

obj = FINN$new(sp = 5, device = "egal") obj$sp

class FINN:
  def __init__(self,
               device: str='cpu',  # 'mps'
               sp: int=5,
               env: int=2,
               which: str="both",
               parGlobal: Optional[np.ndarray]=None, # must be dim [species]
               parGrowth: Optional[np.ndarray]=None, # must be dim [species, 2], first for shade tolerance
               parMort: Optional[np.ndarray]=None, # must be dim [species, 2], first for shade tolerance,
               parReg: Optional[np.ndarray]=None, # must be dim [species]
               parGrowthEnv: Optional[np.ndarray]=None, # must be dim [species, 2], first for shade tolerance
               parMortEnv: Optional[np.ndarray]=None, # must be dim [species, 2], first for shade tolerance,
               parRegEnv: Optional[np.ndarray]=None, # must be dim [species]
               hidden_growth: List[int] = [],
               hidden_mort: List[int]  = [],
               hidden_reg: List[int]  = []
  ):
  # """Initialize the model.
  #
  #       Args:
  #           device (str, optional): The device to use for computation. Supported options are 'cpu', 'cuda', and 'mps'. Defaults to 'cpu'.
  #           sp (int, optional): The number of species. Defaults to 5.
  #           env (int, optional): The number of environmental covariates. Defaults to 2.
  #           parGlobal (Optional[np.ndarray], optional): The global parameters. Must be of dimension [species]. Defaults to None.
  #           parGrowth (Optional[np.ndarray], optional): The growth parameters. Must be of dimension [species, 2], with the first column representing shade tolerance. Defaults to None.
  #           parMort (Optional[np.ndarray], optional): The mortality parameters. Must be of dimension [species, 2], with the first column representing shade tolerance. Defaults to None.
  #           parReg (Optional[np.ndarray], optional): The regeneration parameters. Must be of dimension [species]. Defaults to None.
  #           hidden_growth (List[int], optional): The hidden layers for the growth neural network. Defaults to [].
  #           hidden_mort (List[int], optional): The hidden layers for the mortality neural network. Defaults to [].
  #           hidden_reg (List[int], optional): The hidden layers for the regeneration neural network. Defaults to [].
  #
  #       Returns:
  #           None
  #       """





def continue_fit(self,
                 X: Optional[torch.Tensor]=None,
                 Y: Optional[torch.Tensor]=None,
                 initCohort: CohortMat = None,
                 epochs: int=2,
                 batch_size: int=20,
                 learning_rate: float=0.1,
                 start_time: float=0.5,
                 patches: int=50,
                 response: str="dbh"):

  """Continues the training of the model using the provided data.

        Args:
            X (Optional[torch.Tensor]): The input data. Defaults to None.
            Y (Optional[torch.Tensor]): The target data. Defaults to None.
            epochs (int): The number of training epochs. Defaults to 2.
            batch_size (int): The batch size for training. Defaults to 20.
            learning_rate (float): The learning rate for training. Defaults to 0.1.
            start_time (float): The starting time for training. Defaults to 0.5.
            patches (int): The number of patches. Defaults to 50.
            response (str): The response type. Defaults to "dbh".

        Returns:
            None

        Note:
            This function internally calls the `fit` method to continue the training process.
        """

self.fit(X,
         Y,
         initCohort,
         epochs,
         batch_size,
         learning_rate,
         start_time,
         patches,
         response)

def gradients(self):
  return [p.grad for p in self.optimizer.param_groups[0]["params"]]

@property
def parGrowth(self):
  return self._parGrowth.cpu().data.numpy()

@property
def parMort(self):
  return self._parMort.cpu().data.numpy()

@property
def parReg(self):
  return self._parReg.cpu().data.numpy()

@property
def parGlobal(self):
  return self._parGlobal.cpu().data.numpy()

@property
def GrowthEnv(self):
  return [(lambda p: p.data.cpu().numpy())(p) for p in self.nnGrowthEnv.parameters()]

@property
def MortEnv(self):
  return [(lambda p: p.data.cpu().numpy())(p) for p in self.nnMortEnv.parameters()]

@property
def RegEnv(self):
  return [(lambda p: p.data.cpu().numpy())(p) for p in self.nnRegEnv.parameters()]

@property
def weights(self):
  return [(lambda p: p.data.cpu().numpy())(p) for p in self.parameters]


def set_weights_nnMortEnv(self, w: List[np.ndarray]):
  """Set weights for the neural network regularization environment.

            Args:
                w (List[np.ndarray]): List of numpy arrays containing weights and biases.

            Sets the weights and biases of the linear layers in the neural network regularization environment
            based on the provided numpy arrays. The function iterates through the layers of the neural network
            regularization environment and assigns the weights and biases accordingly.

            Note:
                This function modifies the weights and biases of the neural network regularization environment in-place.

            Returns:
                None
        """
with torch.no_grad():
  counter = 0
for i in range(len(self.nnMortEnv)):
  if type(self.nnMortEnv[i]) is torch.nn.modules.linear.Linear:
  self.nnMortEnv[i].weight = torch.nn.Parameter(torch.tensor(w[counter], dtype=self.nnMortEnv[i].weight.dtype, device=self.nnMortEnv[i].weight.device))
counter+=1
if self.nnMortEnv[i].bias is not None:
  self.nnMortEnv[i].bias = torch.nn.Parameter(torch.tensor(w[counter], dtype=self.nnMortEnv[i].bias.dtype, device=self.nnMortEnv[i].bias.device))
counter+=1

def set_weights_nnGrowthEnv(self, w: List[np.ndarray]):
  """Set weights for the neural network regularization environment.

            Args:
                w (List[np.ndarray]): List of numpy arrays containing weights and biases.

            Sets the weights and biases of the linear layers in the neural network regularization environment
            based on the provided numpy arrays. The function iterates through the layers of the neural network
            regularization environment and assigns the weights and biases accordingly.

            Note:
                This function modifies the weights and biases of the neural network regularization environment in-place.

            Returns:
                None
        """
with torch.no_grad():
  counter = 0
for i in range(len(self.nnGrowthEnv)):
  if type(self.nnGrowthEnv[i]) is torch.nn.modules.linear.Linear:
  self.nnGrowthEnv[i].weight = torch.nn.Parameter(torch.tensor(w[counter], dtype=self.nnGrowthEnv[i].weight.dtype, device=self.nnGrowthEnv[i].weight.device))
counter+=1
if self.nnGrowthEnv[i].bias is not None:
  self.nnGrowthEnv[i].bias = torch.nn.Parameter(torch.tensor(w[counter], dtype=self.nnGrowthEnv[i].bias.dtype, device=self.nnGrowthEnv[i].bias.device))
counter+=1

def set_weights_nnRegEnv(self, w: List[np.ndarray]):
  """Set weights for the neural network regularization environment.

            Args:
                w (List[np.ndarray]): List of numpy arrays containing weights and biases.

            Sets the weights and biases of the linear layers in the neural network regularization environment
            based on the provided numpy arrays. The function iterates through the layers of the neural network
            regularization environment and assigns the weights and biases accordingly.

            Note:
                This function modifies the weights and biases of the neural network regularization environment in-place.

            Returns:
                None
        """

with torch.no_grad():
  counter = 0
for i in range(len(self.nnRegEnv)):
  if type(self.nnRegEnv[i]) is torch.nn.modules.linear.Linear:
  self.nnRegEnv[i].weight = torch.nn.Parameter(torch.tensor(w[counter], dtype=self.nnRegEnv[i].weight.dtype, device=self.nnRegEnv[i].weight.device))
counter+=1
if self.nnRegEnv[i].bias is not None:
  self.nnRegEnv[i].bias = torch.nn.Parameter(torch.tensor(w[counter], dtype=self.nnRegEnv[i].bias.dtype, device=self.nnRegEnv[i].bias.device))
counter+=1


