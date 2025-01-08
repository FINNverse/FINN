#' FINNbase class
#'
#' @description
#' The `FINNBase` class provides core functionalities for building and managing neural networks within the FINN model framework. This class includes methods for constructing neural networks, generating random numbers, and managing model parameters.
FINNBase <- R6::R6Class(
  classname = "FINNbase",
  lock_objects = FALSE,
  lock_class = FALSE,
  public = list(

    #' @field parComp_r (`as.matrix(sp)`)\cr
    #' parComp parameters as R matrix, internal useage
    parComp_r = NULL,

    #' @field parGrowth_r (`as.matrix(sp, ..)`)\cr
    #' parGrowth parameters as R matrix, internal useage
    parGrowth_r = NULL,

    #' @field parMort_r (`as.matrix(sp, ...)`)\cr
    #' parMort parameters as R matrix, internal useage
    parMort_r = NULL,

    #' @field parReg_r (`as.vector(sp)`)\cr
    #' parReg parameters as R vector, internal useage
    parReg_r = NULL,

    #' @field parGrowthEnv_r (`as.list(matrix(), matrix(), ...)`)\cr
    #' parGrowthEnv parameters as R vector, internal useage
    parGrowthEnv_r = NULL,

    #' @field parMortEnv_r (`as.list(matrix(), matrix(), ...)`)\cr
    #' parMortEnv parameters as R vector, internal useage
    parMortEnv_r = NULL,

    #' @field parRegEnv_r (`as.list(matrix(), matrix(), ...)`)\cr
    #' parRegEnv parameters as R vector, internal useage
    parRegEnv_r = NULL,

    #' @description
    #' Build a neural network
    #'
    #' This function builds a neural network with specified input shape, output shape, hidden layers, bias, activation functions, dropout rate, and last activation function.
    #'
    #' @param input_shape integer. Number of predictors.
    #' @param output_shape integer. Number of species.
    #' @param hidden vector. List of hidden layers.
    #' @param bias vector. Boolean values indicating whether to use bias in hidden layers.
    #' @param activation vector. List of activation functions.
    #' @param dropout float. Dropout rate.
    #' @param last_activation character. Last activation function.
    #'
    #' @return torch.nn.modules.container.Sequential. Sequential neural network object.
    build_NN = function(input_shape,
                        output_shape,
                        hidden, # vector
                        bias, # vector
                        activation, # vector
                        dropout,
                        last_activation = "sigmoid") {

      model_list <- list()
      if (length(hidden) != length(activation)) {
        activation <- rep(activation[1], length(hidden))
      }

      if (length(bias) == 1) {
        bias <- rep(bias[1], length(hidden))
      }
      bias <- c(FALSE, bias)

      if (length(hidden) > 0) {
        for (i in 1:length(hidden)) {
          if (i == 1) {
            model_list <- c(model_list, torch::nn_linear(input_shape, hidden[i], bias = TRUE))
          } else {
            model_list <- c(model_list, torch::nn_linear(hidden[i - 1], hidden[i], bias = TRUE))
          }
          if (activation[i] == "relu") model_list <- c(model_list, torch::nn_relu())
          if (activation[i] == "selu") model_list <- c(model_list, torch::nn_selu())
          if (activation[i] == "leakyrelu") model_list <- c(model_list, torch::nn_leaky_relu())
          if (activation[i] == "tanh") model_list <- c(model_list, torch::nn_tanh())
          if (activation[i] == "sigmoid") model_list <- c(model_list, torch::nn_sigmoid())
          if (dropout > 0.0) model_list <- c(model_list, torch::nn_dropout(p = dropout))
        }
      }

      if (length(hidden) > 0) {
        model_list <- c(model_list, torch::nn_linear(hidden[length(hidden)], output_shape, bias = TRUE))
      } else {
        model_list <- c(model_list, torch::nn_linear(input_shape, output_shape, bias = FALSE))
      }
      if (last_activation == "sigmoid") model_list <- c(model_list, torch::nn_sigmoid())
      if (last_activation == "relu") model_list <- c(model_list, torch::nn_relu())
      return(do.call(torch::nn_sequential, model_list))
    },

    #' @description
    #' Generate random numbers from a uniform distribution
    #'
    #' This function generates random numbers from a uniform distribution with specified low and high values and size similar to `np.random.uniform` in Python.
    #'
    #' @param low numeric. Lower bound of the uniform distribution.
    #' @param high numeric. Upper bound of the uniform distribution.
    #' @param size numeric. Size of the output array.
    #'
    #' @return array. A numeric array of random numbers.
    np_runif = function(low, high, size) {
      N <- prod(size) # Calculate total number of values to generate
      array(runif(N, low, high), dim = size) # Generate random numbers and reshape into desired size
    },

    #' @description
    #' Set weights for the Growth Environment Neural Network
    #'
    #' This function assigns the specified weights to the growth environment neural network.
    #'
    #' @param weights list. A list of weight tensors to be set.
    #'
    #' @return None. The weights are set in the neural network.
    set_weights_nnGrowthEnv = function(weights) {
      self$nnGrowthEnv <- self$set_weights(weights, self$nnGrowthEnv)
    },

    #' @description
    #' Set weights for the Mortality Environment Neural Network
    #'
    #' This function assigns the specified weights to the mortality environment neural network.
    #'
    #' @param weights list. A list of weight tensors to be set.
    #'
    #' @return None. The weights are set in the neural network.
    set_weights_nnMortEnv = function(weights) {
      self$nnMortEnv <- self$set_weights(weights, self$nnMortEnv)
    },

    #' @description
    #' Set weights for the Regeneration Environment Neural Network
    #'
    #' This function assigns the specified weights to the regeneration environment neural network.
    #'
    #' @param weights list. A list of weight tensors to be set.
    #'
    #' @return None. The weights are set in the neural network.
    set_weights_nnRegEnv = function(weights) {
      self$nnRegEnv <- self$set_weights(weights, self$nnRegEnv)
    },

    #' @description
    #' Set weights for a specified Neural Network
    #'
    #' This function assigns the specified weights to a given neural network model.
    #'
    #' @param weights list. A list of weight tensors to be set.
    #' @param NN torch.nn.modules.container.Sequential. The neural network model to set weights for.
    #'
    #' @return torch.nn.modules.container.Sequential. The neural network with updated weights.
    set_weights = function(weights, NN) {
      torch::with_no_grad({
        counter <- 1
        for (i in 1:length(NN)) {
          if (inherits(NN[[i]], "nn_linear")) {
            NN$modules[[i]]$parameters$`0.weight`$set_data(weights[[counter]])
            counter <<- counter + 1
            if (!is.null(NN[[i]]$bias)) {
              NN$modules[[i]]$parameters$`0.bias`$set_data(weights[[counter]])
            }
          }
        }
      })
      return(NN)
    },

    #' @description
    #'
    #' Transform torch parameters to R (to save them)
    #'
    parameter_to_r = function() {

      self$parComp_r = to_r(self$parComp, TRUE)
      self$parGrowth_r = to_r(self$parGrowth)
      self$parMort_r = to_r(self$parMort)
      self$parReg_r = to_r(self$parReg, TRUE)
      self$parGrowthEnv_r = lapply(self$parGrowthEnv, function(p) to_r(p))
      self$parMortEnv_r = lapply(self$parMortEnv, function(p) to_r(p))
      self$parRegEnv_r = lapply(self$parRegEnv, function(p) to_r(p))
    },

    #' @description
    #'
    #' Check if all tensors are initialized, if not, reinitialize them. Necessary after restarting r session or after loading the object into the session using `readRDS(..)` or `load(...)`
    #'
    #' @param device can be used to specify a new device (useful for parallelization)
    #'
    #' @return nothing
    check = function(device = NULL) {
      if(is.null(device)) device = self$device_r

      self$device = torch::torch_device(device)
      self$dtype = torch::torch_float32()

      pars = list()
      self$parComp = check_and_recreate(self$parComp, self$parComp_r, dtype=torch::torch_float32(), device=self$device)
      self$parGrowth = check_and_recreate(self$parGrowth, self$parGrowth_r, dtype=torch::torch_float32(), device=self$device)
      self$parMort = check_and_recreate(self$parMort, self$parMort_r, dtype=torch::torch_float32(), device=self$device)
      self$parReg = check_and_recreate(self$parReg, self$parReg_r, dtype=torch::torch_float32(), device=self$device)

      # rebuild NN - if necessary? How to check?
      pointer_check <- tryCatch(torch::as_array(self$nnMortEnv$parameters[[1]]), error = function(e) e)
      if(inherits(pointer_check,"error")){

        self$nnMortEnv = do.call(self$build_NN, self$nnMortConfig)
        self$nnGrowthEnv = do.call(self$build_NN, self$nnGrowthConfig)
        self$nnRegEnv = do.call(self$build_NN, self$nnRegConfig)

        # and set weights
        .null = sapply(1:length(self$nnMortEnv$parameters), function(i) self$nnMortEnv$parameters[[i]]$set_data( torch::torch_tensor(self$parMortEnv_r[[i]], dtype=torch::torch_float32(), device=self$device)))
        .null = sapply(1:length(self$nnGrowthEnv$parameters), function(i) self$nnGrowthEnv$parameters[[i]]$set_data( torch::torch_tensor(self$parGrowthEnv_r[[i]], dtype=torch::torch_float32(), device=self$device)))
        .null = sapply(1:length(self$nnRegEnv$parameters), function(i) self$nnRegEnv$parameters[[i]]$set_data( torch::torch_tensor(self$parRegEnv_r[[i]], dtype=torch::torch_float32(), device=self$device)))

        self$parGrowthEnv = self$nnGrowthEnv$parameters
        .null = lapply(self$parGrowthEnv, function(p) p$requires_grad_(attr(self$parGrowthEnv_r[[1]], "requires_grad")) )
        self$parMortEnv = self$nnMortEnv$parameters
        .null = lapply(self$parMortEnv, function(p) p$requires_grad_(attr(self$parMortEnv_r[[1]], "requires_grad")) )
        self$parRegEnv = self$nnRegEnv$parameters
        .null = lapply(self$parRegEnv, function(p) p$requires_grad_(attr(self$parRegEnv_r[[1]], "requires_grad")))


        if(!attr(self$parMortEnv_r[[1]], "requires_grad")) .n = lapply(self$nnMortEnv$parameters, function(p) p$requires_grad_(FALSE))
        if(!attr(self$parGrowthEnv_r[[1]], "requires_grad")) .n = lapply(self$nnGrowthEnv$parameters, function(p) p$requires_grad_(FALSE))
        if(!attr(self$parRegEnv_r[[1]], "requires_grad")) .n = lapply(self$nnRegEnv$parameters, function(p) p$requires_grad_(FALSE))

        self$update_parameters()

      }

    },

    #' @description
    #'
    #' Update R fields of the parameters
    update_parameters = function() {
      pars = list()
      if(self$parComp$requires_grad) pars = c(pars, parComp = self$parComp)
      if(self$parGrowth$requires_grad) pars = c(pars, parGrowth = self$parGrowth)
      if(self$parReg$requires_grad) pars = c(pars, parReg = self$parReg)
      if(self$parMort$requires_grad) pars = c(pars, parMort = self$parMort)
      if(self$nnRegEnv$parameters[[1]]$requires_grad) pars = c(pars, nnReg = self$nnRegEnv$parameters)
      if(self$nnGrowthEnv$parameters[[1]]$requires_grad) pars = c(pars, nnGrowth = self$nnGrowthEnv$parameters)
      if(self$nnMortEnv$parameters[[1]]$requires_grad) pars = c(pars, nnMort = self$nnMortEnv$parameters)

      if(self$scale_2$requires_grad) pars = c(pars, scale_2 = self$scale_2)
      if(self$scale_4$requires_grad) pars = c(pars, scale_4 = self$scale_4)
      if(self$scale_5$requires_grad) pars = c(pars, scale_5 = self$scale_5)
      if(self$scale_6$requires_grad) pars = c(pars, scale_6 = self$scale_6)
      self$parameters = pars
    },

    #' @description
    #'
    #' Set the raw parameters within their ranges
    #'
    #' @param inputPar which parameter should be set
    #' @param parRange matrix of parameter ranges
    setPars = function(inputPar, parRange){
      # inputPar = torch::torch_tensor(inputPar, requires_grad = TRUE)
      if(is.vector(inputPar)) {
        Npar = 1
        Nsp = length(inputPar)
        inputPar = matrix(inputPar, nrow = Nsp, ncol = Npar)
        parRange = matrix(parRange, nrow = Npar, ncol = 2)
      }else if(is.matrix(inputPar)){
        Npar = ncol(inputPar)
        Nsp = nrow(inputPar)
      }else{
        stop("speciesPars and speciesPars_ranges must contain vectors or a matrices")
      }
      out = matrix(nrow = Nsp, ncol = 0)
      # for(j in 1:Nsp){
      for(i in 1:Npar){
        lower = as.numeric(parRange[i,1,drop=FALSE])
        upper = as.numeric(parRange[i,2,drop=FALSE])
        out <- cbind(out, qlogis((inputPar[,i, drop=FALSE] - lower) / (upper - lower)))
      }
      if(Npar == 1) out = as.vector(out) # TODO replace vectors with matrices as species input everywhere
      # }
      # out <- torch::torch_tensor(out, requires_grad = TRUE, device = "cpu", dtype = "float32")
      out <- torch::torch_tensor(out, requires_grad = TRUE, device = self$device, dtype = self$dtype)
      return(out)
    },

    #' @description
    #'
    #' Get transformed parameters (on the response scale)
    #'
    #' @param internalPar which parameter should be returned (raw internal version)
    #' @param parRange matrix of parameter ranges
    getPars = function(internalPar, parRange) {
      if (is.vector(parRange)) {
        # Case where internalPar is a 1D tensor and parRange is a vector
        Npar <- length(parRange) / 2
        lower <- parRange[1:Npar]
        upper <- parRange[(Npar + 1):(2 * Npar)]

        out <- torch::torch_sigmoid(internalPar) * (upper - lower) + lower
      } else {
        # Case where internalPar is a matrix and parRange is a matrix
        Npar <- ncol(internalPar)
        out <- list()
        for (i in 1:Npar) {
          lower <- parRange[i, 1, drop = FALSE]
          upper <- parRange[i, 2, drop = FALSE]
          out[[i]] <- torch::torch_sigmoid(internalPar[, i, drop = FALSE]) * (upper - lower) + lower
        }
        out <- torch::torch_cat(out, dim = 2L)
      }
      return(out)
    },

    #' @description
    #'
    #' Set function environment to object environment (so that `self$...` can be used within custom process functions)
    #' @param fn function
    set_environment = function(fn) {
      environment(fn) = self$.__enclos_env__
      return(fn)
    }

  ),
  active = list(

    #' @field parGrowthT return transformed parGrowth
    parGrowthT = function() {
      return(self$getPars(self$parGrowth, self$speciesPars_ranges$parGrowth))
    },
    #' @field parMortT return transformed parMort
    parMortT = function() {
      return(self$getPars(self$parMort, self$speciesPars_ranges$parMort))
    },
    #' @field parRegT return transformed parReg
    parRegT = function() {
      return(self$getPars(self$parReg, self$speciesPars_ranges$parReg))
    },
    #' @field parCompT return transformed parComp
    parCompT = function() {
      return(self$getPars(self$parComp, self$speciesPars_ranges$parComp))
    },

    #' @field parGrowthTR return transformed parGrowth
    parGrowthTR = function() {
      return(self$getPars(self$parGrowth, self$speciesPars_ranges$parGrowth) |> as.matrix())
    },
    #' @field parMortTR return transformed parMort
    parMortTR = function() {
      return(self$getPars(self$parMort, self$speciesPars_ranges$parMort)  |> as.matrix())
    },
    #' @field parRegTR return transformed parReg
    parRegTR = function() {
      return(self$getPars(self$parReg, self$speciesPars_ranges$parReg)  |> as.numeric())
    },
    #' @field parCompTR return transformed parComp
    parCompTR = function() {
      return(self$getPars(self$parComp, self$speciesPars_ranges$parComp)  |> as.numeric())
    }
  )
)

