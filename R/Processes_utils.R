#' Create a Process Object
#'
#' This function creates a process object that is used to define a specific ecological process (e.g., growth, mortality, regeneration) within a forest simulation model. The process object can include custom functions, species parameters, environmental initialization, and hidden layers for neural networks.
#'
#' @param formula An optional formula specifying the structure of the model. Default is `NULL`, which results in the formula `~.` being used.
#' @param func A custom function to define the process. This is a required parameter.
#' @param optimize Logical indicating whether the process function should be optimized. Default is `FALSE`.
#' @param initSpecies Initial species parameters for the process. Default is `NULL`.
#' @param initEnv Initial environmental parameters for the process. Default is `NULL`.
#' @param hidden A list specifying the hidden layers for neural network models. Default is an empty list.
#' @param inputNN input dimension for NN, default is inferred from the formula. See details
#' @param outputNN output dimension for NN, default is the number of species. See details
#' @param dispersion_parameter init dispersion parameter, if available (currently only supported for regeneration rate that is based on a negative binomial).
#' @param NN pass custom NN to model
#' @param upper upper boundaries of species parameters
#' @param lower lower boundaries of species parameters
#' @param dropout dropout rate of neural networks
#' @param sample_regeneration sample recruits or not
#'
#' @return A list of class "process" containing the process definition and associated parameters.
#'
#' @examples
#' growth_process <- createProcess(formula = ~temperature + precipitation, func = growthFunction)
#'
#' @export
createProcess = function(formula = NULL, func, initSpecies = NULL, initEnv = NULL, hidden = NULL, optimizeSpecies = FALSE, optimizeEnv = TRUE, inputNN = NULL, outputNN = NULL, dispersion_parameter = 1.0, NN = NULL, upper = NULL, lower = NULL, dropout = 0.0, sample_regeneration = TRUE) {
  out = list()
  if(!is.null(formula)){
    mf = match.call()
    m = match("formula", names(mf))
    if(inherits(mf[3]$formula, "name")) mf[3]$formula = eval(mf[3]$formula, envir = parent.env(environment()))
    formula = stats::as.formula(mf[m]$formula)
    #X = stats::model.matrix(formula, data)
  } else {
    formula = stats::as.formula("~.")
    #X = stats::model.matrix(formula, data)
  }

  out$formula = formula
  out$func = func
  out$custom = TRUE
  if(isNamespace(environment(func))) out$custom = getNamespaceName(environment(func)) != "FINN"
  out$initSpecies = initSpecies
  out$optimizeSpecies = optimizeSpecies
  out$optimizeEnv = optimizeEnv
  out$inputNN = inputNN
  out$outputNN = outputNN
  out$dispersion_parameter = dispersion_parameter
  out$sample_regeneration = TRUE
  out$NN = NN
  out$upper = NULL
  out$lower = NULL
  out$dropout = dropout
  if(!is.null(initEnv)) {
    if(!is.list(initEnv)) initEnv = list(initEnv) # must be a list!
  }
  out$initEnv = initEnv
  out$hidden = hidden

  if(out$custom) cli::cli_text("Custom function detected...")
  if(!is.null(out$hidden)) cli::cli_text("Neural Network detected...")
  class(out) = "process"
  return(out)
}


#' Create a hybrid model Object
#'
#' This function creates a process object that is used to define a specific ecological process (e.g., growth, mortality, regeneration) within a forest simulation model. The process object can include custom functions, species parameters, environmental initialization, and hidden layers for neural networks.
#'
#' @param formula An optional formula specifying the structure of the model. Default is `NULL`, which results in the formula `~.` being used.
#' @param optimize Logical indicating whether the process function should be optimized. Default is `FALSE`.
#' @param hidden A list specifying the hidden layers for neural network models. Default is an empty list.
#' @param dispersion_parameter init dispersion parameter, if available (currently only supported for regeneration rate that is based on a negative binomial).
#' @param NN pass custom NN to model
#' @param dropout dropout rate of neural networks
#'
#' @return A list of class "process" containing the process definition and associated parameters.
#'
#' @examples
#' growth_process <- createProcess(formula = ~temperature + precipitation, func = growthFunction)
#'
#' @export
createHybrid = function(formula = NULL, hidden = NULL, optimize = TRUE, dispersion_parameter = 1.0, NN = NULL, dropout = 0.0, encoder_layers = 1L) {
  out = list()
  if(!is.null(formula)){
    mf = match.call()
    m = match("formula", names(mf))
    if(inherits(mf[3]$formula, "name")) mf[3]$formula = eval(mf[3]$formula, envir = parent.env(environment()))
    formula = stats::as.formula(mf[m]$formula)
    #X = stats::model.matrix(formula, data)
  } else {
    formula = stats::as.formula("~.")
    #X = stats::model.matrix(formula, data)
  }

  out$formula = formula
  out$custom = TRUE
  out$optimizeSpecies = FALSE
  out$optimizeEnv = optimize
  out$dispersion_parameter = dispersion_parameter
  out$NN = NN
  out$dropout = dropout
  out$encoder_layers = encoder_layers
  out$hidden = hidden
  class(out) = "hybrid"
  return(out)
}

#' Extract Environmental Data for a Process
#'
#' This function extracts and transforms environmental data according to a specified process object. The environmental data is formatted into an array suitable for input into the simulation model.
#'
#' @param formula formula object `createProcess`.
#' @param env A data frame containing environmental data, which must include columns "siteID" and "year".
#'
#' @return An array of environmental data formatted according to the process object.
#'
#' @import data.table
#' @examples
#' env_array <- extract_env(growth_process, env_data)
extract_env = function(formula, env) {
  require(data.table)
  process_env = stats::model.matrix(formula, env[,-c("siteID", "year")])
  env_names = colnames(process_env)
  env_array = climateDF2array(climate_dt =  cbind(env[,c("siteID", "year")],process_env), env_vars = env_names)
  return(env_array)
}

backward = function(value, upper, lower) stats::qlogis( t((t(value) - as.numeric(lower)) / (as.numeric(upper) - as.numeric(lower)) ))
forward =  function(value, upper, lower) torch::torch_sigmoid(value) * (upper - lower) + lower

init_species_parameters = function(type, N_species){
  initLight = stats::runif(N_species, 0.45, 0.55)
  if(type == "competition") {
    init =
      cbind(initLight,
            stats::runif(N_species, min = 0.19, 0.21)) # 0.2 corresponds to 0 light at 50m2/ha basal area)
  }
  if(type == "growth") {
    init =
      cbind(initLight,
            stats::runif(N_species, min = 0.045, 0.055))
  }
  if(type == "mortality") {
    init =
      cbind(stats::runif(N_species, min = -0.2, 0.2),
            stats::runif(N_species, min = -0.2, 0.2),
            stats::runif(N_species, min = -0.2, 0.2))
  }
  if(type == "regeneration") {
    init = matrix(initLight, ncol = 1L)
  }
  return(init)
}


get_par_boundary = function(obj, type, upper = TRUE) {
  if(upper) {
    if(is.null(obj$upper)) {
      upper = switch(type,
                     #mortality = { c(0.93, 4.00) },
                     mortality = { c(10.0, 10.0, 10.0) },
                     growth = {c(0.99, 4.0)},
                     regeneration = { 0.99 },
                     competition = { c(0.7, 2.0)})
    } else {
      upper = obj$upper
    }
    return(upper)
  } else {
    if(is.null(obj$lower)) {
      lower = switch(type,
                     mortality = { c(-10.0, -10.0, -10.0) },
                     #mortality = { c(0.011, 0.00) },
                     growth =    {c(0.01, 0.01)},
                     regeneration = { 0.01 },
                     competition = { c(0.3, 0.0)})
    } else {
      lower = obj$lower
    }
    return(lower)
  }
}
