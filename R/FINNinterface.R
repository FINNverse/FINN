#' Create a Process Object
#'
#' This function creates a process object that is used to define a specific ecological process (e.g., growth, mortality, regeneration) within a forest simulation model. The process object can include custom functions, allometric parameters, environmental initialization, and hidden layers for neural networks.
#'
#' @param formula An optional formula specifying the structure of the model. Default is `NULL`, which results in the formula `~.` being used.
#' @param func A custom function to define the process. This is a required parameter.
#' @param optimize Logical indicating whether the process function should be optimized. Default is `FALSE`.
#' @param initAllometric Initial allometric parameters for the process. Default is `NULL`.
#' @param initEnv Initial environmental parameters for the process. Default is `NULL`.
#' @param hidden A list specifying the hidden layers for neural network models. Default is an empty list.
#'
#' @return A list of class "process" containing the process definition and associated parameters.
#'
#' @examples
#' growth_process <- createProcess(formula = ~temperature + precipitation, func = growthFunction)
#'
#' @export
createProcess = function(formula = NULL, func, optimize = FALSE, initAllometric = NULL, initEnv = NULL, hidden = list()) {
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
  out$initAllometric = initAllometric
  out$initEnv = initEnv
  out$hidden = NULL

  if(out$custom) cli::cli_text("Custom function detected...")
  if(!is.null(out$hidden)) cli::cli_text("Neural Network detected...")
  class(out) = "process"
  return(out)
}


#' Extract Environmental Data for a Process
#'
#' This function extracts and transforms environmental data according to a specified process object. The environmental data is formatted into an array suitable for input into the simulation model.
#'
#' @param process A process object created using `createProcess`.
#' @param env A data frame containing environmental data, which must include columns "siteID" and "year".
#'
#' @return An array of environmental data formatted according to the process object.
#'
#' @import data.table
#' @examples
#' env_array <- extract_env(growth_process, env_data)
extract_env = function(process, env) {
  require(data.table)
  process_env = stats::model.matrix(process$formula, env[,-c("siteID", "year")])
  env_names = colnames(process_env)
  env_array = climateDF2array(climate_dt =  cbind(env[,c("siteID", "year")],process_env), env_vars = env_names)
  return(env_array)
}

#' Simulate Forest Dynamics
#'
#' This function simulates forest dynamics, including mortality, growth, regeneration, and competition,
#' across multiple patches and species.
#'
#' @param env A data frame containing the environmental variables used in the simulation.
#' @param mortality A process object or formula defining the mortality process.
#' @param growth A process object or formula defining the growth process.
#' @param regeneration A process object or formula defining the regeneration process.
#' @param competition A function defining the competition process, if applicable.
#' @param height A function or set of parameters defining the height-growth relationship, if applicable.
#' @param patches An integer specifying the number of patches in the simulation. Default is 10.
#' @param patch_size A numeric value representing the size of each patch. Default is 0.1.
#' @param sp An integer specifying the number of species in the simulation. Default is 5.
#' @param init A custom initialization object for the simulation, if available.
#' @param device A character string specifying whether to use `"cpu"` or `"gpu"` for computation. Default is `"cpu"`.
#' @param parallel A logical value indicating whether to run the simulation in parallel. Default is `FALSE`.
#' @param NGPU An integer specifying the number of GPUs to use if `device = "gpu"`. Default is 1.
#'
#' @return A list of predicted values for each patch and species, containing arrays of tree diameters and tree counts.
#'
#' @examples
#' env_data <- data.frame(siteID = 1:5, year = 2001:2005, temp = rnorm(5), precip = rnorm(5))
#' predictions <- simulateForest(env = env_data, patches = 3, sp = 2)
#'
#' @export
simulateForest = function(env,
                          mortalityProcess = NULL,
                          growthProcess = NULL,
                          regenerationProcess = NULL,
                          competitionProcess = NULL,
                          height = NULL,
                          patches = 10L,
                          patch_size = 0.1,
                          sp = 5L,
                          init = NULL,
                          device = c("cpu", "gpu"),
                          parallel = FALSE,
                          NGPU = 1,
                          batchsize = NULL
) {

  out = list()
  device = match.arg(device)
  if(is.null(mortalityProcess)) mortalityProcess = createProcess(formula = ~., func = mortality)
  if(is.null(growthProcess)) growthProcess = createProcess(formula = ~., func = growth)
  if(is.null(regenerationProcess)) regenerationProcess = createProcess(formula = ~., func = regeneration)
  if(is.null(competitionProcess)) competitionProcess = createProcess(formula = ~., func = competition)

  if(inherits(mortalityProcess, "formula")) mortalityProcess = createProcess(formula = mortalityProcess, func = mortality)
  if(inherits(growthProcess, "formula")) growthProcess = createProcess(formula = growthProcess, func = growth)
  if(inherits(regenerationProcess, "formula")) regenerationProcess = createProcess(formula = regenerationProcess, func = regeneration)

  ### Year only ###
  mortality_env = extract_env(mortalityProcess, env)
  growth_env = extract_env(growthProcess, env)
  regeneration_env = extract_env(regenerationProcess, env)


  sites = dim(mortality_env)[1] # TODO check that all env arrays have the same sites and timesteps
  if(is.null(batchsize)) batchsize = sites

  if(is.null(init)) {
    init = CohortMat$new(dims = c(sites, patches, sp),
                         dbh = array(1, dim = c(sites, patches, sp)),
                         trees = array(1, dim = c(sites, patches, sp)),
                         sp = sp)
  }
  device_old = device
  if(device == "gpu") {
    device = "cuda:0"
    if(!torch::cuda_is_available()) {
      cli::cli_text("GPU device not available...")
      device = "cpu"
    }
  }

  if(is.logical(parallel)) {
    if(parallel) NCPUs = parallel::detectCores()-1
    else NCPUs = -999
  }

  if(is.numeric(parallel)) NCPUs = parallel


  model = FINN$new(sp = sp,
                   env = c(dim(mortality_env)[3], dim(growth_env)[3], dim(regeneration_env)[3]),
                   device = device,
                   which = "all" ,
                   hidden_growth = growthProcess$hidden,
                   hidden_mort = mortalityProcess$hidden,
                   hidden_reg = regenerationProcess$hidden,
                   growthFunction = growthProcess$func,
                   mortalityFunction = mortalityProcess$func,
                   regenerationFunction = regenerationProcess$func,
                   competitionFunction = competitionProcess$func,
                   parGrowth = growthProcess$initAllometric,
                   parMort = mortalityProcess$initAllometric,
                   parReg = regenerationProcess$initAllometric,
                   parHeight = height,
                   parGrowthEnv = growthProcess$initEnv,
                   parMortEnv = mortalityProcess$initEnv,
                   parRegEnv = regenerationProcess$initEnv,
                   patch_size_ha = patch_size)

  if(NCPUs < -0.5) {

    predictions = model$predict(dbh = init$dbh,
                                trees = init$trees,
                                species = init$species,
                                env = list(torch::torch_tensor(mortality_env),
                                           torch::torch_tensor(growth_env),
                                           torch::torch_tensor(regeneration_env)),
                                patches = patches,
                                debug = FALSE,
                                verbose = TRUE)

  } else {
    cl = parallel::makeCluster(NCPUs)
    nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))

    batches = cbind(c(1, (seq(1, sites, by = batchsize))[-1]), c((seq(1, sites, by = batchsize)-1)[-1], sites)  )

    parallel::clusterEvalQ(cl, {library(torch);library(FINN)})
    parallel::clusterExport(cl, varlist = c("FINN","CohortMat","model", "init", "mortality_env", "growth_env", "regeneration_env", "batches", "patches", "nodes", "NGPU", "device_old"),envir = environment())
    # parallel::clusterEvalQ(cl, {
    #
    #   # who am I
    #   if(device_old == "gpu") {
    #     myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
    #     dist = cbind(nodes,0:(NGPU-1))
    #     dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
    #     device_hardware = paste0("cuda:",dev)
    #   } else {
    #     device_hardware = "cpu"
    #   }
    #
    #    init$check()
    #    model$check(device = device_hardware)
    #
    #   print(1)
    # })
    predictions =
      parallel::parLapply(cl, 1:nrow(batches), function(i) {

        if(device_old == "gpu") {
          myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
          dist = cbind(nodes,0:(NGPU-1))
          dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
          device_hardware = paste0("cuda:",dev)
        } else {
          device_hardware = "cpu"
        }

        init$check()
        model$check(device = device_hardware)
        pred = model$predict(dbh = init$dbh[batches[i,1]:batches[i,2],,,drop=FALSE],
                            trees = init$trees[batches[i,1]:batches[i,2],,,drop=FALSE],
                            species = init$species[batches[i,1]:batches[i,2],,,drop=FALSE],
                            env = list(torch::torch_tensor(mortality_env[batches[i,1]:batches[i,2],,,drop=FALSE]),
                                       torch::torch_tensor(growth_env[batches[i,1]:batches[i,2],,,drop=FALSE]),
                                       torch::torch_tensor(regeneration_env[batches[i,1]:batches[i,2],,,drop=FALSE])),
                            patches = patches,
                            debug = FALSE,
                            verbose = TRUE)
        return(pred)
      })
    parallel::stopCluster(cl)
  }
  return(predictions)

}





#' Forest informed neural network
#'
#' This function simulates forest dynamics, including mortality, growth, regeneration, and competition,
#' across multiple patches and species.
#'
#' @param data data containing species information (dbh, ba, trees, AL, growth, mortality rates, and regeneration rates)
#' @param env A data frame containing the environmental variables used in the simulation.
#' @param mortality A process object or formula defining the mortality process.
#' @param growth A process object or formula defining the growth process.
#' @param regeneration A process object or formula defining the regeneration process.
#' @param competition A function defining the competition process, if applicable.
#' @param height A function or set of parameters defining the height-growth relationship, if applicable.
#' @param patches An integer specifying the number of patches in the simulation. Default is 10.
#' @param patch_size A numeric value representing the size of each patch. Default is 0.1.
#' @param sp An integer specifying the number of species in the simulation. Default is 5.
#' @param init A custom initialization object for the simulation, if available.
#' @param device A character string specifying whether to use `"cpu"` or `"gpu"` for computation. Default is `"cpu"`.
#' @param bootstrap bootstrap model or not (computationally expensive!)
#' @param parallel A logical value indicating whether to run the bootstrapping in parallel. Default is `FALSE`.
#' @param NGPU An integer specifying the number of GPUs to use if `device = "gpu"`. Default is 1 (only used for parallel bootstrapping).
#' @param ... arguments passed to `simulateForest()`
#'
#' @return A list of predicted values for each patch and species, containing arrays of tree diameters and tree counts.
#'
#' @examples
#' env_data <- data.frame(siteID = 1:5, year = 2001:2005, temp = rnorm(5), precip = rnorm(5))
#' predictions <- simulateForest(env = env_data, patches = 3, sp = 2)
#'
#' @export
finn = function(data = NULL,
                env,
                mortalityProcess = NULL,
                growthProcess = NULL,
                regenerationProcess = NULL,
                competitionProcess = NULL,
                height = NULL,
                patches = 10L,
                patch_size = 0.1,
                sp = 5L,
                init = NULL,
                device = c("cpu", "gpu"),
                bootstrap = NULL,
                parallel = FALSE,
                NGPU = 1,
                batchsize = NULL,
                ...
) {

  out = list()
  device = match.arg(device)
  if(is.null(mortalityProcess)) mortalityProcess = createProcess(formula = ~., func = mortality)
  if(is.null(growthProcess)) growthProcess = createProcess(formula = ~., func = growth)
  if(is.null(regenerationProcess)) regenerationProcess = createProcess(formula = ~., func = regeneration)
  if(is.null(competitionProcess)) competitionProcess = createProcess(formula = ~., func = competition)

  if(inherits(mortalityProcess, "formula")) mortalityProcess = createProcess(formula = mortalityProcess, func = mortality)
  if(inherits(growthProcess, "formula")) growthProcess = createProcess(formula = growthProcess, func = growth)
  if(inherits(regenerationProcess, "formula")) regenerationProcess = createProcess(formula = regenerationProcess, func = regeneration)

  ### Year only ###
  mortality_env = extract_env(mortalityProcess, env)
  growth_env = extract_env(growthProcess, env)
  regeneration_env = extract_env(regenerationProcess, env)


  sites = dim(mortality_env)[1] # TODO check that all env arrays have the same sites and timesteps
  if(is.null(batchsize)) batchsize = sites

  if(is.null(init)) {
    init = CohortMat$new(dims = c(sites, patches, sp),
                         dbh = array(1, dim = c(sites, patches, sp)),
                         trees = array(1, dim = c(sites, patches, sp)),
                         sp = sp)
  }
  device_old = device
  if(device == "gpu") {
    device = "cuda:0"
    if(!torch::cuda_is_available()) {
      cli::cli_text("GPU device not available...")
      device = "cpu"
    }
  }

  if(is.logical(parallel)) {
    if(parallel) NCPUs = parallel::detectCores()-1
    else NCPUs = -999
  }

  if(is.numeric(parallel)) NCPUs = parallel


  model = FINN$new(sp = sp,
                   env = c(dim(mortality_env)[3], dim(growth_env)[3], dim(regeneration_env)[3]),
                   device = device,
                   which = "all" ,
                   hidden_growth = growthProcess$hidden,
                   hidden_mort = mortalityProcess$hidden,
                   hidden_reg = regenerationProcess$hidden,
                   growthFunction = growthProcess$func,
                   mortalityFunction = mortalityProcess$func,
                   regenerationFunction = regenerationProcess$func,
                   competitionFunction = competitionProcess$func,
                   parGrowth = growthProcess$initAllometric,
                   parMort = mortalityProcess$initAllometric,
                   parReg = regenerationProcess$initAllometric,
                   parHeight = height,
                   parGrowthEnv = growthProcess$initEnv,
                   parMortEnv = mortalityProcess$initEnv,
                   parRegEnv = regenerationProcess$initEnv,
                   patch_size_ha = patch_size)

  if(NCPUs < -0.5) {

    predictions = model$predict(dbh = init$dbh,
                                trees = init$trees,
                                species = init$species,
                                env = list(torch::torch_tensor(mortality_env),
                                           torch::torch_tensor(growth_env),
                                           torch::torch_tensor(regeneration_env)),
                                patches = patches,
                                debug = FALSE,
                                verbose = TRUE)

  } else {
    cl = parallel::makeCluster(NCPUs)
    nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))

    batches = cbind(c(1, (seq(1, sites, by = batchsize))[-1]), c((seq(1, sites, by = batchsize)-1)[-1], sites)  )

    parallel::clusterEvalQ(cl, {library(torch);library(FINN)})
    parallel::clusterExport(cl, varlist = c("FINN","CohortMat","model", "init", "mortality_env", "growth_env", "regeneration_env", "batches", "patches", "nodes", "NGPU", "device_old"),envir = environment())
    # parallel::clusterEvalQ(cl, {
    #
    #   # who am I
    #   if(device_old == "gpu") {
    #     myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
    #     dist = cbind(nodes,0:(NGPU-1))
    #     dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
    #     device_hardware = paste0("cuda:",dev)
    #   } else {
    #     device_hardware = "cpu"
    #   }
    #
    #    init$check()
    #    model$check(device = device_hardware)
    #
    #   print(1)
    # })
    predictions =
      parallel::parLapply(cl, 1:nrow(batches), function(i) {

        if(device_old == "gpu") {
          myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
          dist = cbind(nodes,0:(NGPU-1))
          dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
          device_hardware = paste0("cuda:",dev)
        } else {
          device_hardware = "cpu"
        }

        init$check()
        model$check(device = device_hardware)
        pred = model$predict(dbh = init$dbh[batches[i,1]:batches[i,2],,,drop=FALSE],
                             trees = init$trees[batches[i,1]:batches[i,2],,,drop=FALSE],
                             species = init$species[batches[i,1]:batches[i,2],,,drop=FALSE],
                             env = list(torch::torch_tensor(mortality_env[batches[i,1]:batches[i,2],,,drop=FALSE]),
                                        torch::torch_tensor(growth_env[batches[i,1]:batches[i,2],,,drop=FALSE]),
                                        torch::torch_tensor(regeneration_env[batches[i,1]:batches[i,2],,,drop=FALSE])),
                             patches = patches,
                             debug = FALSE,
                             verbose = TRUE)
        return(pred)
      })
    parallel::stopCluster(cl)
  }
  return(predictions)

}

