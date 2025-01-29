#' Forest informed neural network
#'
#' @description
#' \loadmathjax
#' This function is the main interface to the FINN platform. The function can be used to simulate and to estimate (fit) the parameters of the FINN model.
#'data
#'
#' @param data (`data.table`) \cr
#'   Data containing species information, must include the columns siteID, species, year, dbh, ba, trees, AL, growth, mort, and reg.
#' @param env (`data.table`) \cr
#'   Data containing environmental predictors, must include the columns siteID, year, and names of environmental predictors used in the formula argument of the processes.
#' @param disturbance (`data.table`) \cr
#'   Data containing disturbance events, must include the columns siteID, year, and intensity.
#' @param mortalityProcess (`NULL`, `formula`, or `process` created by `createProcess()`) \cr
#'  A process object or formula defining the mortality process. If `NULL`, default `createProcess(...)` object will be used (all envrionmental predictors will be used).
#' @param growthProcess (`NULL`, `formula`, or `process` created by `createProcess()`) \cr
#'  A process object or formula defining the growth process. If `NULL`, default `createProcess(...)` object will be used (all envrionmental predictors will be used).
#' @param regenerationProcess (`NULL`, `formula`, or `process` created by `createProcess()`) \cr
#'  A process object or formula defining the regeneration process. If `NULL`, default `createProcess(...)` object will be used (all envrionmental predictors will be used).
#' @param competitionProcess (`NULL`, `formula`, or `process` created by `createProcess()`) \cr
#'  A function defining the competition process, if applicable. If `NULL`, default `createProcess(...)` object will be used.
#' @param height (`numeric(sp)`) \cr
#'  Set of parameters defining the height-growth relationship, if applicable.
#' @param optimizeHeight (`logical(1)`) \cr
#'  Should the height-growth parameters be estimated or not.
#' @param speciesPars_ranges (`list(parGrowth = list(matrix()), parMort = list(matrix()), parReg = list(matrix()), parComp = list(matrix()))`) \cr
#'  List of boundaries for species parameters.
#' @param patches (`integer(1)`) \cr
#'  An integer specifying the number of patches in the simulation. Default is 10.
#' @param patch_size (`numeric(1)`) \cr
#'  A numeric value representing the size of each patch. Default is 0.1.
#' @param init (`NULL` or `CohortMat` created by `CohortMat$new()`) \cr
#'  A custom initialization object for the simulation, if available.
#' @param batchsize (`NULL` or `integer(1)`) \cr
#'  batchsize (sites can be feed in chunkes to the predict function, can be necessary if GPU memory is limited).
#' @param epochs (`integer(1)`) \cr
#'  Number of iteration (optimization) steps.
#' @param lr (`numeric(1)`) \cr
#'  Learning rate for the first-gradient descent optimizer.
#' @param device (`character(1)`) \cr
#'  A character string specifying whether to use `"cpu"` or `"gpu"` for computation. Default is `"cpu"`.
#' @param bootstrap (`logical(1)` or `integer(1)`) \cr
#'  bootstrap model or not (computationally expensive!).
#' @param parallel (`character(1)`) \cr
#'  A logical value indicating whether to run the bootstrapping in parallel. Default is `FALSE`.
#' @param NGPU (`character(1)`) \cr
#'  An integer specifying the number of GPUs to use if `device = "gpu"`. Default is 1 (only used for parallel bootstrapping).
#' @param weights (`numeric(6`) \cr
#'  Weighting of the 6 errors (ba, trees, AL, growth, mortality, and regeneration).
#' @param file (`character()`) \cr
#'  If weights should be saved after each optimization step (for monitoring), set a path. Default is `NULL`.
#' @param thin (`integer(1)`) \cr
#'  If weights should be save (file argument), save every thin epochs.
#' @param optimizer (`function()`) \cr
#'  Optimizer of the torch package
#' @param ... arguments passed to `simulateForest()`
#'
#'
#' @details
#'
#' While FINN is a platform for modular dynamic forest models, we provide default functions for the four processes, namely `?mortality`, `?growth`, `?regeneration`, and `?competetition`.
#'
#' The idea of finn is to combine mechanistic understanding (the processes) and data-driven approaches. Data-driven means that parameters of the model can be estimated (or informed) by empirical data. This is optional and can be realized in different ways or levels in finn. By default, we assume that the three processes mortality, growth, and regeneration depend on the environment (niche of species), available light (as a consequence of competition), and species specific parameters. The environmental component and the species specific parameters can be estimated or configured. The central function to configure each process is the `?createProcess` function.
#'
#' In the `createProcess()` function, the parameters of the environmental model (by default a linear regression model) can be set by the `initEnv = ...` argument (must be a matrix with n rows for n species and m columns for m environmental predictors (which is set by the formula argument)). Species specific parameters are set by the `initSpecies=...` argument (check help of the default processes for dimensions). If `NULL` the parameters will be sampled from plausible ranges that serve as initial values for the optimization, if inference is desired (which is set by `optimizeSpecies` and `optimizeEnv` arguments). Moreover, custom process functions can be passed to the `createProcess()` function using the `func` argument. However, take caution that all operations in the custom process function are compatible with `torch`, work with multi-dimensional arrays, and are differentiable when inference is wanted. See the section below for more information about custom process functions.
#'
#' ## Custom Process Functions:
#'
#' FINN is based on torch for R which is an library for Deep Learning. The advantage of torch is that all calculations are highly optimized and can be run on the GPU. To fully exploit the potential of GPU, all operations in FINN are based on multi-dimensional arrays, tensors. Usually, the dimensions of data objects within finn are of \mjseqn{[k_sites, m_{patches}, n_{cohorts}]}. For custom process functions to properly work, they must work on these tensors and handle the dimensions correctly (in the future we will provide helper functions to check automatically custom functions):
#'
#' Part of the hybrid modeling example:
#'
#' ```{r, eval=FALSE}
#' growthCustom = function(dbh, species, parGrowth, pred, light, debug = F){
#' # dbh, species, light dimensions -> [sites, patches, cohorts]
#' # parGrowth species parameters, ignore
#' # pred -> [sites, number of environmental predictors]
#' input =
#'     torch::torch_cat(
#'       list(dbh$unsqueeze(4L)/200., # rescale dbh
#'            pred$unsqueeze(2L)$unsqueeze(2L)$`repeat`(c(1L, dbh$shape[2], dbh$shape[3], 1L)),
#'            light$unsqueeze(4L),
#'            torch::nnf_one_hot(species, num_classes = 3L)),
#'       dim = 4L)
#'   growth = self$nnGrowthEnv(input)[,,,1]$exp()
#'   return(growth)
#' }
#' ```
#' Within the custom growth function, the `torch::torch_cat()` function concatenates dbh, light, environment, and species information (one hot encoded species). The concatenation happens in an additional 4th dimension. `$unsqueeze(dim = ..)` adds a new dimension at the specific dimension. For pred, we must add two additional dimensions and then repeat the values along these new dimensions to get the same number of patches and cohorts (`dbh$shape[2], dbh$shape[3]`).
#'
#' Moreover, the custom functions must be differentiable, in short, avoid complex R functions or functions from other packages, do not transform the tensors into R objects (e.g. by using `torch::as_array` or `as.matrix`), and be careful with indexing (we have helper functions for that (TODO!)).
#'
#'
#'
#'
#' @return A list of predicted values for each patch and species, containing arrays of tree diameters and tree counts.
#'
#' @export
finn = function(data = NULL,
                env,
                disturbance = NULL,
                mortalityProcess = NULL,
                growthProcess = NULL,
                regenerationProcess = NULL,
                competitionProcess = NULL,
                height = NULL,
                # optimizeHeight = FALSE,
                speciesPars_ranges = list(
                  parGrowth = rbind(
                    c(0.01, 0.99),
                    c(0.01, 4)
                  ),
                  parMort = rbind(
                    c(0.011,0.93),
                    c(0, 4)
                  ),
                  parReg = c(0.01, 0.99),
                  # parHeight = c(0.3, 0.7),
                  parComp = rbind(
                    c(0.3, 0.7),
                    c(0, 2)
                  ),
                  parGrowthEnv = rbind(
                    c(-1, 1),
                    c(-1, 1)
                  ),
                  parMortEnv = rbind(
                    c(-2, 2),
                    c(-2, 2)
                  ),
                  parRegEnv = rbind(
                    c(-2, 2),
                    c(-2, 2)
                  )),
                patches = 100L,
                patch_size = 0.1,
                init = NULL,
                batchsize = NULL,
                epochs = 20L,
                lr = 0.1,
                device = c("cpu", "gpu"),
                bootstrap = NULL,
                parallel = FALSE,
                NGPU = 1,
                weights = c(0.05, 0.50, 3.00, 0.50, 3.00, 2.00),
                file = NULL,
                thin = 100,
                sp = NULL,
                optimizer = optim_ignite_adam,
                ...
) {

  limits = get_mort_thresholds(base_steepness = 5)
  if(speciesPars_ranges$parMort[1,1] < limits[1]) {
    cat('Mortality species parameter lower limit will be adjusted to improve numerical stability\n')
    speciesPars_ranges$parMort[1,1] = limits[1]
  }
  if(speciesPars_ranges$parMort[1,2] > limits[2]) {
    cat('Mortality species parameter upper limit will be adjusted to improve numerical stability\n')
    speciesPars_ranges$parMort[1,2] = limits[2]
  }

  if(is.null(data)) {
    cli::cli_text("No data provided, simulations will be generated...")
    return(simulateForest(env=env,
                          mortalityProcess = mortalityProcess,
                          growthProcess = growthProcess,
                          regenerationProcess = regenerationProcess,
                          competitionProcess = competitionProcess,
                          height = height,
                          patches = patches,
                          patch_size = patch_size,
                          init = init, device = device,
                          parallel = parallel,
                          NGPU = NGPU,
                          batchsize = batchsize,
                          ...))
  }

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


  if(is.null(sp)){
    message(paste0("No number of species ('sp') provided, assuming sp = ", length(unique(data$species)), " species from data..."))
    sp = length(unique(data$species))
  }
  # if(!growthProcess$optimizeSpecies) if(nrow(parGrowth) != sp) stop("Number of species in parGrowth does not match sp")
  # if(!mortalityProcess$optimizeSpecies) if(nrow(parMort) != sp) stop("Number of species in parMort does not match sp")
  # if(!regenerationProcess$optimizeSpecies) if(length(parReg) != sp) stop("Number of species in parReg does not match sp")
  # if(!growthProcess$optimizeEnv) if(nrow(parGrowthEnv) != sp) stop("Number of species in parGrowthEnv does not match sp")
  # if(!mortalityProcess$optimizeEnv) if(nrow(parMortEnv) != sp) stop("Number of species in parMortEnv does not match sp")
  # if(!regenerationProcess$optimizeEnv) if(nrow(parRegEnv) != sp) stop("Number of species in parRegEnv does not match sp")
  # if(!optimizeHeight) if(!is.null(height) && length(height) != sp) stop("Number of species in height does not match sp")
  out$sp = sp
  out$lr = lr


  #### Prepare response ####
  # TODO: improve!!!
  options(na.action='na.pass')
  response = list(
    dbh = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+dbh), data[data$species==i,])), along = 3L),
    ba = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+ba), data[data$species==i,])), along = 3L),
    trees = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+trees), data[data$species==i,])), along = 3L),
    AL = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+AL), data[data$species==i,])), along = 3L),
    growth = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+growth), data[data$species==i,])), along = 3L),
    mort = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+mort), data[data$species==i,])), along = 3L),
    reg = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+reg), data[data$species==i,])), along = 3L)
  )
  options(na.action='na.omit')

  year_sequence = which(levels(as.factor(env$year)) %in% levels(as.factor(data$year)), arr.ind = TRUE)

  disturbance_T = NULL
  if(!is.null(disturbance)) {
    disturbance = extract_env(list(formula=~0+intensity), disturbance)
    disturbance_T = torch::torch_tensor(disturbance, dtype=torch::torch_float32(), device="cpu")
  }

  response_T = torch::torch_cat(lapply(response, function(x) torch::torch_tensor(x, dtype=torch::torch_float32(), device="cpu")$unsqueeze(4)), 4)

  sites = dim(mortality_env)[1] # TODO check that all env arrays have the same sites and timesteps
  if(is.null(batchsize)) batchsize = sites

  if(batchsize > sites) batchsize = sites

  out$batchsize = batchsize

  if(is.null(init)) {
    init = CohortMat$new(dims = c(sites, patches, sp),
                         dbh = array(1, dim = c(sites, patches, sp)),
                         trees = array(1, dim = c(sites, patches, sp)),
                         sp = sp)
  } else {
    if(init$sp != sp) stop(paste("sp in cohort", init$sp, "does not match sp from data", sp))
    patches = dim(init$species_r)[2]
    init$check()
  }

  out$patches = patches
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


  runEnvMort = TRUE
  runEnvGrowth = TRUE
  runEnvReg = TRUE

  input_dimensions = c(dim(mortality_env)[3], dim(growth_env)[3], dim(regeneration_env)[3])
  if(!is.null(mortalityProcess$inputNN)) {
    input_dimensions[1] = mortalityProcess$inputNN
    runEnvMort = FALSE
  }
  if(!is.null(growthProcess$inputNN)) {
    input_dimensions[2] = growthProcess$inputNN
    runEnvGrowth = FALSE
  }
  if(!is.null(regenerationProcess$inputNN)) {
    input_dimensions[3] = regenerationProcess$inputNN
    runEnvReg = FALSE
  }

  out$input_dimensions = input_dimensions
  output = c(sp, sp, sp)
  if(!is.null(mortalityProcess$outputNN)) output[1] = mortalityProcess$outputNN
  if(!is.null(growthProcess$outputNN)) output[2] = growthProcess$outputNN
  if(!is.null(regenerationProcess$outputNN)) output[3] = regenerationProcess$outputNN
  out$output = output


  model = FINNModel$new(sp = sp,
                        env = input_dimensions,
                        output = output,
                        device = device,
                        hidden_growth = growthProcess$hidden,
                        hidden_mort = mortalityProcess$hidden,
                        hidden_reg = regenerationProcess$hidden,
                        runEnvGrowth = runEnvGrowth,
                        runEnvMort = runEnvMort,
                        runEnvReg = runEnvReg,
                        growthFunction = growthProcess$func,
                        mortalityFunction = mortalityProcess$func,
                        regenerationFunction = regenerationProcess$func,
                        competitionFunction = competitionProcess$func,
                        parGrowth = growthProcess$initSpecies,
                        parMort = mortalityProcess$initSpecies,
                        parReg = regenerationProcess$initSpecies,
                        parComp = competitionProcess$initSpecies,
                        parThetaReg = regenerationProcess$dispersion_parameter,
                        # parHeight = height,
                        parGrowthEnv = growthProcess$initEnv,
                        parMortEnv = mortalityProcess$initEnv,
                        parRegEnv = regenerationProcess$initEnv,
                        speciesPars_ranges = speciesPars_ranges,
                        patch_size_ha = patch_size)

  #### set parameters for optimization #####
  if(!mortalityProcess$optimizeSpecies) model$parMort$requires_grad_(FALSE)
  if(!growthProcess$optimizeSpecies) model$parGrowth$requires_grad_(FALSE)
  if(!regenerationProcess$optimizeSpecies) model$parReg$requires_grad_(FALSE)
  if(!competitionProcess$optimizeSpecies) model$parComp$requires_grad_(FALSE)
  if(!mortalityProcess$optimizeEnv) .n = lapply(model$nnMortEnv$parameters, function(p) p$requires_grad_(FALSE))
  if(!growthProcess$optimizeEnv) .n = lapply(model$nnGrowthEnv$parameters, function(p) p$requires_grad_(FALSE))
  if(!regenerationProcess$optimizeEnv) .n = lapply(model$nnRegEnv$parameters, function(p) p$requires_grad_(FALSE))
  # update parameter list
  model$parameter_to_r()
  model$update_parameters()



  if(is.null(bootstrap)) {

    model$fit(initCohort = init,
              X = list(torch::torch_tensor(mortality_env),
                       torch::torch_tensor(growth_env),
                       torch::torch_tensor(regeneration_env)),
              Y = response_T,
              disturbance = disturbance_T,
              patches = patches,
              batch_size = batchsize,
              epochs = epochs,
              learning_rate = lr,
              update_step = 1L,
              weights = weights,
              year_sequence = year_sequence,
              file = file,
              thin = thin,
              optimizer = optimizer)

    out$model = model
    out$init = init
    out$growthProcess = growthProcess
    out$mortalityProcess = mortalityProcess
    out$regenerationProcess = regenerationProcess
    out$competition = competitionProcess
    out$response = response
    out$disturbance = disturbance
    out$data = data
    out$env = list(mortality_env = mortality_env,
                   growth_env = growth_env,
                   regeneration_env = regeneration_env)
    class(out) = "finnModel"
    return(out)

  } else {

    if(NCPUs > 0.5) {

      cl = parallel::makeCluster(NCPUs)
      nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))
      parallel::clusterEvalQ(cl, {library(torch);library(FINN)})
      parallel::clusterExport(cl, varlist = ls(environment()) ,envir = environment())
      models_list  =
        parallel::parLapply(cl, 1:bootstrap, function(i) {

          indices = sample.int(sites, sites, replace = TRUE)

          if(device_old == "gpu") {
            myself = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
            dist = cbind(nodes,0:(NGPU-1))
            dev = as.integer(as.numeric(dist[which(dist[,1] %in% myself, arr.ind = TRUE), 2]))
            device_hardware = paste0("cuda:",dev)
          } else {
            device_hardware = "cpu"
          }

          response_T = torch::torch_cat(lapply(response, function(x) torch::torch_tensor(x, dtype=torch::torch_float32(), device="cpu")$unsqueeze(4)), 4)

          init$check()
          tmp_model = FINNModel$new(sp = sp,
                                    env = input_dimensions,
                                    output = output,
                                    device = device_hardware,
                                    hidden_growth = growthProcess$hidden,
                                    hidden_mort = mortalityProcess$hidden,
                                    hidden_reg = regenerationProcess$hidden,
                                    runEnvGrowth = runEnvGrowth,
                                    runEnvMort = runEnvMort,
                                    runEnvReg = runEnvReg,
                                    growthFunction = growthProcess$func,
                                    mortalityFunction = mortalityProcess$func,
                                    regenerationFunction = regenerationProcess$func,
                                    competitionFunction = competitionProcess$func,
                                    parGrowth = growthProcess$initSpecies,
                                    parMort = mortalityProcess$initSpecies,
                                    parReg = regenerationProcess$initSpecies,
                                    parComp = competitionProcess$initSpecies,
                                    # parHeight = height,
                                    parGrowthEnv = growthProcess$initEnv,
                                    parMortEnv = mortalityProcess$initEnv,
                                    parRegEnv = regenerationProcess$initEnv,
                                    speciesPars_ranges = speciesPars_ranges,
                                    patch_size_ha = patch_size)

          if(!mortalityProcess$optimizeSpecies) tmp_model$parMort$requires_grad_(FALSE)
          if(!growthProcess$optimizeSpecies) tmp_model$parGrowth$requires_grad_(FALSE)
          if(!regenerationProcess$optimizeSpecies) tmp_model$parReg$requires_grad_(FALSE)
          if(!competitionProcess$optimizeSpecies) tmp_model$parComp$requires_grad_(FALSE)
          if(!mortalityProcess$optimizeEnv) .n = lapply(tmp_model$nnMortEnv$parameters, function(p) p$requires_grad_(FALSE))
          if(!growthProcess$optimizeEnv) .n = lapply(tmp_model$nnGrowthEnv$parameters, function(p) p$requires_grad_(FALSE))
          if(!regenerationProcess$optimizeEnv) .n = lapply(tmp_model$nnRegEnv$parameters, function(p) p$requires_grad_(FALSE))
          # update parameter list
          tmp_model$parameter_to_r()
          tmp_model$update_parameters()

          tmp_model$fit(initCohort = init,
                        X = list(torch::torch_tensor(mortality_env[indices,,,drop=FALSE]),
                                 torch::torch_tensor(growth_env[indices,,,drop=FALSE]),
                                 torch::torch_tensor(regeneration_env[indices,,,drop=FALSE])),
                        Y = response_T,
                        disturbance = disturbance_T,
                        patches = patches,
                        batch_size = batchsize,
                        epochs = epochs,
                        learning_rate = lr,
                        update_step = 1L,
                        weights = weights,
                        optimizer = optimizer)

          out$model = tmp_model
          out$indices = indices
          out$init = init
          out$growthProcess = growthProcess
          out$mortalityProcess = mortalityProcess
          out$regenerationProcess = regenerationProcess
          out$competition = competitionProcess
          out$response = response
          out$disturbance = disturbance
          out$data = data
          out$env = list(mortality_env = mortality_env,
                         growth_env = growth_env,
                         regeneration_env = regeneration_env)
          class(out) = "finnModel"
          return(out)
        })
      parallel::stopCluster(cl)
      out$models_list = models_list
      .n = lapply(models_list, function(m) {
        m$model$check(device = device)
        m$init$check()
      })
      class(out) = "finnModelBootstrap"
    } else {
      models_list =
        lapply(1:bootstrap, function(i) {
          indices = sample.int(sites, sites, replace = TRUE)
          init$check()
          tmp_model = FINNModel$new(sp = sp,
                                    env = input_dimensions,
                                    output = output,
                                    device = device,
                                    hidden_growth = growthProcess$hidden,
                                    hidden_mort = mortalityProcess$hidden,
                                    hidden_reg = regenerationProcess$hidden,
                                    runEnvGrowth = runEnvGrowth,
                                    runEnvMort = runEnvMort,
                                    runEnvReg = runEnvReg,
                                    growthFunction = growthProcess$func,
                                    mortalityFunction = mortalityProcess$func,
                                    regenerationFunction = regenerationProcess$func,
                                    competitionFunction = competitionProcess$func,
                                    parGrowth = growthProcess$initSpecies,
                                    parMort = mortalityProcess$initSpecies,
                                    parReg = regenerationProcess$initSpecies,
                                    parComp = competitionProcess$initSpecies,
                                    # parHeight = height,
                                    parGrowthEnv = growthProcess$initEnv,
                                    parMortEnv = mortalityProcess$initEnv,
                                    parRegEnv = regenerationProcess$initEnv,
                                    speciesPars_ranges = speciesPars_ranges,
                                    patch_size_ha = patch_size)

          if(!mortalityProcess$optimizeSpecies) tmp_model$parMort$requires_grad_(FALSE)
          if(!growthProcess$optimizeSpecies) tmp_model$parGrowth$requires_grad_(FALSE)
          if(!regenerationProcess$optimizeSpecies) tmp_model$parReg$requires_grad_(FALSE)
          if(!competitionProcess$optimizeSpecies) tmp_model$parComp$requires_grad_(FALSE)
          if(!mortalityProcess$optimizeEnv) .n = lapply(tmp_model$nnMortEnv$parameters, function(p) p$requires_grad_(FALSE))
          if(!growthProcess$optimizeEnv) .n = lapply(tmp_model$nnGrowthEnv$parameters, function(p) p$requires_grad_(FALSE))
          if(!regenerationProcess$optimizeEnv) .n = lapply(tmp_model$nnRegEnv$parameters, function(p) p$requires_grad_(FALSE))
          # update parameter list
          tmp_model$parameter_to_r()
          tmp_model$update_parameters()

          tmp_model$fit(initCohort = init,
                        X = list(torch::torch_tensor(mortality_env[indices,,,drop=FALSE]),
                                 torch::torch_tensor(growth_env[indices,,,drop=FALSE]),
                                 torch::torch_tensor(regeneration_env[indices,,,drop=FALSE])),
                        Y = response_T,
                        disturbance = disturbance_T,
                        patches = patches,
                        batch_size = batchsize,
                        epochs = epochs,
                        learning_rate = lr,
                        update_step = 1L,
                        weights = weights,
                        optimizer = optimizer)

          out$model = tmp_model
          out$indices = indices
          out$init = init
          out$growthProcess = growthProcess
          out$mortalityProcess = mortalityProcess
          out$regenerationProcess = regenerationProcess
          out$competition = competitionProcess
          out$response = response
          out$disturbance = disturbance
          out$data = data
          out$env = list(mortality_env = mortality_env,
                         growth_env = growth_env,
                         regeneration_env = regeneration_env)
          class(out) = "finnModel"
          return(out)
        })
      out$models_list = models_list
      class(out) = "finnModelBootstrap"
    }
  }
  return(out)

}



#' Predict Forest Dynamics
#'
#' @param object object of class 'finnModel' created by `finn()`
#' @param init new `cohortMat$new()` object, if `NULL`, the initial cohortMat of the object is used
#' @param env new environment to simulate for
#' @param ... currently ignored
#'
#' @export
predict.finnModel = function(object, init = NULL, env = NULL, disturbance = NULL, ...) {
  object$model$check()
  object$init$check()

  if(is.null(env)) {
    mortality_env = object$env$mortality_env
    growth_env = object$env$growth_env
    regeneration_env = object$env$regeneration_env
  } else {
    mortality_env = extract_env(object$mortalityProcess, env)
    growth_env = extract_env(object$growthProcess, env)
    regeneration_env = extract_env(object$regenerationProcess, env)
  }

  disturbance_T = NULL
  if(!is.null(disturbance)) {
    disturbance = extract_env(list(formula=~0+intensity), disturbance)
    disturbance_T = torch::torch_tensor(disturbance, dtype=torch::torch_float32(), device="cpu")
  } else {
    disturbance = object$disturbance
    if(!is.null(disturbance)) {
      disturbance_T = torch::torch_tensor(disturbance, dtype=torch::torch_float32(), device="cpu")
    }
  }


  if(is.null(init)) {
    init = object$init
  }

  X = list(torch::torch_tensor(mortality_env),
           torch::torch_tensor(growth_env),
           torch::torch_tensor(regeneration_env))

  predictions =
    object$model$predict(
      dbh = init$dbh,
      trees = init$trees,
      species = init$species,
      env = X,
      disturbance = disturbance_T,
      patches = init$dbh$shape[2],
      debug = FALSE
    )

  return(list(wide = pred2DF(predictions, format = "wide"),
              long = pred2DF(predictions, format = "long")
  ))
}


#' Continue Training a `finnModel`
#'
#' This function continues the training process for an existing `finnModel` object. It fits the model for a specified number of epochs using the provided learning rate and batch size, updating the model's weights accordingly.
#'
#' @param object An object of class `finnModel`, typically obtained from the `finn()` function. This object contains the model, initial conditions, and data required for training.
#' @param epochs An integer specifying the number of epochs to train the model. Defaults to 20.
#' @param lr A numeric value specifying the learning rate for the training process. If `NULL`, the learning rate from the `object` is used.
#' @param batchsize An integer specifying the batch size to use during training. If `NULL`, the batch size from the `object` is used.
#'
#' @return The function returns an updated object of class `finnModel` with the model's training continued for the specified number of epochs.
#'
#' @details The function prepares the response variables and disturbance data, converts them into tensors, and then proceeds with the training process using the `fit` method of the `finnModel` object. The training is performed on the CPU.
#'
#' @examples
#' \dontrun{
#'   model <- finn(data, ...)
#'   model <- continue_fit(model, epochs = 50, lr = 0.01, batchsize = 32)
#' }
#'
#' @export
continue_fit = function(object, epochs = 20L, lr = NULL, batchsize = NULL, weights = NULL) {
  object$model$check()
  object$init$check()

  sp = object$sp

  options(na.action='na.pass')
  response = list(
    dbh = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+dbh), object$data[object$data$species==i,])), along = 3L),
    ba = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+ba), object$data[object$data$species==i,])), along = 3L),
    trees = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+trees), object$data[object$data$species==i,])), along = 3L),
    AL = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+AL), object$data[object$data$species==i,])), along = 3L),
    growth = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+growth), object$data[object$data$species==i,])), along = 3L),
    mort = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+mort), object$data[object$data$species==i,])), along = 3L),
    reg = abind::abind(lapply(1:sp, function(i) extract_env(list(formula=~0+reg), object$data[object$data$species==i,])), along = 3L)
  )
  options(na.action='na.omit')

  disturbance = object$disturbance
  disturbance_T = NULL
  if(!is.null(disturbance)) {
    disturbance = extract_env(list(formula=~0+intensity), disturbance)
    disturbance_T = torch::torch_tensor(disturbance, dtype=torch::torch_float32(), device="cpu")
  }

  response_T = torch::torch_cat(lapply(response, function(x) torch::torch_tensor(x, dtype=torch::torch_float32(), device="cpu")$unsqueeze(4)), 4)

  if(is.null(lr)) lr = object$lr
  if(is.null(batchsize)) batchsize=object$batchsize

  year_sequence = object$model$year_sequence

  object$model$fit(initCohort = object$init,
                   X = list(torch::torch_tensor(object$env$mortality_env),
                            torch::torch_tensor(object$env$growth_env),
                            torch::torch_tensor(object$env$regeneration_env)),
                   Y = response_T,
                   disturbance = disturbance_T,
                   patches = object$patches,
                   batch_size = batchsize,
                   epochs = epochs,
                   learning_rate = lr,
                   update_step = 1L,
                   weights = weights,
                   year_sequence = year_sequence)
  return(invisible(object))
}

#' Print finnModel
#'
#' print finnModel object created by `finn()`
#'
#' @param object object of class finnModel
#' @param ... not implemented
#'
#' @export
print.finnModel = function(object, ...) {
  object$model$check()
  cat("Estimated Parameters: \n")
  if(attr(object$model$parGrowth_r, "requires_grad")) {
    cat("\nGrowth Species Parameter: \n")
    print(object$model$parGrowthTR)
  }
  if(attr(object$model$parMort_r, "requires_grad")) {
    cat("\nMort Species Parameter: \n")
    print(object$model$parMortTR)
  }
  if(attr(object$model$parReg_r, "requires_grad")) {
    cat("\nGrowth Species Parameter: \n")
    print(object$model$parRegTR)
  }
  if(attr(object$model$parComp_r, "requires_grad")) {
    cat("\nGrowth Species Parameter: \n")
    print(object$model$parCompTR)
  }
  if(attr(object$model$parGrowthEnv_r[[1]], "requires_grad")) {
    cat("\nGrowth Environment Parameter: \n")
    if(length(object$model$parGrowthEnv)) print(object$model$parGrowthEnv[[1]] |> as.matrix())
  }
  if(attr(object$model$parMortEnv_r[[1]], "requires_grad")) {
    cat("\nMort Environment Parameter: \n")
    if(length(object$model$parMortEnv)) print(object$model$parMortEnv[[1]] |> as.matrix())
  }
  if(attr(object$model$parRegEnv_r[[1]], "requires_grad")) {
    cat("\nReg Environment Parameter: \n")
    if(length(object$model$parRegEnv)) print(object$model$parRegEnv[[1]] |> as.matrix())
  }
}


#' Summary of finnModel
#'
#' summary of finnModel object created by `finn()`
#'
#' @param object object of class finnModel
#' @param ... not implemented
#'
#' @export
summary.finnModel = function(object, ...) {
  print(object)
}


#' finnModel plot
#'
#' Plotting finnModel outputs (predictions or loss).
#'
#' @param x a model fitted by `finnModel`
#' @param ... ignored
#'
#' @return
#'
#' ggplot2 object
#'
#' @import graphics
#' @export
plot.finnModel = function(x,  plot = c("predictions", "loss"),...) {

}




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
#'
#' @return A list of class "process" containing the process definition and associated parameters.
#'
#' @examples
#' growth_process <- createProcess(formula = ~temperature + precipitation, func = growthFunction)
#'
#' @export
createProcess = function(formula = NULL, func, initSpecies = NULL, initEnv = NULL, hidden = NULL, optimizeSpecies = FALSE, optimizeEnv = TRUE, inputNN = NULL, outputNN = NULL, dispersion_parameter = NULL) {
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
#' @param sp An integer specifying the number of species in the simulation. Default is NULL.
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
                          disturbance = NULL,
                          mortalityProcess = NULL,
                          growthProcess = NULL,
                          regenerationProcess = NULL,
                          competitionProcess = NULL,
                          speciesPars_ranges = list(
                            parGrowth = rbind(
                              c(0.01, 0.99),
                              c(0, 0.5)
                            ),
                            parMort = rbind(
                              c(0.011,0.93),
                              c(0, 4)
                            ),
                            parReg = c(0.01, 0.99),
                            parComp = rbind(
                              c(0, 1),
                              c(0, 2)
                            )),
                            # parHeight = c(0.3, 0.7)
                          # height = NULL,
                          patches = 100L,
                          patch_size = 0.1,
                          sp = NULL,
                          init = NULL,
                          device = c("cpu", "gpu"),
                          parallel = FALSE,
                          NGPU = 1,
                          batchsize = NULL,
                          debug = FALSE
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

  disturbance_T = NULL
  if(!is.null(disturbance)) {
    disturbance = extract_env(list(formula=~0+intensity), disturbance)
    disturbance_T = torch::torch_tensor(disturbance, dtype=torch::torch_float32(), device="cpu")
  }


  sites = dim(mortality_env)[1] # TODO check that all env arrays have the same sites and timesteps
  if(is.null(batchsize)) batchsize = sites

  if(is.null(init)) {
    if(is.null(sp)) stop("Number of species (argument 'sp') must be specified!")
    init = CohortMat$new(dims = c(sites, patches, sp),
                         dbh = array(1, dim = c(sites, patches, sp)),
                         trees = array(1, dim = c(sites, patches, sp)),
                         sp = sp)
  }else{
    if(is.null(sp)) warning(
      paste("Number of species (argument 'sp') not specified, using the number of species in the initial cohort matrix: sp =", init$sp)
      )
    sp = init$sp
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

  runEnvMort = TRUE
  runEnvGrowth = TRUE
  runEnvReg = TRUE

  input_dimensions = c(dim(mortality_env)[3], dim(growth_env)[3], dim(regeneration_env)[3])
  if(!is.null(mortalityProcess$inputNN)) {
    input_dimensions[1] = mortalityProcess$inputNN
    runEnvMort = FALSE
  }
  if(!is.null(growthProcess$inputNN)) {
    input_dimensions[2] = growthProcess$inputNN
    runEnvGrowth = FALSE
  }
  if(!is.null(regenerationProcess$inputNN)) {
    input_dimensions[3] = regenerationProcess$inputNN
    runEnvReg = FALSE
  }

  out$input_dimensions = input_dimensions

  output = c(sp, sp, sp)
  if(!is.null(mortalityProcess$outputNN)) output[1] = mortalityProcess$outputNN
  if(!is.null(growthProcess$outputNN)) output[2] = growthProcess$outputNN
  if(!is.null(regenerationProcess$outputNN)) output[3] = regenerationProcess$outputNN
  out$output = output

  model = FINNModel$new(sp = sp,
                        # Can be forcefully overwritten if the NNs should be used for something else e.g. hbyrid modeling, replace one of the process function by a NN!
                        env = input_dimensions,
                        output = output,
                        device = device,
                        hidden_growth = growthProcess$hidden,
                        hidden_mort = mortalityProcess$hidden,
                        hidden_reg = regenerationProcess$hidden,
                        runEnvGrowth = runEnvGrowth,
                        runEnvMort = runEnvMort,
                        runEnvReg = runEnvReg,
                        growthFunction = growthProcess$func,
                        mortalityFunction = mortalityProcess$func,
                        regenerationFunction = regenerationProcess$func,
                        competitionFunction = competitionProcess$func,
                        parGrowth = growthProcess$initSpecies,
                        parMort = mortalityProcess$initSpecies,
                        parReg = regenerationProcess$initSpecies,
                        parThetaReg = regenerationProcess$dispersion_parameter,
                        # parHeight = height,
                        parComp = competitionProcess$initSpecies,
                        parGrowthEnv = growthProcess$initEnv,
                        parMortEnv = mortalityProcess$initEnv,
                        parRegEnv = regenerationProcess$initEnv,
                        speciesPars_ranges = speciesPars_ranges,
                        patch_size_ha = patch_size)



  if(NCPUs < -0.5) {
    predictions = model$predict(dbh = init$dbh,
                                trees = init$trees,
                                species = init$species,
                                env = list(torch::torch_tensor(mortality_env),
                                           torch::torch_tensor(growth_env),
                                           torch::torch_tensor(regeneration_env)),
                                disturbance = disturbance_T,
                                patches = patches,
                                debug = debug,
                                verbose = TRUE)


    out = list(long = pred2DF(predictions, "long"), wide = pred2DF(predictions, "wide"))

  } else {
    cl = parallel::makeCluster(NCPUs)
    nodes = unlist(parallel::clusterEvalQ(cl, paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')))

    batches = cbind(c(1, (seq(1, sites, by = batchsize))[-1]), c((seq(1, sites, by = batchsize)-1)[-1], sites)  )

    parallel::clusterEvalQ(cl, {library(torch);library(FINN)})
    parallel::clusterExport(cl, varlist = c("FINNModel","CohortMat","model", "init", "mortality_env", "growth_env", "regeneration_env", "batches", "patches", "nodes", "NGPU", "device_old"),envir = environment())

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
        long = pred2DF(pred,format = "long")
        wide = pred2DF(pred,format = "wide")
        long$site$siteID = (batches[i,1]:batches[i,2])[long$site$siteID]
        wide$site$siteID = (batches[i,1]:batches[i,2])[wide$site$siteID]
        tmp_out = list(long = long, wide = wide)
        return(tmp_out)
      })
    parallel::stopCluster(cl)
    out = list(long = list(site=data.table::rbindlist(lapply(predictions, function(p) p$long$site))),
               wide = list(site=data.table::rbindlist(lapply(predictions, function(p) p$wide$site)))
               )
  }

  out$model = model
  return(out)
}






