#' finn module
#'
#' @export
finn = nn_module(
  "finn",
  initialize = function(
    N_species,
    mortality_process = NULL,
    growth_process = NULL,
    regeneration_process = NULL,
    competition_process = NULL,
    recruits_dbh = 1.0
  ) {
    self$N_species = N_species
    self$recruits_dbh = recruits_dbh
    private$add_process(mortality_process, "mortality")
    private$add_process(growth_process, "growth")
    private$add_process(regeneration_process, "regeneration")
    private$add_process(competition_process, "competition")

  },

  #' @description
  #' Predicts the growth, mortality, and regeneration of trees based on the given inputs.
  #'
  #' The `forward` method generates predictions for tree growth, mortality, and regeneration for the specified species across different environmental conditions. It uses the initialized model parameters and can handle optional input tensors like diameter at breast height (dbh), number of trees, and species. If these are not provided, they will be initialized internally.
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
  forward = function(dbh = NULL,
                     trees = NULL,
                     species = NULL,
                     env = NULL,
                     y = NULL,
                     disturbance = NULL,
                     start_time = 1L,
                     pred_growth = NULL,
                     pred_mort = NULL,
                     pred_reg = NULL,
                     patches = 100L,
                     debug = FALSE,
                     update_step = 1L,
                     verbose = TRUE,
                     year_sequence = NULL){


    # if no cohorts exist initialize empty cohort array
    if(is.null(dbh)){
      cohorts = CohortMat(dims = c(env$shape[1], patches, self$sp), sp = self$sp)
      cohorts$to(device = self$device) # TODO: device must be set
      trees = cohorts$trees
      species = cohorts$species
      dbh = cohorts$dbh
    }

    # # repeat env if same env for the three processes
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
    sp = self$N_species

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

    # send data to correct devices (and dtypes) TODO: unncessary?
    dbh = torch_tensor(dbh, dtype=self$dtype, device=self$device)
    trees = torch_tensor(trees, dtype=self$dtype, device=self$device)
    species = torch_tensor(species, dtype=torch_int64(), device=self$device)

    # create cohort ids
    cohort_ids = torch_tensor(array(
      1:(prod(species$shape)+1),
      dim = species$shape), dtype=torch_int32(), device = self$device
    )

    # init Result tensors
    Result = lapply(1:7,function(tmp) torch::torch_zeros(list(sites, time, self$N_species), device=self$device))
    names(Result) =  c("dbh","ba", "trees", "growth", "mort", "reg", "r_mean_ha")

    # if debug modus, create empty lists
    if(debug) {
      Raw_cohort_results = list()
      Raw_cohort_ids = list()
      Raw_patch_results = list()
    }

    # if simulation mode, environmental predictions can be made a priori (it would interrupt the gradients in inference mode)
    if(is.null(y)) {
      if(!inherits(self$process_growth, "hybrid")) predGrowthGlobal = self$nn_growth(env[["growth"]])
      if(!inherits(self$process_mortality, "hybrid")) predMortGlobal = self$nn_mortality(env[["mort"]])
      if(!inherits(self$process_regeneration, "hybrid")) predRegGlobal = self$nn_regeneration(env[["reg"]])
    }

    # create process bar
    if(verbose) cli::cli_progress_bar(format = "Year: {cli::pb_current}/{cli::pb_total} {cli::pb_bar} ETA: {cli::pb_eta} ", total = time, clear = FALSE)
    for(i in 1:time){

      # In inference mode, make env predictions in each time step (to get the gradients)
      # otherwise, just take the i-th prediction
      if(!is.null(y)) {
        if(!inherits(self$process_growth, "hybrid")) pred_growth = self$nn_growth(env[["growth"]][,i,])
        if(!inherits(self$process_mortality, "hybrid")) pred_mort = self$nn_mortality(env[["mort"]][,i,])
        if(!inherits(self$process_regeneration, "hybrid")) pred_reg = self$nn_regeneration(env[["reg"]][,i,])
      } else {
        if(!inherits(self$process_growth, "hybrid")) pred_growth = predGrowthGlobal[,i,]
        if(!inherits(self$process_mortality, "hybrid")) pred_mort = predMortGlobal[,i,]
        if(!inherits(self$process_regeneration, "hybrid")) pred_reg = predRegGlobal[,i,]
      }

      # empty rate objects/tensors
      light = torch_zeros(list(sites, time,  dbh$shape[3]), device=self$device)
      g = torch_zeros(list(sites, time, dbh$shape[3]), device=self$device)
      m = torch_zeros(list(sites, time, dbh$shape[3]), device=self$device)
      r = torch_zeros(list(sites, time, dbh$shape[3]), device=self$device)
      if(debug) trees_before = torch::torch_zeros_like(g)

      # detach previous cohort objects (to interrupt the gradients)

      # if(!is.null(y)) {
      #   #for(j in 1:3) Result[[j]] = Result[[j]]$detach()
      #   # if period_length = NA on all sites, it is assumed that we need only yearly gradients
      #   if(as.numeric(y[,,,7]$isnan()$bitwise_not()$sum()) < 0.5) {
      #     dbh=dbh$detach()
      #     trees=trees$detach()
      #     species=species$detach()
      #     cohort_ids=cohort_ids$detach()
      #   }
      #
      # }

      # dbh=dbh$detach()
      # trees=trees$detach()
      # species=species$detach()
      # cohort_ids=cohort_ids$detach()

      # Apply disturbance
      if(!is.null(disturbance)) {
        trees = trees*disturbances_tens[,i,]$unsqueeze(3L)
      }

      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      ## Demographic processes ####
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
      # check if there is at least one cohort alive (otherwise just jump directly to regeneration?)
      if(dbh$shape[3] > 0.5){
        # calculate available light for each cohort
        ### Competition 1 ####
        light = self$competition_func(
          dbh = dbh,
          species = species,
          trees = trees,
          # parHeight = parHeight,
          parComp = self$par_competition,
          h = NULL,
          # minLight = self$minLight,
          patch_size_ha = self$patch_size_ha,
          ba = NULL,
          cohortHeights = NULL
        )

        if (light$gt(1.0)$sum()$item() > 0) {
          stop("Light values > 1")
        }

        ## Growth ####
        if(!inherits(self$process_growth, "hybrid")) pred = index_species(pred_growth, species)
        else pred = env[[1]][,i,]

        g = self$growth_func(
          dbh = dbh,
          species = species,
          parGrowth = self$par_growth,
          pred = pred,
          light = light,
          trees = trees
        )
        self$g = g
        dbh_growth = dbh*g
        dbh = dbh + dbh_growth

        ## Competition 2 ####
        light = self$competition_func(
          dbh = dbh,
          species = species,
          trees = trees,
          parComp = self$par_competition,
          h = NULL,
          # minLight = self$minLight,
          patch_size_ha = self$patch_size_ha,
          ba = NULL,
          cohortHeights = NULL
        )
        # cat("Second section  B2\n")

        ## Mortality ####
        if(!inherits(self$process_mortality, "hybrid")) pred = index_species(pred_mort, species)
        else pred = env[[2]][,i,]

        m = self$mortality_func(
          dbh = dbh,
          species = species,
          trees = trees + 0.001,
          parMort = self$par_mortality,
          pred = pred,
          light = light
        )

        trees_dead = binomial_from_gamma(torch::torch_clamp(trees+trees$le(0.5)$float()+0.01, min = 1.0) , torch::torch_clamp(m, 0.01, 0.99))*trees$ge(0.5)$float()
        trees_dead = trees_dead + trees_dead$round()$detach() - trees_dead$detach()
        trees_before = trees

        #.unsqueeze(3) # TODO check!
        #trees$sub_(m)$clamp_(min = 0.0)
        trees = torch_clamp(trees - trees_dead, min = 0) #### TODO if trees = 0 then NA...prevent!
      }

      ### Competition 3 ####
      # start reg
      AL_reg = self$competition_func( # must have dimension = n species in last dim
        dbh = dbh,
        species = species,
        trees = trees,
        parComp = self$par_competition,
        h = 1,
        # minLight = self$minLight,
        patch_size_ha = self$patch_size_ha,
        ba = NULL,
        cohortHeights = NULL
      )

      ### Regeneration ####
      if(!inherits(self$process_regeneration, "hybrid"))pred = pred_reg
      else pred = env[["reg"]][,i,]

      r_mean_ha = self$regeneration_func(species = species,
                                            parReg = self$par_regeneration[,1],
                                            pred = pred,
                                            light = AL_reg)

      r_mean_patch = r_mean_ha*self$patch_size_ha

      if(self$sample_regeneration){
        r = rnbinom_torch(r_mean_patch, self$par_theta_recruits)
      }else if(!self$sample_regeneration){
        r = r_mean_patch
      }

      r = r + r$round()$detach() - r$detach()

      # New recruits
      dbh_new = ((r-1+0.1)/1e-3)$sigmoid() + (self$recruits_dbh - 1.0) # TODO: check!!! --> when r 0 dann dbh = 0, ansonsten dbh = 1 dbh[r==0] = 0
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
        samples = vector("list", 3)
        names(samples) = c(
          "g*trees_before",
          "m*trees_before",
          "trees_before"
        )
        samples[[1]] = g*trees_before
        samples[[2]] = m*trees_before
        samples[[3]] = trees_before # original number of trees

        Results_tmp = replicate(length(samples), torch_zeros_like(Result[[1]][,i,]))
        # calculate the sum of all elements in samples over all patches for each species, patch and site
        tmp_res = aggregate_results(species, samples, Results_tmp, aggregation = "sum")

        # Aggregation [sites, patches, cohorts] --> [sites, species]:
        if(as.numeric(trees$sum() > 0)) {
          mask = tmp_res[[3]]$gt(0.5)
          # Growth
          Result[[4]][,i,][mask] = Result[[4]][,i,][mask]+tmp_res[[1]][mask]/tmp_res[[3]][mask] # summe dbh_growth / summe trees  #/cohort_counts[alive_species])
          # Mort
          Result[[5]][,i,][mask] = Result[[5]][,i,][mask]+tmp_res[[2]][mask]/tmp_res[[3]][mask] # summe dbh_growth / summe trees  #/cohort_counts[alive_species])
        }
      }
      # reg extra
      ## Regeneration count
      tmp_res = aggregate_results(species_new, list(r), list(torch::torch_zeros(Result[[1]][,i,]$shape[1], sp, device = self$device )))
      Result[[6]][,i,] = Result[[6]][,i,]$add(tmp_res[[1]])/patches

      ## Regeneration rate mean
      tmp_res = aggregate_results(species_new, list(r_mean_ha), list(torch::torch_zeros(Result[[1]][,i,]$shape[1], sp, device = self$device )))
      r_mean_ha = tmp_res[[1]]/patches
      Result[[7]][,i,] = r_mean_ha
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

        #})
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
        )
        Raw_cohort_ids[[i]] = torch::as_array(cohort_ids$cpu())
      }

      loss = torch_zeros(7L, device = self$device)

      ##### Rework #######
      if(i > 0 && dbh$shape[3] != 0 && !is.null(y) && (i %% update_step == 0)) {
        if(i %in% year_sequence) {
          tmp_index = which(year_sequence %in% i, arr.ind = TRUE)
          # #dbh
          loss[1] = self$loss_dbh_func(y[, tmp_index,,1], Result[[1]][,i,] )
          # ba
          loss[2] = self$loss_ba_func(y[, tmp_index,,2], Result[[2]][,i,] )
          # counts
          loss[3] = self$loss_trees_func(y[,tmp_index,,3], Result[[3]][,i,])


          # growth rates - check for NA in period_length, if not, then accumulate gradients?
          # currently we assume that they are constant over sites
          accumulate_gradients = y[,tmp_index,,7] |> as.matrix()
          period = unique(accumulate_gradients[,1])
          if(!is.na(period)) {
            # loss[1] = self$loss_dbh_func(y[, tmp_index,,1], Result[[1]][,(i-period+1):(i),]$mean(2) )
            # # ba
            # loss[2] = self$loss_ba_func(y[, tmp_index,,2], Result[[2]][,(i-period+1):(i),]$mean(2) )
            # # counts
            # loss[3] = self$loss_trees_func(y[,tmp_index,,3], Result[[3]][,(i-period+1):(i),]$mean(2))

            #browser()
            loss[4] = self$loss_growth_func(y[,tmp_index,,4], Result[[4]][,(i-period+1):(i),]$mean(2) )
            # mort rates
            loss[5] = self$loss_mortality_func(y[,tmp_index,,5], Result[[5]][,(i-period+1):(i),]$mean(2))
            # reg rates ha
            loss[6] = self$loss_regeneration_func(y[,tmp_index,,6], Result[[7]][,(i-period+1):(i),]$sum(2) )#(Result[[7]][,(i-period+1),] - Result[[7]][,i,])$clamp(min = 0.0)  )
            self$obs_rec = y[,tmp_index,,6] |> as.matrix()
            self$pred_rec = Result[[7]][,(i-period+1):(i),]$sum(2) |> as.matrix()
            self$loss_raw = as.numeric(loss)
            loss$sum()$backward()
            for(j in 1:7) Result[[j]] = Result[[j]]$detach()
          } else {

            # loss[1] = self$loss_dbh_func(y[, tmp_index,,1], Result[[1]][,i,] )
            # # ba
            # loss[2] = self$loss_ba_func(y[, tmp_index,,2], Result[[2]][,i,] )
            # # counts
            # loss[3] = self$loss_trees_func(y[,tmp_index,,3], Result[[3]][,i,])

            loss[4] = self$loss_growth_func(y[,tmp_index,,4], Result[[4]][,i,])
            # mort rates
            loss[5] = self$loss_mortality_func(y[,tmp_index,,5], Result[[5]][,i,])
            # reg rates ha
            loss[6] = self$loss_regeneration_func(y[,tmp_index,,6], Result[[7]][,i,])
            self$loss_raw = as.numeric(loss)
            loss$sum()$backward()
            for(j in 1:7) Result[[j]] = Result[[j]]$detach()
          }

          dbh=dbh$detach()
          trees=trees$detach()
          species=species$detach()
          cohort_ids=cohort_ids$detach()

        }
        # for(j in 1:7) Result[[j]] = Result[[j]]$detach()
      }

      if(!is.null(y)) {
        #for(j in 1:3) Result[[j]] = Result[[j]]$detach()
        # if period_length = NA on all sites, it is assumed that we need only yearly gradients
        if(as.numeric(y[,,,7]$isnan()$bitwise_not()$sum()) < 0.5) for(j in 1:7) Result[[j]] = Result[[j]]$detach()

      }
      loss$detach_()

      if(verbose) cli::cli_progress_update()

    }

    #browser()

    names(Result) =  c("dbh","ba", "trees", "growth", "mort", "reg", "r_mean_ha")
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

  simulate = function(env,
                      disturbance = NULL,
                      patches = 100L,
                      patch_size = 0.1,
                      init_cohort = NULL,
                      batchsize = NULL,
                      device = c("cpu", "gpu"),
                      debug = FALSE
                      ) {
    device = match.arg(device)
    if(device == "gpu") device="cuda:0"
    self$device = device
    self$patch_size_ha = patch_size

    envs = private$extract_env_method(env)

    if(!is.null(disturbance)) {
      disturbance = extract_env(~0+intensity, disturbance)
    }

    # init networks if not yet done
    private$create_nn(self$process_mortality, "mortality", dim(envs$mortality_env)[3])
    private$create_nn(self$process_growth, "growth", dim(envs$growth_env)[3])
    private$create_nn(self$process_regeneration, "regeneration", dim(envs$regeneration_env)[3])

    if(is.null(batchsize)) batchsize = nrow(envs[[1]])
    sites =  nrow(envs[[1]])

    if(is.null(init_cohort)) {
      sp = self$N_species
      self$init_cohort = CohortMat(dims = c(sites, patches, sp),
                                   dbh = array(1, dim = c(sites, patches, sp)),
                                   trees = array(1, dim = c(sites, patches, sp)),
                                   sp = sp)
    } else {
      self$init_cohort = init_cohort
    }
    self$to(device = device)
    if(is.null(disturbance)) {
      data = tensor_dataset(torch_tensor(envs[[1]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[2]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[3]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_arange(1,dim(envs[[1]])[1])$to(dtype = torch_int64() , device=torch_device('cpu'))
      )
    } else {
      data = tensor_dataset(torch_tensor(envs[[1]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[2]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[3]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_arange(1, dim(envs[[1]])[1] )$to(dtype = torch_int64() , device=torch_device('cpu')),
                            torch_tensor(disturbance, dtype=self$dtype, device=torch_device('cpu'))
      )
    }

    self$eval()

    DataLoader = torch::dataloader(data, batch_size=batchsize, shuffle=FALSE, num_workers=0, pin_memory=TRUE, drop_last=FALSE)

    predictions_batch = list()
    coro::loop(for (b in DataLoader) {
      x_mort = b[[1]]$to(device = self$device, non_blocking=TRUE)
      x_growth = b[[2]]$to(device = self$device, non_blocking=TRUE)
      x_reg = b[[3]]$to(device = self$device, non_blocking=TRUE)
      ind = b[[4]]$to(device = self$device, non_blocking=TRUE)
      dist = NULL
      if(!is.null(disturbance)) dist = b[[5]]$to(device = self$device, non_blocking=TRUE)

      trees = self$init_cohort$trees$to(device = self$device, non_blocking=TRUE)[ind,]
      species = self$init_cohort$species$to(device = self$device, non_blocking=TRUE)[ind,]
      dbh = self$init_cohort$dbh$to(device = self$device, non_blocking=TRUE)[ind,]
      pred_tmp = self$forward(dbh = dbh,
                              trees = trees,
                              species = species,
                              env = list(x_mort, x_growth, x_reg),
                              disturbance = dist,
                              verbose = TRUE,
                              debug = debug)
      pred = pred_tmp[[1]]
      predictions_batch = append(predictions_batch, list(list(long = pred2DF(list(Predictions = pred), "long"), wide = pred2DF(list(Predictions = pred), "wide"))))

      })
    predictions =
      list(long = list(site = rbindlist( lapply(predictions_batch, function(X) X$long$site ))),
           wide = list(site = rbindlist( lapply(predictions_batch, function(X) X$wide$site )))
      )
    if(debug){
      predictions = list(wide = pred2DF(pred_tmp, format = "wide"), long = pred2DF(pred_tmp, format = "long"))
    }
    return(predictions)

  },


  fit = function(data = NULL,
                 env,
                 disturbance = NULL,
                 patches = 100L,
                 patch_size = 0.1,
                 init_cohort = NULL,
                 epochs = 20L,
                 lr = 0.01,
                 # # 1 -> dbh, 2 -> ba, 3 -> trees, 4 -> growth rates, 5 -> mort rates, 6 -> reg rates
                 loss = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "mse", regeneration = "nbinom"), #
                 weights = rep(1, 6),
                 optimizer = optim_ignite_adam,
                 batchsize = NULL,
                 device = c("cpu", "gpu"),
                 update_step = 1L,
                 start_time = 1L,
                 plot_progress = TRUE,
                 folder = NULL,
                 checkpoints = 100L,
                 shuffle = TRUE,
                 record_gradients = FALSE,
                 ...) {

    old_par = par(no.readonly = TRUE)
    if(!any(loss %in% c("mse", "gaussian", "poisson", "nbinom"))) stop("Loss not supported")

    if(is.null(data)) {
      print("No data. Switching into simulation modus...")
      return(self$predict(env = env, patches = patches, patch_size = patch_size, batchsize=batchsize, device=device))
    }

    # setup loss functions
    private$create_loss_functions(loss, weights)

    # setup data
    device = match.arg(device)
    if(device == "gpu") device="cuda:0"
    self$device = device
    self$patch_size_ha = patch_size

    options(na.action='na.pass')
    sp = self$N_species
    if(!"period_length" %in% colnames(data)) data$period_length = NA_real_
    response = list(
      dbh = abind::abind(lapply(1:self$N_species, function(i) extract_env(~0+dbh, data[data$species==i,])), along = 3L),
      ba = abind::abind(lapply(1:sp, function(i) extract_env(~0+ba, data[data$species==i,])), along = 3L),
      trees = abind::abind(lapply(1:sp, function(i) extract_env(~0+trees, data[data$species==i,])), along = 3L),
      growth = abind::abind(lapply(1:sp, function(i) extract_env(~0+growth, data[data$species==i,])), along = 3L),
      mort = abind::abind(lapply(1:sp, function(i) extract_env(~0+mort, data[data$species==i,])), along = 3L),
      reg = abind::abind(lapply(1:sp, function(i) extract_env(~0+reg, data[data$species==i,])), along = 3L),
      period_length = abind::abind(lapply(1:sp, function(i) extract_env(~0+period_length, data[data$species==i,])), along = 3L)
    )
    options(na.action='na.omit')

    Y = torch::torch_cat(lapply(response, function(x) torch::torch_tensor(x, dtype=torch::torch_float32(), device="cpu")$unsqueeze(4)), 4)

    envs = private$extract_env_method(env)

    if(!is.null(disturbance)) {
      disturbance = extract_env(~0+intensity, disturbance)
    }

    year_sequence = which(levels(as.factor(env$year)) %in% levels(as.factor(data$year)), arr.ind = TRUE)

    # init networks if not yet done
    private$create_nn(self$process_mortality, "mortality", dim(envs$mortality_env)[3])
    private$create_nn(self$process_growth, "growth", dim(envs$growth_env)[3])
    private$create_nn(self$process_regeneration, "regeneration", dim(envs$regeneration_env)[3])

    if(is.null(batchsize)) batchsize = nrow(envs[[1]])
    sites =  nrow(envs[[1]])

    if(is.null(init_cohort)) {
      sp = self$N_species
      self$init_cohort = CohortMat(dims = c(sites, patches, sp),
                                   dbh = array(1, dim = c(sites, patches, sp)),
                                   trees = array(1, dim = c(sites, patches, sp)),
                                   sp = sp)
    } else {
      self$init_cohort = init_cohort
      if(self$init_cohort$sp != sp) stop(paste("sp in cohort", self$init_cohort$sp, "does not match sp from data", sp))
      patches = dim(self$init_cohort$species_r)[2]
    }

    if(device == "gpu") {
      device = "cuda:0"
      if(!torch::cuda_is_available()) {
        cli::cli_text("GPU device not available...")
        device = "cpu"
      }
    }

    self$to(device = device)
    if(is.null(disturbance)) {
      data = tensor_dataset(torch_tensor(envs[[1]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[2]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[3]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(Y[,,,], dtype=self$dtype, device=torch_device('cpu')),
                            torch_arange(1,dim(envs[[1]])[1])$to(dtype = torch_int64() , device=torch_device('cpu'))
      )
    } else {
      data = tensor_dataset(torch_tensor(envs[[1]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[2]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(envs[[3]], dtype=self$dtype, device=torch_device('cpu')),
                            torch_tensor(Y[,,,], dtype=self$dtype, device=torch_device('cpu')),
                            torch_arange(1, dim(envs[[1]])[1] )$to(dtype = torch_int64() , device=torch_device('cpu')),
                            torch_tensor(disturbance, dtype=self$dtype, device=torch_device('cpu'))
      )
    }

    DataLoader = torch::dataloader(data, batch_size=batchsize, shuffle=shuffle, num_workers=0, pin_memory=TRUE, drop_last=TRUE)

    self$history = list()
    self$gradients = list()

    self$optimizer = optimizer(self$parameters, lr = lr, ...)

    cli::cli_progress_bar(format = "Epoch: {cli::pb_current}/{cli::pb_total} {cli::pb_bar} ETA: {cli::pb_eta} DBH: {dbh_l} BA: {ba_l} Trees: {trees_l} g: {g_l} m: {m_l} r: {r_l}", total = epochs, clear = FALSE)
    if(is.null(year_sequence)) year_sequence = 1:envs[[1]]$shape[2]

    self$train()

    for(epoch in 1:epochs){
      counter = 1
      coro::loop(for (b in DataLoader) {
        batch_loss = matrix(NA, nrow = 10000, ncol = 7L)
        x_mort =   b[[1]]$to(device = self$device, non_blocking=TRUE)
        x_growth = b[[2]]$to(device = self$device, non_blocking=TRUE)
        x_reg =    b[[3]]$to(device = self$device, non_blocking=TRUE)
        y =        b[[4]]$to(device = self$device, non_blocking=TRUE)
        ind =      b[[5]]$to(device = self$device, non_blocking=TRUE)
        dist = NULL
        if(!is.null(disturbance))
          dist =  b[[6]]$to(device = self$device, non_blocking=TRUE)

        self$optimizer$zero_grad()

        if (is.null(self$init_cohort)) {
          cohorts = CohortMat(dims = c(x$shape[1], patches, self$sp),
                                  sp = sp,
                                  device=self$device)
          trees = cohorts$trees
          species = cohorts$species
          dbh = cohorts$dbh
        } else {
          trees = self$init_cohort$trees$to(device = self$device, non_blocking=TRUE)[ind,]
          species = self$init_cohort$species$to(device = self$device, non_blocking=TRUE)[ind,]
          dbh = self$init_cohort$dbh$to(device = self$device, non_blocking=TRUE)[ind,]
        }

        pred_tmp = self$forward(dbh = dbh,
                                trees = trees,
                                species = species,
                                env = list(x_mort, x_growth, x_reg),
                                disturbance = dist,
                                start_time = start_time,
                                y = y,
                                update_step = update_step,
                                verbose = FALSE,
                                year_sequence = year_sequence)

        pred = pred_tmp[[1]]
        loss = pred_tmp[[2]]
        .null = torch::nn_utils_clip_grad_norm_(self$parameters, 2.0)

        self$optimizer$step()

        batch_loss[counter, ] =  as.numeric(loss$data()$cpu())
        counter <<- counter + 1
      })

      # too expensive with many parameters
      if(epoch %% checkpoints == 0) {
        self$param_history = c(self$param_history, list(lapply(self$parameters, function(p) as.matrix(p$cpu()))))
        if(record_gradients) {
          self$gradients = c(self$gradients, list(lapply(self$parameters, function(p) {
            if(prod(p$grad$shape) > 0.5) return(p$grad |> as.matrix())
          })))
        }
      }

      bl = colMeans(batch_loss, na.rm = TRUE)
      bl = round(bl, 5)
      dbh_l = bl[1]
      ba_l = bl[2]
      trees_l = bl[3]
      g_l = bl[4]
      m_l = bl[5]
      r_l = bl[6]

      #cat("Epoch: ", epoch, "Loss: ", bl, "\n")
      self$history[[epoch]] = colMeans(batch_loss, na.rm = TRUE)



      if(plot_progress) {
        losses = do.call(rbind, self$history)
        losses = apply(losses, 2, scales::rescale)

        par(mfrow = c(1,1))
          e = tryCatch({

             matplot(losses[,1:6], type = "l", lty = 1, col = c("#FA4D19", "#49777A", "#42DA5D","#FA8A19", "#000000", "#1992FA"),
                     lwd = 1.5,
                     las = 1, xlab = "Epochs", ylab = "Loss in %")
             legend("topright", pch = 15, bty = "n", col = c("#FA4D19", "#49777A", "#42DA5D","#FA8A19", "#000000", "#1992FA"),
                    legend = c("dbh", "ba", "trees", "growth", "mort", "reg"))


          # plot(1:epochs, y = losses[,i],
          #      xlab = "Epochs",
          #      ylab = paste0("Loss ", c("dbh", "ba", "trees", "growth", "mort", "regeneration")[i]),
          #      type = "l")
          }, error = function(e) e)
      }

      if(!is.null(folder)) {
        if(!dir.exists(folder)) dir.create(folder)
        if(epoch %% checkpoints == 0) torch::torch_save(self$clone(deep=FALSE), path = paste0(folder, "/", "epoch_", epoch, "model.pt"))
      }

      cli::cli_progress_update()

    }
    cli::cli_progress_done()


    if(torch::cuda_is_available()) torch::cuda_empty_cache()

    # ignore debugging method
    self$pred = list(long = pred2DF(list(Predictions = pred), "long"), wide = pred2DF(list(Predictions = pred), "wide"))
    do.call(par, old_par)

  },

  private = list(
    create_loss_functions = function(loss, weights) {
      for(l in 1:6) {
          tmp_loss = loss[l]
          tmp_loss_name = names(loss)[l]
          if(tmp_loss == "mse") {
            func =local({
              local_l = l
              local_weights = weights
              function(true, pred) {
                mask = true$isnan()$bitwise_not()
                if(as.logical(mask$max()$data())) return(torch::nnf_mse_loss(pred[mask], true[mask])$mean()*local_weights[local_l])
                else return(0.0)
              }
            })
          } else if(tmp_loss == "gaussian") {
            self[[paste0("par_loss_",tmp_loss_name, "_scale")]] = torch::nn_parameter(torch_tensor(0.55))
            func = local({
              local_tmp_loss_name = tmp_loss_name
              local_l = l
              local_weights = weights
              function(true, pred) {
                sigma_raw = self[[paste0("par_loss_",local_tmp_loss_name, "_scale")]]
                sigma = torch_log(1.0+torch_exp(sigma_raw) )
                mask = true$isnan()$bitwise_not()
                if(as.logical(mask$max()$data())) {
                  return(torch::distr_normal(pred[mask], scale = sigma)$log_prob(true[mask])$negative()$mean()*local_weights[local_l] + 0.001*(sigma**2))
                } else { return(0.0) }
              }
            })
          } else if(tmp_loss == "poisson") {
            func = local({
                local_l = l
                local_weights = weights
                function(true, pred) {
                  mask = true$isnan()$bitwise_not()
                  if(as.logical(mask$max()$data())) return(torch::distr_poisson(pred[mask])$log_prob(true[mask])$negative()$mean()*local_weights[local_l])
                  else return(0.0)
                }
              })
          } else if(tmp_loss == "nbinom") {
            # Trees can be also set to nbinom, recruits always (?) have theta, so create theta only for trees, otherwise use self$par_theta_recruits
            if(tmp_loss_name == "trees") {
              self[[paste0("par_loss_",tmp_loss_name, "_theta")]] = torch::nn_parameter(torch_tensor(0.5413))
              func =
                local({
                  local_tmp_loss_name = tmp_loss_name
                  local_l = l
                  local_weights = weights
                  function(true, pred) {
                    mask = true$isnan()$bitwise_not()
                    theta = 1.0/torch::nnf_softplus( self[[paste0("par_loss_",local_tmp_loss_name, "_theta")]] )
                    if(as.logical(mask$max()$data())) return(dnbinom_torch(pred[mask]+0.001, true[mask], theta)$mean()*local_weights[local_l])
                    else return(0.0)
                  }
                })
            } else if(tmp_loss_name == "regeneration") {
              func =
                local({
                  local_l = l
                  local_weights = weights
                  function(true, pred) {
                    mask = true$isnan()$bitwise_not()
                    theta = self$par_theta_recruits
                    if(as.logical(mask$max()$data())) return(dnbinom_torch(pred[mask]+0.001, true[mask], theta)$mean()*local_weights[local_l])
                    else return(0.0)
                  }
                })
            }
          }
          self[[paste0("loss_", tmp_loss_name, "_func")]] = func
      }
    },

    create_nn = function(obj, type, inputs) {
      hybrid = inherits(obj, "hybrid")
      if(is.null(self[[paste0("nn_", type)]])) {
        if(!hybrid) {
          nn =
            if(is.null(obj$NN)) build_NN(input_shape = inputs, output_shape = self$N_species, bias = TRUE, activation = "selu", hidden = obj$hidden, dropout = obj$dropout, last_activation = "linear")
            else nn

          if(!is.null(obj$initEnv)) {
            for(i in 1:length(nn$parameters)) nn$parameters[[i]]$set_data( obj$initEnv[[i]] )
          }

        } else {
          if(type == "mortality") {
            if(obj$transformer) {
              nn = hybrid_transformer(num_species = self$N_species,
                                      num_env_vars = inputs,
                                      dgtl_embedder_dim = 4L,
                                      max_len = 500L,
                                      emb_dim=obj$emb_dim,
                                      num_heads=1L,
                                      num_layers=obj$encoder_layers,
                                      dropout=obj$dropout,
                                      dim_feedforward = obj$dim_feedforward)
            } else {
               nn = hybrid_DNN(num_species = self$N_species,
                               num_env_vars = inputs+1, +1, # because of growth!
                               emb_dim=obj$emb_dim,
                               dropout=0.1,
                               hidden = obj$hidden)
            }
          }

          if(type == "growth") {
            if(obj$transformer) {
            nn = hybrid_transformer(num_species = self$N_species,
                                    num_env_vars = inputs,
                                    dgtl_embedder_dim = 3L,
                                    max_len = 500L,
                                    emb_dim=obj$emb_dim,
                                    num_heads=1L,
                                    num_layers=obj$encoder_layers,
                                    dropout=obj$dropout,
                                    dim_feedforward = obj$dim_feedforward)
            } else {
              nn = hybrid_DNN(num_species = self$N_species,
                              num_env_vars = inputs,
                              emb_dim=obj$emb_dim,
                              dropout=0.1,
                              hidden = obj$hidden)
            }

          }
        }

        if(!obj$optimizeEnv) {
          for(p in nn$parameters) p$requires_grad_(FALSE)
        }
        self[[paste0("nn_", type)]] = nn
      }
    },

    extract_env_method = function(env) {
      mortality_env = extract_env(self$mortality_formula, env)
      growth_env = extract_env(self$growth_formula, env)
      regeneration_env = extract_env(self$regeneration_formula, env)
      return(list(mortality_env=mortality_env, growth_env=growth_env, regeneration_env=regeneration_env))
    },

    add_process = function(obj, type = c("mortality", "growth", "regeneration", "competition")) {
      type = match.arg(type)

      hybrid = FALSE
      if(inherits(obj, "hybrid")) {
        hybrid = TRUE
        func = switch(type, growth= { growth_hybrid }, mortality = { mortality_hybrid })
        obj$func = func

      }

      if(is.null(obj)) {
        func = switch(type, mortality = { mortality }, growth = { growth }, regeneration = {regeneration}, competition = { competition })
        obj = createProcess(func = func)
      }

      private$setup_species_parameters(obj, type, hybrid)

      # Train env model or not
      self[[paste0("env_", type, "_optimized")]] = obj$optimizeEnv
      self[[paste0("process_", type)]] = obj

      # process specific parameters, e.g. theta for nbinom sampling of the recruits
      if(type == "regeneration") {
        self$sample_regeneration = obj$sample_regeneration
        self$par_theta_recruits = obj$dispersion_parameter
      }

    },
    setup_species_parameters = function(obj, type, hybrid) {
      self[[paste0(type, "_func")]] = private$set_environment(obj$func)
      self[[paste0(type, "_formula")]] = obj$formula
      self$register_buffer(paste0("par_", type, "_upper"), torch::torch_tensor(get_par_boundary(obj, type, upper = TRUE)))
      self$register_buffer(paste0("par_", type, "_lower"), torch::torch_tensor(get_par_boundary(obj, type, upper = FALSE)))

      self[[paste0("par_", type, "_optimized")]] = obj$optimizeSpecies

      # to keep things simple, in case of hybrid, we will just ignore the species parameters for now....TODO!

      # if null, random initialisierung
      # else recalclate from init values the required values
      if(is.null(obj$initSpecies)) {
        self[[paste0("par_", type)]] = init_species_parameters(type, self$N_species)
      } else {
        self[[paste0("par_", type)]] = obj$initSpecies
      }


    },
    set_environment = function(fn) {
      environment(fn) = self$.__enclos_env__
      return(fn)
    }

  ),
  active = list(
    # TODO: reduce code redundancy!
    par_mortality = function(value) {
      if(missing(value)) {
        return(forward(self$par_mortality_unconstrained, self$par_mortality_upper, self$par_mortality_lower))
      } else {
        init = torch::torch_tensor( backward(value,  self$par_mortality_upper, self$par_mortality_lower) )
        if(self$par_mortality_optimized) self$par_mortality_unconstrained = torch::nn_parameter(init)
        else self$register_buffer("par_mortality_unconstrained", init) # if not trainable, the parameter should be created via register_buffer
      }
    },
    par_growth = function(value) {
      if(missing(value)) {
        return(forward(self$par_growth_unconstrained, self$par_growth_upper, self$par_growth_lower))
      } else {
        init = torch::torch_tensor( backward(value,  self$par_growth_upper, self$par_growth_lower) )
        if(self$par_growth_optimized) self$par_growth_unconstrained = torch::nn_parameter(init)
        else self$register_buffer("par_growth_unconstrained", init)
      }
    },
    par_regeneration = function(value) {
      if(missing(value)) {
        return(forward(self$par_regeneration_unconstrained, self$par_regeneration_upper, self$par_regeneration_lower))
      } else {
        init = torch::torch_tensor( backward(value,  self$par_regeneration_upper, self$par_regeneration_lower) )
        if(self$par_regeneration_optimized) self$par_regeneration_unconstrained = torch::nn_parameter(init)
        else self$register_buffer("par_regeneration_unconstrained", init)
      }
    },
    par_competition = function(value){
      if(missing(value)) {
        return(forward(self$par_competition_unconstrained, self$par_competition_upper, self$par_competition_lower))
      } else {
        init = torch::torch_tensor( backward(value,  self$par_competition_upper, self$par_competition_lower) )
        if(self$par_competition_optimized) self$par_competition_unconstrained = torch::nn_parameter(init)
        else self$register_buffer("par_competition_unconstrained", init)
      }
    },

    par_theta_recruits = function(value) {
      if(missing(value)) {
        return(1.0/(torch::nnf_softplus(self$par_theta_recruits_raw)+0.0001))
      } else {
        value = 1.0 / value - 0.0001
        self$par_theta_recruits_raw = torch::nn_parameter( torch_tensor(log( exp(value) - 1.0 )) )
      }
    },

    par_mortality_r = function() {
      return( as.matrix(self$par_mortality))
    },
    par_growth_r = function() {
      return( as.matrix(self$par_growth))
    },
    par_regeneration_r = function() {
      return( as.matrix(self$par_regeneration))
    },
    par_competition_r = function() {
      return( as.matrix(self$par_competition))
    },
    parameters_r = function() {
      return(lapply(self$parameters, function(p) as.matrix(p)))
    }
  )
)


# gg = createProcess(func = growth, optimizeSpecies = TRUE)
#
# m = finn(N_species = 5, growth_process = gg, mortality_process = createProcess(func = mortality))
# library(data.table)
# sim_env_data <- data.table(expand.grid(siteID = 1:10, year = 1:5, temp = 0.9, precip = 0.7))
#
# m = finn(N_species = 5, growth_process = gg, mortality_process = createProcess(func = mortality))
# t = m$simulate(env = sim_env_data, batchsize = 10L)
#
# data = t$wide$site
# m = finn(N_species = 5, growth_process = gg, mortality_process = createProcess(func = mortality))
# m$fit(data = data, env = sim_env_data, lr = 0.001, epochs = 20L, folder = "checkpoints", checkpoints = 5L)




# period_length

# plot(t$wide$site$dbh)
#
# list(long = list(site = rbindlist( lapply(t, function(X) X$long$site ))),
#      wide = list(site = rbindlist( lapply(t, function(X) X$wide$site )))
#      )
# rbindlist( lapply(t, function(X) X$long$site ))
# rbindlist( lapply(t, function(X) X$wide$site ))
#
# nn2 = torch::torch_load( torch::torch_serialize(m) )
# nn2$simulate(sim_env_data)

#
# m$fit(X, Y, .....)
#
# m$fit(...)
#
# # von den Prozessen, brauchen wir:
# # function
# # parameters
# # data model
# # formula
# # limits/parameter ranges der parameter?!
#
#
#
#
#
#
#
# ff = FINN(Nspecies = 5)
# ff$mortality
#
# #
# A = 5
# ff = function() A
# library(torch)
# Test = nn_module("Test",
#                  initialize = function(ff = NULL) {
#                    self$ff = ff
#                    self$nn = nn_sequential(nn_linear(5, 5))
#                    self$par = nn_parameter(torch_tensor(5.0, requires_grad = TRUE))
#                    self$register_buffer("g", torch_tensor(5.0))
#                    self[["Test"]] = nn_parameter(torch_tensor(1.0))
#                    self$egal = NULL
#                    self$optimizer = torch::optim_ignite_adagrad(params = self$par)
#
#                  },
#                  forward = function(x) torch::torch_serialize(self$clone()),
#                  active = list(
#                    testPar = function(value) {
#                      if(missing(value)) {return(self$par+5.0)}
#                      else self$par = nn_parameter(torch_tensor(value))
#                    }
#                  )
# )
# #
# nn = Test(ff)
# torch_save(nn, path = "test.pt")
# nn$forward()
# class(nn)
#
# nn2 = torch_load("test.pt")
# nn2$parameters

# nn$parameters
# nn$test_NN
#
# nn$to(device = "mps")
#
# nn$optimizer$param_groups
# nn$test2 = nn_parameter(torch_tensor(3.0))
#
# nn$testPar = 10.0
# nn$testPar
#
# nn$parameters
# nn$g =torch_tensor(3.0, requires_grad = TRUE)
# nn$g
#
# nn$ff()
# nn$par
#
# nn2 = torch::torch_load( torch::torch_serialize(m) )
# nn2$parameters
# nn2$optimizer$param_groups
#
# torch::torch_save(nn, "T.RDS")
# nn2 = torch_load("T.RDS")
# nn2$ff()
# nn2$par
#
#
# binomial_coefficient = function(n, x){
#   return(torch::torch_exp(torch::torch_lgamma((n + 1.0)) -
#                             torch::torch_lgamma((x + 1.0)) -
#                             torch::torch_lgamma((n - x + 1.0))))
# }
#
# binomial_likelihood = function(x, n, p) {
#   binom_coeff = binomial_coefficient(n, x)
#   likelihood = binom_coeff * (p ** x) * ((1 - p) ** (n - x))
#   return(likelihood$log())
# }
# binomial_likelihood(7, 10, 0.3)
# dbinom(7, 10, 0.3, log = TRUE)
