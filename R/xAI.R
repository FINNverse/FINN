#' Accumulated local effect plots
#'
#' @description
#' Calculates accumulated local effects (ALE) for the three processes
#'
#' @param model (`finn_class`)\cr Model object of class `finn_class`.
#' @param env (`data.table|data.frame`)\cr Environmental covariates for which the ALE should be calculated.
#' @param init_cohort (`CohortMat`)\cr Initial cohort of class `CohortMat`, if `NULL` then simulated from bare ground.
#' @param ... \cr Not supported yet.
#'
#' @export
ALE = function(model, env, init_cohort = NULL, ...) {

  env_dt = as.data.table(env)

  model$raw_g = NULL
  model$raw_m = NULL
  model$raw_r = NULL

  model$record_raws = TRUE
  sim = model |> simulateForest(env = env_dt, init_cohort = init_cohort)
  model$record_raws = FALSE

  out = list()

  years = env_dt$year |> unique() |> length()
  sitesTotal = env_dt$year |> unique() |> length()
  patches = dim(model$raw_g[[1]])[2]
  patch = 1
  time = 1

  for(process in c("mortality", "growth", "regeneration")) {
    print(process)
    df = data.frame()

    for(timeInd in 1:years) {
      time = (env_dt$year |> unique())[timeInd]
      for(siteInd in 1:sitesTotal) {
        site = (env_dt$siteID |> unique())[siteInd]
        for(patch in 1:patches) {

          #### env
          env_tmp = as.matrix(model.matrix(model$process_growth$formula, data = env_dt[siteID==site & year == time][1,-(1:2),drop=F]))
          cols_env = colnames(env_tmp)
          tmp =
            if(process == "growth") {
              model$raw_g[[timeInd]][siteInd, patch,,]
            } else if(process == "mortality") {
              model$raw_m[[timeInd]][siteInd, patch,,]
            } else if(process == "regeneration") {
              model$raw_r[[timeInd]][siteInd, patch,,]
            }
          if(is.vector(tmp)) tmp = matrix(tmp, nrow = 1L)
          tmp = cbind(tmp, do.call(rbind, lapply(1:nrow(tmp), function(i) env_tmp)))
          tmp = data.frame(tmp)
          if(process == "growth") colnames(tmp) = c("dbh", "light", "trees", "species", cols_env)
          if(process == "mortality") colnames(tmp) = c("dbh","growth", "light", "trees", "species", cols_env)
          if(process == "regeneration") colnames(tmp) = c("light", cols_env)
          tmp$time = time
          if(!process == "regeneration") tmp = tmp[tmp$trees > 0.5,]
          df = rbind(df, tmp)
        }
      }
    }
    df = as.data.table(df)
    if(!process == "regeneration") df = df[df$dbh > 0.5,]

    df_results_ale_process = data.frame()
    internal_df = data.frame()

    if(!process == "regeneration") {

      for(i in unique(df$species)) {
        SPECIES = i
        df_tmp = df[df$species == SPECIES,]
        if(nrow(df_tmp) > 1) {

          sites = nrow(df_tmp)
          dbh = torch::torch_tensor(array(df_tmp$dbh, dim = c(sites, 1, 1)), dtype = torch::torch_float32(), device = model$device, requires_grad = TRUE)
          if(process == "mortality") growthV = torch::torch_tensor(array(df_tmp$growth, dim = c(sites, 1, 1)), dtype = torch::torch_float32(), device = model$device, requires_grad = TRUE)
          trees = torch::torch_tensor(array(df_tmp$trees, dim = c(sites, 1, 1)), dtype = torch::torch_float32(), device = model$device, requires_grad = TRUE)
          light = torch::torch_tensor(array(df_tmp$light, dim = c(sites, 1, 1)), dtype = torch::torch_float32(), device = model$device, requires_grad = TRUE)
          species = torch::torch_tensor(array(SPECIES, dim = c(sites, 1, 1)), dtype = torch::torch_long(), device = model$device, requires_grad = FALSE)
          pred = torch_tensor(as.matrix(df_tmp[,..cols_env]), device = model$device, requires_grad = TRUE)
          pred2 = pred
          ## TODO
          if(process == "growth") {
            if(!inherits(model$process_growth, "hybrid")) {
              pred2 = model$nn_growth(pred2)
              pred2 = FINN::index_species(pred2, species)
            }
          }
          if(process == "mortality") {
            if(!inherits(model$process_mortality, "hybrid")) {
              pred2 = model$nn_mortality(pred2)
              pred2 = FINN::index_species(pred2, species)
            }
          }
          process_output =
            if(process == "growth") {
              model$growth_func(dbh = dbh, trees = trees, light = light, species = species, pred = pred2, parGrowth = model$par_growth)
            } else if(process == "mortality") {
              model$mortality_func(dbh = dbh, trees = trees, light = light, species = species, pred = pred2, parMort = model$par_mortality, growth = growthV)
            }

          vars =
            if(process == "mortality") {
              list(dbh, growthV, light)
            } else {
              list(dbh, light)
            }
          stand_vars_grad = sapply(vars, function(var) torch::autograd_grad(process_output[,1,1], inputs = var, retain_graph = TRUE,
                                                                            grad_outputs = torch_ones_like(process_output[,1,1]))[[1]][,1,1] |> as.matrix() )
          env_grad = torch::autograd_grad(process_output[,1,1], inputs = pred, retain_graph = TRUE, grad_outputs = torch_ones_like(process_output[,1,1]))[[1]]|> as.matrix()
          grads = cbind(stand_vars_grad, env_grad)
          df_tmp[,time:=NULL]
          df_tmp[,species:=NULL]
          df_tmp[,trees:=NULL]
          ales = ALE_ce(as.data.frame(df_tmp), grads, predictions = process_output[,1,1] |> as.numeric())
          ales$species = i
          copy_df = as.data.frame(df_tmp)
          copy_df$pred = process_output[,1,1] |> as.numeric()
          copy_df$species = i
          internal_df = rbind(internal_df, copy_df)
          df_results_ale_process = rbind(df_results_ale_process, ales)
        }
      }
    } else {
      df_tmp = df

      if(nrow(df_tmp) > 1) {

        sites = nrow(df_tmp)
        light = torch::torch_tensor(array(df_tmp$light, dim = c(sites, 1, 1)), dtype = torch::torch_float32(), device = model$device, requires_grad = TRUE)
        pred = torch_tensor(as.matrix(df_tmp[,..cols_env]), device = model$device, requires_grad = TRUE)
        pred2 = pred
        if(process == "regeneration") {
          if(!inherits(model$process_mortality, "hybrid")) {
            pred2 = model$nn_regeneration(pred2)
          }
        }

        process_output =
          if(process == "regeneration") {
            species = torch::torch_tensor(array(1:model$N_species, dim = c(sites, 1, model$N_species)), dtype = torch::torch_long(), device = model$device, requires_grad = FALSE)
            model$regeneration_func(light = light, species = species, parReg = model$par_regeneration[,1], pred = pred2)
          }

        vars = list(light)
        df_results_ale_process = data.frame()
        df_tmp[,time:=NULL]
        for(J in 1:model$N_species) {
          stand_vars_grad = sapply(vars, function(var) torch::autograd_grad(process_output[,1,J], inputs = var, retain_graph = TRUE,
                                                                            grad_outputs = torch_ones_like(process_output[,1,J]))[[1]][,1,1] |> as.matrix() )

          env_grad = torch::autograd_grad(process_output[,1,J], inputs = pred, retain_graph = TRUE, grad_outputs = torch_ones_like(process_output[,1,1]))[[1]]|> as.matrix()
          grads = cbind(stand_vars_grad, env_grad)

          ales = ALE_ce(as.data.frame(df_tmp), grads, predictions = process_output[,1,J] |> as.numeric())
          ales$species = J
          df_results_ale_process = rbind(df_results_ale_process, ales)
        }
      }

    }
    out[[process]] = df_results_ale_process
    #out[[paste0(process, "internal", collapse = "_")]] = internal_df
  }
  class(out) = c("FINNale")
  return(out)
}



# ALE_ce = function(X, ce, predictions = NULL, center = FALSE) {
#   ales =
#     lapply(1:ncol(X), function(i) {
#       effs = ce[,i]
#       data = X[,i]
#       ord = order(data)
#
#       x_sorted <- data[ord]
#       g_sorted <- effs[ord]
#
#       if(!is.null(predictions)) {
#         intercept = predictions[ord][1]
#       } else {
#         intercept = 0
#       }
#
#       ale_vals <- cumsum( (g_sorted[-1] + g_sorted[-length(g_sorted)]) / 2 * diff(x_sorted) )
#       if(!center) ale_vals <- ale_vals - mean(ale_vals)
#       ale_vals = ale_vals + intercept
#       data.frame(x = x_sorted[-1], ale = ale_vals, var = colnames(X)[i])
#     })
#   return(do.call(rbind, ales))
# }
ALE_ce <- function(X, ce, predictions = NULL, center = FALSE) {
  stopifnot(nrow(X) == nrow(ce), ncol(X) == ncol(ce))
  vars <- colnames(X)

  ales <- lapply(seq_len(ncol(X)), function(j) {
    xj <- X[, j]
    gj <- ce[, j]

    # 1) unique x values
    ux <- sort(unique(xj))
    n_ux <- length(ux)

    if (n_ux == 1L) {
      return(data.frame(
        x   = ux,
        ale = 0,
        var = if (!is.null(vars)) vars[j] else j
      ))
    }
    g_mean <- vapply(ux, function(val) {
      idx <- which(xj == val)
      mean(gj[idx], na.rm = TRUE)
    }, numeric(1))

    # trapezoidal integration
    dx      <- diff(ux)
    g_mid   <- (g_mean[-1] + g_mean[-length(g_mean)]) / 2
    increments <- g_mid * dx

    ale_vals <- c(0, cumsum(increments)) #+ mean(predictions)
    ale_vals = ale_vals - mean(ale_vals) + mean(predictions)
    data.frame(
      x   = ux,
      ale = ale_vals,
      var = if (!is.null(vars)) vars[j] else j
    )
  })
  do.call(rbind, ales)
}

