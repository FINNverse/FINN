ALE_ce = function(X, ce) {
  ales =
    lapply(1:ncol(X), function(i) {
      effs = ce[,i]
      data = X[,i]
      ord = order(data)

      x_sorted <- data[ord]
      g_sorted <- effs[ord]

      # simple trapezoidal integration
      ale_vals <- cumsum( (g_sorted[-1] + g_sorted[-length(g_sorted)]) / 2 * diff(x_sorted) )
      ale_vals <- ale_vals - mean(ale_vals)
      data.frame(x = x_sorted[-1], ale = ale_vals, var = colnames(X)[i])
    })
  return(do.call(rbind, ales))
}

.finn_ce_key = function(env, init_cohort, sim_seed, env_autoscale, fit_id) {
  env_fp = if(is.null(env)) "NULL" else {
    num = vapply(env, is.numeric, logical(1))
    paste(nrow(env), ncol(env), paste(names(env), collapse = ","),
          format(sum(as.matrix(env[, num, with = FALSE]), na.rm = TRUE), digits = 15), sep = "|")
  }
  ic_fp = if(is.null(init_cohort)) "bare" else paste(class(init_cohort), collapse = ",")
  paste(env_fp, ic_fp, if(is.null(sim_seed)) "noseed" else sim_seed,
        isTRUE(env_autoscale), if(is.null(fit_id)) "nofit" else fit_id, sep = "##")
}

.finn_conditional_effects = function(model, env, init_cohort = NULL, env_autoscale = TRUE, sim_seed = 42L) {

  if(!is.null(sim_seed)) FINN.seed(sim_seed)

  env_dt   = as.data.table(env)
  scaling  = model$env_scaling
  env_pred = if(isTRUE(env_autoscale) && !is.null(scaling)) apply_env_scaling(env_dt, scaling) else env_dt

  model$raw_g = NULL
  model$raw_m = NULL
  model$raw_r = NULL
  model$record_raws = TRUE
  if(!isTRUE(env_autoscale)) model$env_scaling = NULL
  sim = model$simulate(env_dt, init_cohort = init_cohort)
  model$env_scaling = scaling
  model$record_raws = FALSE

  # sorted, not order of appearance: model$raw_* are recorded in the simulation's
  # own (sorted) site and year order, so an unsorted `env` paired each stand state
  # with another site's or year's predictors (#72, B35)
  env_years  = sort(unique(env_dt$year))
  env_sites  = sort(unique(env_dt$siteID))
  years      = length(env_years)
  sitesTotal = length(env_sites)
  patches    = dim(model$raw_g[[1]])[2]

  processes = list()

  for(process in c("mortality", "growth", "regeneration")) {
    raws = switch(process, growth = model$raw_g, mortality = model$raw_m, regeneration = model$raw_r)
    if (length(raws) < years) {
      # a user's process function in a model saved before FINN recorded the inputs of
      # every process itself (or read back with torch_load(), which keeps the methods
      # it was saved with): no inputs to differentiate at, so this process is left out
      message(sprintf("conditional effects: the %s process recorded no inputs during the simulation (a custom process function in a model saved with an earlier FINN version); it is left out.", process))
      next
    }
    proc_formula = switch(process,
                          growth       = model$process_growth$formula,
                          mortality    = model$process_mortality$formula,
                          regeneration = model$process_regeneration$formula)
    df = data.frame()

    for(timeInd in 1:years) {
      time = env_years[timeInd]
      for(siteInd in 1:sitesTotal) {
        site = env_sites[siteInd]
        for(patch in 1:patches) {
          env_tmp  = as.matrix(model.matrix(proc_formula, data = env_pred[siteID==site & year == time][1,-(1:2),drop=FALSE]))
          cols_env = colnames(env_tmp)
          tmp = switch(process,
                       growth       = model$raw_g[[timeInd]][siteInd, patch,,],
                       mortality    = model$raw_m[[timeInd]][siteInd, patch,,],
                       regeneration = model$raw_r[[timeInd]][siteInd, patch,,])
          if(is.vector(tmp)) tmp = matrix(tmp, nrow = 1L)
          tmp = cbind(tmp, do.call(rbind, lapply(1:nrow(tmp), function(i) env_tmp)))
          tmp = data.frame(tmp)
          if(process == "growth")       colnames(tmp) = c("dbh", "light", "trees", "species", cols_env)
          if(process == "mortality")    colnames(tmp) = c("dbh","growth", "light", "trees", "species", cols_env)
          if(process == "regeneration") colnames(tmp) = c("light", cols_env)
          tmp$time = time
          if(!process == "regeneration") tmp = tmp[tmp$trees > 0.5,]
          df = rbind(df, tmp)
        }
      }
    }
    df = as.data.table(df)
    if(!process == "regeneration") df = df[df$dbh > 0.5,]

    env_vars = setdiff(cols_env, "(Intercept)")
    sp_list  = list()

    if(!process == "regeneration") {
      for(i in unique(df$species)) {
        SPECIES = i
        df_tmp  = df[df$species == SPECIES,]
        if(nrow(df_tmp) > 1) {
          sites = nrow(df_tmp)
          dbh   = torch::torch_tensor(array(df_tmp$dbh,   dim = c(sites, 1, 1)), dtype = torch::torch_float32(), device = model$device, requires_grad = TRUE)
          if(process == "mortality") growthV = torch::torch_tensor(array(df_tmp$growth, dim = c(sites, 1, 1)), dtype = torch::torch_float32(), device = model$device, requires_grad = TRUE)
          trees = torch::torch_tensor(array(df_tmp$trees, dim = c(sites, 1, 1)), dtype = torch::torch_float32(), device = model$device, requires_grad = TRUE)
          light = torch::torch_tensor(array(df_tmp$light, dim = c(sites, 1, 1)), dtype = torch::torch_float32(), device = model$device, requires_grad = TRUE)
          species = torch::torch_tensor(array(SPECIES, dim = c(sites, 1, 1)), dtype = torch::torch_long(), device = model$device, requires_grad = FALSE)
          pred = torch_tensor(as.matrix(df_tmp[,..cols_env]), device = model$device, requires_grad = TRUE)
          pred2 = pred
          if(process == "growth" && !inherits(model$process_growth, "hybrid")) {
            pred2 = model$nn_growth(pred2); pred2 = index_species(pred2, species)
          }
          if(process == "mortality" && !inherits(model$process_mortality, "hybrid")) {
            pred2 = model$nn_mortality(pred2); pred2 = index_species(pred2, species)
          }
          process_output = switch(process,
                                  growth    = model$growth_func(dbh = dbh, trees = trees, light = light, species = species, pred = pred2, parGrowth = model$par_growth),
                                  mortality = model$mortality_func(dbh = dbh, trees = trees, light = light, species = species, pred = pred2, parMort = model$par_mortality, growth = growthV))

          vars = if(process == "mortality") list(dbh, growthV, light) else list(dbh, light)
          stand_vars_grad = sapply(vars, function(var) torch::autograd_grad(process_output[,1,1], inputs = var, retain_graph = TRUE,
                                                                            grad_outputs = torch_ones_like(process_output[,1,1]))[[1]][,1,1] |> as.matrix())
          env_grad = torch::autograd_grad(process_output[,1,1], inputs = pred, retain_graph = TRUE, grad_outputs = torch_ones_like(process_output[,1,1]))[[1]] |> as.matrix()
          grads = cbind(stand_vars_grad, env_grad)

          df_tmp[,time:=NULL]; df_tmp[,species:=NULL]; df_tmp[,trees:=NULL]
          X = as.data.frame(df_tmp)
          colnames(grads) = colnames(X)

          sp_list[[as.character(i)]] = list(X = X, grads = grads, pred = as.numeric(process_output[,1,1]),
                                            cols_env = cols_env, env_vars = env_vars)
        }
      }
    } else {
      df_tmp = df
      if(nrow(df_tmp) > 1) {
        sites = nrow(df_tmp)
        light = torch::torch_tensor(array(df_tmp$light, dim = c(sites, 1, 1)), dtype = torch::torch_float32(), device = model$device, requires_grad = TRUE)
        pred  = torch_tensor(as.matrix(df_tmp[,..cols_env]), device = model$device, requires_grad = TRUE)
        pred2 = pred
        if(!inherits(model$process_regeneration, "hybrid")) pred2 = model$nn_regeneration(pred2)
        species = torch::torch_tensor(array(1:model$N_species, dim = c(sites, 1, model$N_species)), dtype = torch::torch_long(), device = model$device, requires_grad = FALSE)
        process_output = model$regeneration_func(light = light, species = species, parReg = model$par_regeneration[,1], pred = pred2)

        df_tmp[,time:=NULL]
        X = as.data.frame(df_tmp)
        for(J in 1:model$N_species) {
          stand_vars_grad = sapply(list(light), function(var) torch::autograd_grad(process_output[,1,J], inputs = var, retain_graph = TRUE,
                                                                                   grad_outputs = torch_ones_like(process_output[,1,J]))[[1]][,1,1] |> as.matrix())
          env_grad = torch::autograd_grad(process_output[,1,J], inputs = pred, retain_graph = TRUE, grad_outputs = torch_ones_like(process_output[,1,1]))[[1]] |> as.matrix()
          grads = cbind(stand_vars_grad, env_grad)
          colnames(grads) = colnames(X)
          sp_list[[as.character(J)]] = list(X = X, grads = grads, pred = as.numeric(process_output[,1,J]),
                                            cols_env = cols_env, env_vars = env_vars)
        }
      }
    }
    processes[[process]] = sp_list
  }

  out = list(processes = processes, env = env_dt, env_autoscale = env_autoscale, sim_seed = sim_seed, key = NULL)
  class(out) = "FINNconditionalEffects"
  out
}

.finn_train_ic = function(model) {
  ic = model$train_init_cohort
  # a model read back with torch_load() holds the cohorts with dead tensors:
  # rebuilt from their R copy (#77)
  if(is.null(ic)) NULL else restore_init_cohort(ic$cohort)
}

.finn_get_ce = function(model, env = NULL, init_cohort = NULL, env_autoscale = TRUE, sim_seed = 42L) {
  # the key is a fingerprint (row and column count, names, column sums), so two
  # different inputs can share it; the cache is therefore used only for the
  # default inputs, where there is nothing to collide with (#72, B36)
  own_inputs = !is.null(env) || !is.null(init_cohort)
  if(is.null(env))         env         = model$train_env
  if(is.null(env))         stop("No `env` supplied and no training env is cached. Pass `env=` or fit the model first.")
  if(is.null(init_cohort)) init_cohort = .finn_train_ic(model)

  if(own_inputs) return(.finn_conditional_effects(model, env, init_cohort, env_autoscale, sim_seed))

  key = .finn_ce_key(as.data.table(env), init_cohort, sim_seed, env_autoscale, model$fit_id)
  if(!is.null(model$conditional_effects) && identical(model$conditional_effects$key, key)) {
    return(model$conditional_effects)
  }
  ce = .finn_conditional_effects(model, env, init_cohort, env_autoscale, sim_seed)
  ce$key = key
  model$conditional_effects = ce
  model$ale = NULL
  ce
}

#' Accumulated local effect plots
#'
#' @description
#' Calculates accumulated local effects (ALE) for the three processes
#'
#' @param model (`finn_class`)\cr Model object of class `finn_class`.
#' @param env (`data.table|data.frame`)\cr Environmental covariates for which the ALE should be
#'   calculated. If `NULL` (default) the training `env` cached by [fit()] is used.
#' @param init_cohort (`CohortMat`)\cr Initial cohort of class `CohortMat`. If `NULL` (default)
#'   the training init_cohort cached by [fit()] is used (or bare ground if none was cached).
#' @param env_autoscale (`logical(1)`)\cr If `TRUE` (default) `env` is assumed to be on the
#'   raw (unscaled) scale and the model's stored `env_scaling` is applied internally before the
#'   effects are computed, mirroring how the model was fitted (see [fit()]'s `env_autoscale`).
#'   Set `FALSE` if `env` is already on the model scale, or for a model fitted without autoscaling.
#' @param sim_seed (`integer(1)`)\cr Seed applied via [FINN.seed()] before the state-harvesting
#'   simulation, so the (stochastic) conditional effects / ALE are reproducible and cacheable.
#'   `NULL` disables seeding.
#' @param plot (`logical(1)`)\cr If `TRUE` (default) an ALE plot is drawn via [plot.FINNale()]
#'   (rows = processes, columns = environmental predictors, one coloured line per species).
#' @param process (`character`|`NULL`)\cr If given (one of `"growth"`, `"mortality"`,
#'   `"regeneration"`) only that process is plotted. `NULL` (default) plots all three.
#' @param scale (`logical(1)`)\cr Only for the plot (`plot = TRUE`): if `TRUE` each
#'   drawn curve is divided by the SD of its process x species rate, yielding
#'   dimensionless, comparable effects (the curve's variance then equals the
#'   Sobol-style normalised importance). The returned tables are not changed; they
#'   carry the SD as `rate_sd`. Default `FALSE`.
#' @param ... Not supported yet.
#'
#' @return A list of class `FINNale` with one table per process (`$growth`,
#'   `$mortality`, `$regeneration`). Each table gives the accumulated local
#'   effect (`ale`, on the scale of the process rate) of every driver (`var`)
#'   across its observed range (`x`), per `species`, and `rate_sd`, the SD of
#'   that process x species rate. When `plot = TRUE` the effects are also drawn.
#' @seealso [plot.FINNale()] to draw the curves, [conditionalEffects()] and
#'   [averageConditionalEffects()] for the local slopes behind them,
#'   [summary.finn_class()] and [feature_importance()] for variable importance
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # a small model fitted on three Oregon FIA sites (?fia_obs_dt, ?fit); the NN
#' # variant is the same call with `hidden = c(10, 10)` on the process
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(
#'   N_species            = Nsp, recruits_dbh = 12.9,
#'   growth_process       = createProcess(~ temp + prec, FINN::growth,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE),
#'   mortality_process    = createProcess(~ temp + prec, FINN::mortality,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE),
#'   regeneration_process = createProcess(~ temp + prec, FINN::regeneration,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06, plot_progress = FALSE)
#'
#' # accumulated local effects of every driver, per process and species. `env` and
#' # `init_cohort` default to the ones fit() cached, so the curves are read over
#' # exactly the range the model was trained on.
#' ale <- ALE(m, plot = FALSE)
#' head(ale$growth)
#'
#' # the curves of one process, one panel per driver and one line per species;
#' # `scale = TRUE` divides each drawn curve by the SD of its process x species
#' # rate, so that the panels of different processes become comparable
#' ALE(m, process = "growth", scale = TRUE)
#' @export
ALE = function(model, env = NULL, init_cohort = NULL, env_autoscale = TRUE, sim_seed = 42L, plot = TRUE, process = NULL, scale = FALSE, ...) {

  ce = .finn_get_ce(model, env, init_cohort, env_autoscale, sim_seed)

  # a NULL key means the effects were computed for inputs the caller passed, which
  # are not cached; two of those would otherwise match on NULL == NULL (#72, B36)
  if(!is.null(ce$key) && !is.null(model$ale) && identical(attr(model$ale, "key"), ce$key)) {
    out = model$ale
  } else {
    out = list()
    for(proc in names(ce$processes)) {
      sp_list    = ce$processes[[proc]]
      df_results = data.frame()
      for(sp in names(sp_list)) {
        s = sp_list[[sp]]
        if(is.null(s) || nrow(s$X) < 2L) next
        ales = ALE_ce(s$X, s$grads, predictions = s$pred)
        ales$species = as.integer(sp)
        rate_sd = stats::sd(s$pred, na.rm = TRUE)
        ales$rate_sd = if(is.finite(rate_sd) && rate_sd > 0) rate_sd else 1
        df_results = rbind(df_results, ales)
      }
      out[[proc]] = df_results
    }
    class(out) = c("FINNale")
    attr(out, "key") = ce$key
    attr(out, "env_vars") = unique(unlist(lapply(ce$processes, function(sp) if(length(sp)) sp[[1]]$env_vars else NULL)))
    if(!is.null(ce$key)) model$ale = out
  }

  # the species names the model keeps, for the legend
  attr(out, "species_names") = species_names_lookup(model$N_species, model = model)
  if(isTRUE(plot)) tryCatch(plot(out, process = process, scale = scale), error = function(e) warning("ALE plot could not be drawn: ", conditionMessage(e), call. = FALSE))
  invisible(out)
}


#' Plot ALE curves of a FINN model
#'
#' @description
#' Plots the accumulated local effect (ALE) curves produced by [ALE()] as a
#' grid: one row per demographic process (growth, mortality, regeneration) and
#' one column per environmental predictor. When several species are present,
#' their curves are overlaid in the same panel as one coloured line each.
#'
#' @param x an object of class `FINNale` (returned by [ALE()]).
#' @param process (`character`|`NULL`)\cr restrict the plot to one process (`"growth"`,
#'   `"mortality"` or `"regeneration"`); `NULL` (default) plots all three.
#' @param scale (`logical(1)`)\cr if `TRUE`, divide each curve by the SD of its
#'   process x species rate so the effects are dimensionless and comparable (a shared
#'   y-axis is then used); if `FALSE` (default) the raw ALE is shown with free y-axes.
#' @param species_names Species names for the legend (see [plot.finn_class()]);
#'   `NULL` (default): the names the model keeps from [fit()].
#' @param ... currently ignored.
#' @return invisibly, the `ggplot` object.
#'
#' @seealso [ALE()]
#'
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # a small model fitted on two Oregon FIA sites (see ?fit), then its ALE curves
#' sites <- 1:2
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(N_species = Nsp, recruits_dbh = 12.9,
#'           growth_process       = createProcess(~ temp + prec, growth),
#'           mortality_process    = createProcess(~ temp + prec, mortality),
#'           regeneration_process = createProcess(~ temp + prec, regeneration_capped))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 1L, patches = 4, patch_size = 0.06,
#'     plot_progress = FALSE)
#'
#' ale <- ALE(m, plot = FALSE)
#' plot(ale, process = "growth")
#' # all three processes on a common, dimensionless scale
#' plot(ale, scale = TRUE)
#'
#' @method plot FINNale
#' @export
plot.FINNale = function(x, process = NULL, scale = FALSE, species_names = NULL, ...) {
  env_vars = attr(x, "env_vars")
  procs    = if(is.null(process)) names(x) else intersect(names(x), process)

  parts = lapply(procs, function(proc) {
    d = x[[proc]]
    if(is.null(d) || !nrow(d)) return(NULL)
    d$process = proc
    d
  })
  df = do.call(rbind, parts)
  if(is.null(df) || !nrow(df)) { message("No ALE curves to plot."); return(invisible(NULL)) }

  if(!is.null(env_vars)) df = df[df$var %in% env_vars, , drop = FALSE]
  if(!nrow(df)) { message("No environmental predictors to plot."); return(invisible(NULL)) }

  df$process = factor(df$process, levels = intersect(c("growth", "mortality", "regeneration"), unique(df$process)))
  labs = if (!is.null(species_names) || is.null(attr(x, "species_names")))
    species_names_lookup(max(df$species), species_names) else attr(x, "species_names")
  sp_levels = sort(unique(df$species))
  df$species = factor(labs[as.character(df$species)], levels = labs[as.character(sp_levels)])

  ylab         = "ALE"
  facet_scales = "free"
  if(isTRUE(scale)) {
    if(!"rate_sd" %in% names(df)) {
      warning("no rate_sd stored on this object; re-run ALE() to scale.", call. = FALSE)
    } else {
      dt = data.table::as.data.table(df)
      dt[, ale := (ale - mean(ale)) / rate_sd, by = list(process, var, species)]
      df           = as.data.frame(dt)
      ylab         = "scaled ALE (rate SD)"
      facet_scales = "free_x"
    }
  }

  facets = if(length(levels(df$process)) == 1L) {
    ggplot2::facet_wrap(~ var, scales = facet_scales)
  } else {
    ggplot2::facet_grid(process ~ var, scales = facet_scales)
  }

  p = ggplot2::ggplot(df, ggplot2::aes(x = x, y = ale, colour = species, group = species)) +
    ggplot2::geom_rug(data = df, mapping = ggplot2::aes(x = x), sides = "b", alpha = 0.15, inherit.aes = FALSE) +
    ggplot2::geom_line() +
    facets +
    ggplot2::labs(x = "environmental predictor value", y = ylab, colour = "species") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.line = ggplot2::element_line())

  print(p)
  invisible(p)
}


#' Conditional effects of a FINN model
#'
#' @description
#' Computes (and caches on the model) the conditional effects — the local
#' derivatives of each demographic process (growth, mortality, regeneration)
#' with respect to its inputs, per process and per species. This is the shared
#' primitive behind [ALE()], [averageConditionalEffects()] and the analytical
#' variable importance in [summary.finn_class()].
#'
#' @param model (`finn_class`)\cr fitted model.
#' @param env (`data.table|data.frame`|`NULL`)\cr env covariates; `NULL` uses the cached training env.
#' @param init_cohort (`CohortMat`|`NULL`)\cr initial cohort; `NULL` uses the cached training init_cohort.
#' @param env_autoscale (`logical(1)`)\cr see [ALE()].
#' @param sim_seed (`integer(1)`)\cr seed for the state-harvesting simulation (see [ALE()]).
#' @return an object of class `FINNconditionalEffects` (also cached on `model$conditional_effects`):
#'   a list with one element per process (`growth`, `mortality`, `regeneration`), each a list with
#'   one element per species holding the inputs the process was evaluated at (`X`), the local
#'   derivatives with respect to them (`grads`) and the predicted rate (`pred`).
#' @seealso [averageConditionalEffects()] for one number per driver, [ALE()] for the curves
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # a small model fitted on three Oregon FIA sites (?fia_obs_dt, ?fit); the NN
#' # variant is the same call with `hidden = c(10, 10)` on the process
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(
#'   N_species            = Nsp, recruits_dbh = 12.9,
#'   growth_process       = createProcess(~ temp + prec, FINN::growth,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE),
#'   mortality_process    = createProcess(~ temp + prec, FINN::mortality,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE),
#'   regeneration_process = createProcess(~ temp + prec, FINN::regeneration,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06, plot_progress = FALSE)
#'
#' # the local slope of each process with respect to each of its inputs, at the
#' # states the fitted model passes through: a list per process and species
#' ce <- conditionalEffects(m)
#' names(ce)
#' # species 1, growth: the inputs, the slopes and the predicted rate
#' str(ce$growth[[1]], max.level = 1)
#' @export
conditionalEffects = function(model, env = NULL, init_cohort = NULL, env_autoscale = TRUE, sim_seed = 42L) {
  .finn_get_ce(model, env, init_cohort, env_autoscale, sim_seed)
}


#' Average conditional effects of a FINN model
#'
#' @description
#' Summarises the conditional effects into a per process x species x variable
#' average marginal effect: the mean local derivative (`mean_effect`, an
#' approximate linear effect). Derived cheaply from the cached conditional
#' effects.
#'
#' @inheritParams conditionalEffects
#' @param env_only (`logical(1)`)\cr if `TRUE` (default) report only the environmental
#'   predictors of each process; if `FALSE` also include the stand-state inputs
#'   (dbh, light, growth) each process depends on.
#' @return a named list (one entry per process) of data.frames with columns
#'   `species, variable, mean_effect`.
#' @seealso [conditionalEffects()], [summary.finn_class()], which prints these effects
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # a small model fitted on three Oregon FIA sites (?fia_obs_dt, ?fit); the NN
#' # variant is the same call with `hidden = c(10, 10)` on the process
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(
#'   N_species            = Nsp, recruits_dbh = 12.9,
#'   growth_process       = createProcess(~ temp + prec, FINN::growth,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE),
#'   mortality_process    = createProcess(~ temp + prec, FINN::mortality,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE),
#'   regeneration_process = createProcess(~ temp + prec, FINN::regeneration,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06, plot_progress = FALSE)
#'
#' # one number per driver and process: the average slope over the observed range,
#' # which is the linear-model coefficient a mechanistic process would report and
#' # the closest thing a neural-network process has to one
#' ace <- averageConditionalEffects(m)
#' head(ace$growth)
#' # with the stand-state inputs of each process (dbh, light, ...) as well
#' head(averageConditionalEffects(m, env_only = FALSE)$mortality)
#' @export
averageConditionalEffects = function(model, env = NULL, init_cohort = NULL, env_autoscale = TRUE, sim_seed = 42L, env_only = TRUE) {
  ce  = .finn_get_ce(model, env, init_cohort, env_autoscale, sim_seed)
  res = list()
  for(process in names(ce$processes)) {
    sp_list = ce$processes[[process]]
    rows = list()
    for(sp in names(sp_list)) {
      s = sp_list[[sp]]
      if(is.null(s) || nrow(s$X) < 1L) next
      cols = colnames(s$grads)
      keep = if(isTRUE(env_only)) intersect(cols, s$env_vars) else cols
      for(v in keep) {
        g = s$grads[, v]
        rows[[length(rows) + 1L]] = data.frame(species = as.integer(sp), variable = v, mean_effect = mean(g, na.rm = TRUE))
      }
    }
    res[[process]] = if(length(rows)) do.call(rbind, rows) else data.frame()
  }
  class(res) = "FINNconditionalEffectsSummary"
  res
}


.importance_ale_var = function(model, env = NULL, init_cohort = NULL, env_autoscale = TRUE, sim_seed = 42L, env_only = TRUE, scale = TRUE) {
  ce  = .finn_get_ce(model, env, init_cohort, env_autoscale, sim_seed)
  ale = ALE(model, env, init_cohort, env_autoscale, sim_seed, plot = FALSE)
  res = list()
  for(process in names(ale)) {
    dfp = ale[[process]]
    if(is.null(dfp) || !nrow(dfp)) { res[[process]] = data.frame(); next }
    env_vars = unique(unlist(lapply(ce$processes[[process]], function(s) s$env_vars)))
    rate_var = vapply(ce$processes[[process]], function(s) { v = stats::var(s$pred, na.rm = TRUE); if(!is.finite(v) || v == 0) 1 else v }, numeric(1))
    dt  = as.data.table(dfp)
    imp = dt[, .(importance = if(.N > 1L) stats::var(ale) else 0), by = list(species, var)]
    if(isTRUE(scale)) imp[, importance := importance / rate_var[as.character(species)]]
    imp[!is.finite(importance), importance := 0]
    if(isTRUE(env_only)) imp = imp[var %in% env_vars]
    data.table::setnames(imp, "var", "variable")
    res[[process]] = as.data.frame(imp[order(species, -importance)])
  }
  res
}


#' Permutation feature importance for FINN demographic rates
#'
#' @description
#' Permutes each environmental predictor of a process and measures how strongly
#' the permutation shifts that process's predicted rate, **per process and per
#' species**. Unlike the analytical ALE-variance importance, this re-simulates
#' the model, so it captures the full dynamical response (feedback through the
#' stand state) rather than the process response function alone. It does NOT use
#' the conditional-effects cache.
#'
#' @details
#' Two scorings via `method`:
#'  * `"rmse"` — RMSE between the unpermuted and permuted rate, in units of that
#'    species' rate SD. Unbounded; larger = more important.
#'  * `"sobol"` — total-effect estimator `0.5 * mean(MSE_shift) / Var(rate)`;
#'    dimensionless (only bounded in `[0, 1]` under independent predictors — FINN's
#'    climate predictors are usually correlated, so treat it as relative).
#'
#' Predictors are read per-process from the model formulas, so processes with
#' different formulas get different variable sets. Common random numbers
#' (`sim_seed`, applied to the torch simulation RNG) make the stochastic
#' mortality/regeneration draws shared across the reference and permuted runs,
#' so a driver with no effect returns ~0.
#'
#' @param model (`finn_class`)\cr fitted model.
#' @param env (`data.table|data.frame`|`NULL`)\cr env covariates; `NULL` uses the cached training env.
#' @param init_cohort (`CohortMat`|`NULL`)\cr init cohort; `NULL` uses the cached training init_cohort.
#' @param nperm (`integer(1)`)\cr number of permutation replicates (default 20).
#' @param method (`character(1)`)\cr `"rmse"` or `"sobol"`.
#' @param seed (`integer(1)`|`NULL`)\cr R RNG seed controlling which permutations are drawn.
#' @param sim_seed (`integer(1)`|`NULL`)\cr torch seed for common random numbers across runs.
#' @param env_autoscale (`logical(1)`)\cr see [ALE()]. `TRUE` (default) leaves
#'   [`predict()`][predict.finn_class] to
#'   rescale raw `env` internally; `FALSE` treats `env` as already on the model scale.
#' @param ... passed to [`predict()`][predict.finn_class] (e.g. `patches`, `patch_size`).
#' @return a named list (one per process) of data.frames with columns
#'   `species, variable, importance`, sorted within species.
#' @seealso [summary.finn_class()] (`importance = "permutation"`), [ALE()]
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # a small model fitted on three Oregon FIA sites (?fia_obs_dt, ?fit); the NN
#' # variant is the same call with `hidden = c(10, 10)` on the process
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(
#'   N_species            = Nsp, recruits_dbh = 12.9,
#'   growth_process       = createProcess(~ temp + prec, FINN::growth,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE),
#'   mortality_process    = createProcess(~ temp + prec, FINN::mortality,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE),
#'   regeneration_process = createProcess(~ temp + prec, FINN::regeneration,
#'                                        optimizeSpecies = TRUE, optimizeEnv = TRUE))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06, plot_progress = FALSE)
#'
#' # permutation importance: how much the simulated rate of a process shifts
#' # when a driver is shuffled between sites and years. `nperm` is small here
#' # to keep the example quick; use more for a real reading.
#' fi <- feature_importance(m, nperm = 2L, seed = 1L)
#' head(fi$growth)
#' # `method = "sobol"` gives the dimensionless Sobol-style score instead
#' @export
feature_importance = function(model, env = NULL, init_cohort = NULL, nperm = 20L, method = c("rmse", "sobol"), seed = NULL, sim_seed = 42L, env_autoscale = TRUE, ...) {

  method = match.arg(method)
  if(!is.null(seed)) set.seed(seed)
  if(is.null(env))         env         = model$train_env
  if(is.null(env))         stop("No `env` supplied and no training env is cached. Pass `env=` or fit the model first.")
  if(is.null(init_cohort)) init_cohort = .finn_train_ic(model)

  # Permute over distinct rows: exact copies of a site-year (the FIA env tables held
  # every row four times up to 0.2.0) would receive different permuted values, of
  # which the environment array used to keep one at random; climateDF2array()
  # now rejects any duplicate (#35).
  env      = unique(data.table::as.data.table(data.table::copy(env)))
  env_vars = setdiff(names(env), c("siteID", "year"))

  scaling = model$env_scaling
  if(!isTRUE(env_autoscale)) model$env_scaling = NULL
  on.exit(model$env_scaling <- scaling, add = TRUE)

  procs = list(
    growth       = list(formula = model$growth_formula,       col = "growth"),
    mortality    = list(formula = model$mortality_formula,    col = "mort"),
    regeneration = list(formula = model$regeneration_formula, col = "reg")
  )
  get_predictors = function(formula, env_vars) {
    if(is.null(formula)) return(character(0))
    if("." %in% all.names(formula)) return(env_vars)
    intersect(all.vars(formula), env_vars)
  }
  crn = function() if(!is.null(sim_seed)) torch::torch_manual_seed(sim_seed)

  crn()
  ref = data.table::as.data.table(predict(model, env, init_cohort = init_cohort, ...)$wide$site)
  species_levels = sort(unique(ref$species))

  out = list()
  for(pname in names(procs)) {
    rate_col = procs[[pname]]$col
    preds    = get_predictors(procs[[pname]]$formula, env_vars)
    if(length(preds) == 0L) { out[[pname]] = data.frame(); next }

    rate_var_sp = vapply(species_levels, function(sp) {
      v = stats::var(ref[[rate_col]][ref$species == sp], na.rm = TRUE)
      if(!is.finite(v) || v == 0) 1 else v
    }, numeric(1))
    names(rate_var_sp) = as.character(species_levels)

    rows = list()
    for(v in preds) {
      mse = matrix(0, nrow = length(species_levels), ncol = nperm, dimnames = list(as.character(species_levels), NULL))
      for(k in seq_len(nperm)) {
        tmp = data.table::copy(env)
        tmp[[v]] = sample(tmp[[v]])
        crn()
        p = data.table::as.data.table(predict(model, tmp, init_cohort = init_cohort, ...)$wide$site)
        d = (ref[[rate_col]] - p[[rate_col]])^2
        agg = tapply(d, ref$species, function(x) mean(x, na.rm = TRUE))
        mse[names(agg), k] = as.numeric(agg)
      }
      for(sp in as.character(species_levels)) {
        m_sp = mse[sp, ]
        imp  = switch(method,
                      rmse  = mean(sqrt(m_sp)) / sqrt(rate_var_sp[[sp]]),
                      sobol = 0.5 * mean(m_sp) / rate_var_sp[[sp]])
        rows[[length(rows) + 1L]] = data.frame(species = as.integer(sp), variable = v, importance = imp)
      }
    }
    df = do.call(rbind, rows)
    out[[pname]] = df[order(df$species, -df$importance), ]
  }
  class(out) = "FINNfeatureImportance"
  out
}


#' Summarise a fitted FINN model
#'
#' @description
#' Collects what is needed to judge a fitted model, in this order:
#'
#' 1. **Skill**: per response, R2 (`1 - SSE/SST`) first, then RMSE, Pearson r,
#'    the bias, the ratio of the mean prediction to the mean observation (is the
#'    magnitude right?; for mortality weighted by the trees at risk, where the
#'    observations carry `n_at_risk` and `n_died`) and n; for mortality also the
#'    tree-level AUC and the ratio of predicted to observed deaths. On the
#'    training data and, when [evaluate()] or the argument `holdout` of [fit()]
#'    stored one, on the holdout. The observations are compared with the
#'    simulation aggregated the way [fit()] scores it (see [comparePlot()]).
#'    Without a stored evaluation, `summary()` runs [evaluate()] on the training
#'    data first (one simulation) and stores it on the model.
#' 2. **Convergence**: by how many percent each loss term changed over the last
#'    `epochs_back` epochs (the mean of the last `epochs_back` epochs against
#'    the mean of the `epochs_back` before), and by how many percent each
#'    parameter moved over the last `epochs_back` epochs (relative L2 change
#'    over the species, `||theta_now - theta_then|| / ||theta_then||`): every
#'    species-specific parameter of a process by its name and on the scale it
#'    acts on (e.g. "mortality: growth"), every network and other parameter
#'    block as a whole. The print names the largest change of each. [fit()]
#'    keeps the parameters of the last 51 epochs for this, so `epochs_back` can
#'    be at most 50; for a model fitted before this version the parameter
#'    snapshots of the fit (`model$param_history`, every `checkpoints` epochs)
#'    are used instead. A fit has settled when both are small; network weights
#'    can keep moving while the loss stays flat.
#' 3. **Bare ground**, when [evaluate()] stored a bare-ground run: per site,
#'    the total basal area and the leading species with its share, as the mean
#'    over the last tenth of the run.
#' 4. **Environmental importance and average conditional effects**, per process
#'    and species, from the model's conditional effects (computed once and
#'    cached, see [ALE()]). Where they cannot be computed for a model, a
#'    warning says why and the other sections are still returned.
#'
#' `print()` shows it; `print(s, epochs_back = 10)` redoes the convergence
#' check over another window without recomputing anything. [plot.finn_class()]
#' draws the fitted parameters themselves.
#'
#' @param object (`finn_class`)\cr fitted model.
#' @param env,init_cohort (`NULL` or data)\cr interpretation data for the
#'   importance and the conditional effects; `NULL` uses the cached training
#'   data (see [ALE()]).
#' @param importance (`character(1)`)\cr `"ale"` (default) = analytical ALE-variance importance
#'   (cheap, from the cache); `"permutation"` = re-simulating permutation importance
#'   (see [feature_importance()]), cached on the model under the same input key.
#' @param env_autoscale (`logical(1)`)\cr see [ALE()].
#' @param sim_seed (`integer(1)`)\cr seed for the simulation (see [ALE()]).
#' @param nperm (`integer(1)`)\cr replicates for `importance = "permutation"`.
#' @param scale (`logical(1)`)\cr for `importance = "ale"`, divide `Var(ALE)` by the
#'   process x species rate variance so the importances are dimensionless and comparable
#'   across processes and species (Sobol-style). Default `TRUE`.
#' @param epochs_back (`integer(1)`)\cr the window of the convergence check, in
#'   epochs. Default `50`.
#' @param ... `summary()`: passed to [feature_importance()] with
#'   `importance = "permutation"` (and from there to
#'   [`predict()`][predict.finn_class]). `print()`: unused.
#' @return `summary()`: an object of class `summary.finn_class`, a list with
#'   `skill` (per set, overall), `skill_species` (per set and species),
#'   `convergence` (the loss history and the recent parameters), `bare_ground`,
#'   `importance`, `average_conditional_effects`, `method` (the `importance`
#'   option), `importance_label`, `environment_error` (`NULL`, or why item 4
#'   could not be computed) and `epochs_back`. `print()`: `x`, invisibly.
#'
#' @seealso [evaluate()], [comparePlot()], [plot.finn_class()], [ALE()],
#'   [averageConditionalEffects()], [feature_importance()]
#'
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # a small model fitted on two Oregon FIA sites (see ?fit)
#' sites <- 1:2
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' FINN.seed(1)
#' m <- finn(N_species = Nsp, recruits_dbh = 12.9,
#'           growth_process       = createProcess(~ temp + prec, growth),
#'           mortality_process    = createProcess(~ temp + prec, mortality),
#'           regeneration_process = createProcess(~ temp + prec, regeneration_capped))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 3L, patches = 4, patch_size = 0.06,
#'     plot_progress = FALSE)
#'
#' s <- summary(m)
#' s
#' # the convergence check over the last two epochs instead
#' print(s, epochs_back = 2)
#' # the tables are also returned
#' s$skill$train
#'
#' @method summary finn_class
#' @export
summary.finn_class = function(object, env = NULL, init_cohort = NULL, importance = c("ale", "permutation"), env_autoscale = TRUE, sim_seed = 42L, nperm = 20L, scale = TRUE, epochs_back = 50L, ...) {
  importance = match.arg(importance)

  # 1. skill, from the stored evaluation (run on the training data if missing)
  if (is.null(object$evaluation) && !is.null(object$train_data)) {
    message("Evaluating on the training data (one simulation); evaluate() can add a holdout and a bare-ground run.")
    # a failed comparison must not cost the rest of the summary
    tryCatch(evaluate(object), error = function(e)
      warning(sprintf("summary(): the model could not be evaluated on its training data (%s); the skill section is left out.",
                      conditionMessage(e)), call. = FALSE))
  }
  ev = object$evaluation
  sets = intersect(c("train", "holdout"), names(ev))
  skill = stats::setNames(lapply(sets, function(s) ev[[s]]$skill), sets)
  skill_species = stats::setNames(lapply(sets, function(s) ev[[s]]$skill_species), sets)

  # 2. convergence: the stored history, judged at print time
  loss = if (length(object$history)) do.call(rbind, object$history) else NULL
  # the mechanistic parameters by name, mapped to the scale they act on
  # (lower + (upper - lower) * plogis(x)), so their change reads in their own units
  mech = list()
  for (type in c("growth", "mortality", "regeneration", "competition")) {
    lo = object[[paste0("par_", type, "_lower")]]; up = object[[paste0("par_", type, "_upper")]]
    if (is.null(lo) || is.null(up)) next
    lo = as.numeric(torch::as_array(lo$cpu())); up = as.numeric(torch::as_array(up$cpu()))
    mech[[paste0("par_", type, "_unconstrained")]] =
      list(type = type, lower = lo, upper = up, labels = parameter_labels(object, type, length(lo)))
  }
  # snapshots every `checkpoints` epochs (model$param_history): the fallback for a
  # model that keeps no per-epoch parameters (fitted before this version)
  snap_every = if (length(object$param_history) && length(object$history))
    max(1L, as.integer(round(length(object$history) / length(object$param_history)))) else NA_integer_
  convergence = list(loss = loss, params = object$param_recent, epochs = object$param_recent_epochs,
                     snapshots = object$param_history, snap_every = snap_every, mechanistic = mech)

  # 3. bare ground
  bare = NULL
  if (!is.null(ev$bare_ground)) {
    d = data.table::as.data.table(ev$bare_ground$long$site)[variable == "ba"]
    last = max(d$year)
    end = d[year > last - max(1L, floor(0.1 * last)), list(ba = mean(value)), by = list(siteID, species)]
    sp_lab = species_names_lookup(object$N_species, model = object)
    bare = end[, list(total_ba = sum(ba), leading_species = sp_lab[as.character(species[which.max(ba)])],
                      leading_share = max(ba) / sum(ba)), by = siteID]
    bare[, years := last]
  }

  # 4. importance and average conditional effects; if they cannot be computed
  # for this model, the sections above are still returned and printed
  env_section = tryCatch({
  if(importance == "ale") {
    vi       = .importance_ale_var(object, env, init_cohort, env_autoscale, sim_seed, scale = scale)
    vi_label = if(isTRUE(scale)) "Analytical ALE importance (rate-normalised)" else "Analytical ALE-variance importance"
  } else {
    env_res = if(is.null(env)) object$train_env else env
    ic_res  = if(is.null(init_cohort)) .finn_train_ic(object) else init_cohort
    key = .finn_ce_key(data.table::as.data.table(env_res), ic_res, sim_seed, env_autoscale, object$fit_id)
    if(!is.null(object$perm_importance) && identical(object$perm_importance$key, key)) {
      vi = object$perm_importance$result
    } else {
      vi = feature_importance(object, env, init_cohort, nperm = nperm, sim_seed = sim_seed, env_autoscale = env_autoscale, ...)
      object$perm_importance = list(key = key, result = vi)
    }
    vi_label = sprintf("Permutation importance (nperm = %d)", nperm)
  }
  ace = averageConditionalEffects(object, env, init_cohort, env_autoscale, sim_seed)
  list(vi = vi, ace = ace, vi_label = vi_label, error = NULL)
  }, error = function(e) {
    warning(sprintf("summary(): the environmental importance and conditional effects could not be computed for this model (%s); the other sections are shown.",
                    conditionMessage(e)), call. = FALSE)
    list(vi = NULL, ace = NULL, vi_label = NULL, error = conditionMessage(e))
  })

  structure(list(skill = skill, skill_species = skill_species, convergence = convergence,
                 bare_ground = bare, importance = env_section$vi,
                 average_conditional_effects = env_section$ace,
                 method = importance, importance_label = env_section$vi_label,
                 species_names = species_names_lookup(object$N_species, model = object),
                 environment_error = env_section$error, epochs_back = as.integer(epochs_back)),
            class = "summary.finn_class")
}

# Percent change of every loss term over the last k epochs, and of every
# parameter block over the last k (at most 50) stored epochs.
convergence_table = function(conv, k) {
  k = as.integer(k)
  out = list(loss = NULL, params = NULL, k_loss = NA_integer_, k_params = NA_integer_)
  L = conv$loss
  if (!is.null(L) && nrow(L) > k) {
    n = nrow(L)
    now  = colMeans(L[(n - k + 1):n, , drop = FALSE], na.rm = TRUE)
    then = colMeans(L[max(1, n - 2 * k + 1):(n - k), , drop = FALSE], na.rm = TRUE)
    terms = c("dbh", "ba", "trees", "growth", "mortality", "regeneration", "prior penalty")
    nm = if (ncol(L) == length(terms)) terms else paste0("term ", seq_len(ncol(L)))
    keep = is.finite(now) & is.finite(then) & then != 0
    out$loss = data.frame(term = nm[keep], loss_now = now[keep],
                          change_pct = 100 * (now[keep] - then[keep]) / abs(then[keep]))
    out$k_loss = k
  }
  # the per-epoch record of the last 51 epochs, else the snapshots of the fit
  P = conv$params; step = 1L; src = "epoch"
  if (length(P) < 2L && length(conv$snapshots) >= 2L && is.finite(conv$snap_every)) {
    P = conv$snapshots; step = conv$snap_every; src = "snapshot"
  }
  if (length(P) >= 2L) {
    kk = min(max(1L, as.integer(round(k / step))), length(P) - 1L)
    a = P[[length(P)]]; b = P[[length(P) - kk]]
    rel = function(x, y) {
      den = sqrt(sum(y^2))
      if (!is.finite(den) || den == 0) return(NA_real_)
      100 * sqrt(sum((x - y)^2)) / den
    }
    rows = list()
    for (nm in intersect(names(a), names(b))) {
      mk = conv$mechanistic[[nm]]
      if (!is.null(mk) && is.matrix(a[[nm]]) && ncol(a[[nm]]) == length(mk$lower)) {
        # one row per named parameter, on the scale it acts on (over all species)
        for (j in seq_len(ncol(a[[nm]]))) {
          on_scale = function(z) mk$lower[j] + (mk$upper[j] - mk$lower[j]) * stats::plogis(z[, j])
          rows[[length(rows) + 1L]] = data.frame(parameter = sprintf("%s: %s", mk$type, mk$labels[j]),
                                                 kind = "species parameter",
                                                 change_pct = rel(on_scale(a[[nm]]), on_scale(b[[nm]])))
        }
      } else {
        rows[[length(rows) + 1L]] = data.frame(parameter = nm,
                                               kind = if (grepl("^nn_", nm)) "network weights" else "other",
                                               change_pct = rel(a[[nm]], b[[nm]]))
      }
    }
    pr = do.call(rbind, rows)
    out$params = pr[order(-pr$change_pct), , drop = FALSE]
    out$k_params = kk * step
    out$params_source = src
  }
  out
}

#' @rdname summary.finn_class
#' @param x A `summary.finn_class` object.
#' @param species (`logical(1)`)\cr also print the skill per species. Default `FALSE`.
#' @export
print.summary.finn_class = function(x, epochs_back = x$epochs_back, species = FALSE, ...) {
  cat("FINN model summary\n==================\n")

  fmt_skill = function(s) {
    d = as.data.frame(s)
    cols = intersect(c("variable", "species", "r2", "rmse", "r", "bias", "ratio", "n", "auc", "deaths_ratio"), names(d))
    d = d[, cols, drop = FALSE]
    names(d) = c(variable = "response", species = "species", r2 = "R2", rmse = "RMSE", r = "r", bias = "bias",
                 ratio = "pred/obs", n = "n", auc = "AUC", deaths_ratio = "deaths pred/obs")[cols]
    d
  }
  cat("\n## Skill (observed vs simulated, aggregated as fit() scores)\n")
  if (!length(x$skill)) {
    cat("  (no evaluation: run evaluate() on a model fitted with this version)\n")
  } else for (s in names(x$skill)) {
    cat(sprintf("\n%s:\n", if (s == "train") "Training data" else "Holdout"))
    print(fmt_skill(x$skill[[s]]), row.names = FALSE, digits = 3)
    if (isTRUE(species)) {
      cat(sprintf("\n%s, per species:\n", if (s == "train") "Training data" else "Holdout"))
      ss = data.table::copy(x$skill_species[[s]])
      if (!is.null(x$species_names)) ss[, species := x$species_names[as.character(species)]]
      print(fmt_skill(ss), row.names = FALSE, digits = 3)
    }
  }

  cv = convergence_table(x$convergence, epochs_back)
  cat("\n## Convergence\n")
  pct = function(v) sprintf("%.2f %%", v)
  if (is.null(cv$loss)) {
    cat(sprintf("  loss: fewer than %d + 1 epochs recorded\n", epochs_back))
  } else {
    l = cv$loss; i = which.max(abs(l$change_pct))
    cat(sprintf("\nLoss: the largest change is %s, %s over the last %d epochs (mean of the last %d epochs against the %d before).\n",
                l$term[i], pct(l$change_pct[i]), cv$k_loss, cv$k_loss, cv$k_loss))
    l$change = pct(l$change_pct); l$change_pct = NULL
    print(l, row.names = FALSE, digits = 3)
  }
  if (is.null(cv$params)) {
    cat("  parameters: none stored (a model fitted before this version without parameter snapshots, or a single epoch)\n")
  } else {
    p = cv$params
    sp = p[p$kind == "species parameter", , drop = FALSE]
    cat(sprintf("\nParameters, change over the last %d epochs (relative L2 over the species, %%; from the %s):\n",
                cv$k_params, if (identical(cv$params_source, "snapshot")) "parameter snapshots of the fit" else "parameters of every epoch"))
    if (nrow(sp)) cat(sprintf("  the species parameter that moves most: %s, %s\n", sp$parameter[1], pct(sp$change_pct[1])))
    cat(sprintf("  the block that moves most overall: %s, %s\n", p$parameter[1], pct(p$change_pct[1])))
    p$change = pct(p$change_pct); p$change_pct = NULL
    print(utils::head(p, 12L), row.names = FALSE)
    if (nrow(p) > 12L) cat(sprintf("  ... %d more in summary(m)$convergence\n", nrow(p) - 12L))
  }

  if (!is.null(x$bare_ground)) {
    cat("\n## Bare ground (mean of the last tenth of the run)\n")
    print(as.data.frame(x$bare_ground), row.names = FALSE, digits = 3)
  }

  wide = function(df, value_col) {
    if(is.null(df) || !nrow(df)) return(df)
    dt = data.table::as.data.table(df)
    w  = data.table::dcast(dt, variable ~ species, value.var = value_col)
    sp = setdiff(names(w), "variable")
    # species by name where the model knows them, else sp1, sp2, ...
    nm = if (!is.null(x$species_names) && !all(grepl("^species [0-9]+$", x$species_names[sp])))
      x$species_names[sp] else paste0("sp", sp)
    data.table::setnames(w, sp, nm)
    spc = nm
    w   = w[order(-rowMeans(abs(as.matrix(w[, ..spc])), na.rm = TRUE))]
    as.data.frame(w)
  }
  cat("\n## Environment\n")
  if (!is.null(x$environment_error)) {
    cat(sprintf("  (not available for this model: %s)\n", x$environment_error))
    return(invisible(x))
  }
  for(process in c("growth", "mortality", "regeneration")) {
    cat(sprintf("\n### %s\n", toupper(process)))
    cat(sprintf("\n%s (species in columns):\n", x$importance_label))
    vp = x$importance[[process]]
    if(is.null(vp) || !nrow(vp)) { cat("  (no environmental predictors)\n") } else { print(wide(vp, "importance"), row.names = FALSE, digits = 3) }
    cat("\nAverage conditional effects (mean; species in columns):\n")
    ap = x$average_conditional_effects[[process]]
    if(is.null(ap) || !nrow(ap)) { cat("  (no environmental predictors)\n") } else { print(wide(ap, "mean_effect"), row.names = FALSE, digits = 3) }
  }
  invisible(x)
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

