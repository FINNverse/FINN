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

  years      = env_dt$year   |> unique() |> length()
  sitesTotal = env_dt$siteID |> unique() |> length()
  patches    = dim(model$raw_g[[1]])[2]

  processes = list()

  for(process in c("mortality", "growth", "regeneration")) {
    proc_formula = switch(process,
                          growth       = model$process_growth$formula,
                          mortality    = model$process_mortality$formula,
                          regeneration = model$process_regeneration$formula)
    df = data.frame()

    for(timeInd in 1:years) {
      time = (env_dt$year |> unique())[timeInd]
      for(siteInd in 1:sitesTotal) {
        site = (env_dt$siteID |> unique())[siteInd]
        for(patch in 1:patches) {
          env_tmp  = as.matrix(model.matrix(proc_formula, data = env_pred[siteID==site & year == time][1,-(1:2),drop=F]))
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
            pred2 = model$nn_growth(pred2); pred2 = FINN::index_species(pred2, species)
          }
          if(process == "mortality" && !inherits(model$process_mortality, "hybrid")) {
            pred2 = model$nn_mortality(pred2); pred2 = FINN::index_species(pred2, species)
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
  if(is.null(ic)) NULL else ic$cohort
}

.finn_get_ce = function(model, env = NULL, init_cohort = NULL, env_autoscale = TRUE, sim_seed = 42L) {
  if(is.null(env))         env         = model$train_env
  if(is.null(env))         stop("No `env` supplied and no training env is cached. Pass `env=` or fit the model first.")
  if(is.null(init_cohort)) init_cohort = .finn_train_ic(model)

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
#' @param scale (`logical(1)`)\cr If `TRUE` each curve is divided by the SD of its
#'   process x species rate, yielding dimensionless, comparable effects (the curve's
#'   variance then equals the Sobol-style normalised importance). Default `FALSE`.
#' @param ... Not supported yet.
#'
#' @return A list with one table per process (e.g. `$growth`, `$mortality`,
#'   `$regeneration`). Each table gives the accumulated local effect (`ale`) of
#'   every driver (`var`) across its observed range (`x`), per `species`. When
#'   `plot = TRUE` the effects are also drawn.
#' @export
ALE = function(model, env = NULL, init_cohort = NULL, env_autoscale = TRUE, sim_seed = 42L, plot = TRUE, process = NULL, scale = FALSE, ...) {

  ce = .finn_get_ce(model, env, init_cohort, env_autoscale, sim_seed)

  if(!is.null(model$ale) && identical(attr(model$ale, "key"), ce$key)) {
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
    model$ale = out
  }

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
#' @param ... currently ignored.
#' @return invisibly, the `ggplot` object.
#' @method plot FINNale
#' @export
plot.FINNale = function(x, process = NULL, scale = FALSE, ...) {
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
  df$species = factor(df$species)

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
#' @return an object of class `FINNconditionalEffects` (also cached on `model$conditional_effects`).
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
#' @param env_autoscale (`logical(1)`)\cr see [ALE()]. `TRUE` (default) leaves `predict()` to
#'   rescale raw `env` internally; `FALSE` treats `env` as already on the model scale.
#' @param ... passed to [predict()] (e.g. `patches`, `patch_size`).
#' @return a named list (one per process) of data.frames with columns
#'   `species, variable, importance`, sorted within species.
#' @export
feature_importance = function(model, env = NULL, init_cohort = NULL, nperm = 20L, method = c("rmse", "sobol"), seed = NULL, sim_seed = 42L, env_autoscale = TRUE, ...) {

  method = match.arg(method)
  if(!is.null(seed)) set.seed(seed)
  if(is.null(env))         env         = model$train_env
  if(is.null(env))         stop("No `env` supplied and no training env is cached. Pass `env=` or fit the model first.")
  if(is.null(init_cohort)) init_cohort = .finn_train_ic(model)

  env      = data.table::as.data.table(data.table::copy(env))
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
#' Prints, per process and per species, the environmental variable importance
#' and the average conditional effects. Both are derived from the model's
#' conditional effects, which are computed once and cached — so if [ALE()] was
#' already run (or a previous `summary()`), the default path needs no further
#' simulation.
#'
#' @param object (`finn_class`)\cr fitted model.
#' @param env,init_cohort (`NULL` or data)\cr interpretation data; `NULL` uses the cached
#'   training data (see [ALE()]).
#' @param importance (`character(1)`)\cr `"ale"` (default) = analytical ALE-variance importance
#'   (cheap, from the cache); `"permutation"` = re-simulating permutation importance
#'   (see [feature_importance()]), cached on the model under the same input key.
#' @param env_autoscale (`logical(1)`)\cr see [ALE()].
#' @param sim_seed (`integer(1)`)\cr seed for the simulation (see [ALE()]).
#' @param nperm (`integer(1)`)\cr replicates for `importance = "permutation"`.
#' @param scale (`logical(1)`)\cr for `importance = "ale"`, divide `Var(ALE)` by the
#'   process x species rate variance so the importances are dimensionless and comparable
#'   across processes and species (Sobol-style). Default `TRUE`.
#' @param ... passed through (e.g. to [predict()] for the permutation option).
#' @return invisibly, a list with `importance`, `average_conditional_effects`, and `method`.
#' @method summary finn_class
#' @export
summary.finn_class = function(object, env = NULL, init_cohort = NULL, importance = c("ale", "permutation"), env_autoscale = TRUE, sim_seed = 42L, nperm = 20L, scale = TRUE, ...) {
  importance = match.arg(importance)

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

  wide = function(df, value_col) {
    if(is.null(df) || !nrow(df)) return(df)
    dt = data.table::as.data.table(df)
    w  = data.table::dcast(dt, variable ~ species, value.var = value_col)
    sp = setdiff(names(w), "variable")
    data.table::setnames(w, sp, paste0("sp", sp))
    spc = paste0("sp", sp)
    w   = w[order(-rowMeans(abs(as.matrix(w[, ..spc])), na.rm = TRUE))]
    as.data.frame(w)
  }

  cat("FINN model summary\n==================\n")
  for(process in c("growth", "mortality", "regeneration")) {
    cat(sprintf("\n### %s\n", toupper(process)))
    cat(sprintf("\n%s (species in columns):\n", vi_label))
    vp = vi[[process]]
    if(is.null(vp) || !nrow(vp)) { cat("  (no environmental predictors)\n") } else { print(wide(vp, "importance"), row.names = FALSE, digits = 3) }
    cat("\nAverage conditional effects (mean; species in columns):\n")
    ap = ace[[process]]
    if(is.null(ap) || !nrow(ap)) { cat("  (no environmental predictors)\n") } else { print(wide(ap, "mean_effect"), row.names = FALSE, digits = 3) }
  }

  invisible(list(importance = vi, average_conditional_effects = ace, method = importance))
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

