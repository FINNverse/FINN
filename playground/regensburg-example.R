library(FINN)
library(torch)
library(data.table)
library(ggplot2)

# the data read in here is created by the script
# "fit-fia-data/FIA_dataprep_vignette.R"
env_scales_dt <- fread("vignettes/fit-fia-data/env_scales_dt.csv")
tree_dt <- fread("vignettes/fit-fia-data/tree_dt.csv")
env_unscaled_dt <- fread("vignettes/fit-fia-data/env_unscaled_dt.csv")
env_dt <- fread("vignettes/fit-fia-data/env_dt.csv")
env_scales_dt <- fread("vignettes/fit-fia-data/env_scales_dt.csv")
obs_dt <- fread("vignettes/fit-fia-data/obs_dt.csv")
full_obs_dt <- fread("vignettes/fit-fia-data/full_obs_dt.csv")
full_tree_dt <- fread("vignettes/fit-fia-data/full_tree_dt.csv")
species_dt <- fread("vignettes/fit-fia-data/species_dt.csv")


init_trees <- full_tree_dt[year == 0,.(siteID, patchID, year, species, species_name, dbh, trees, living)]
init_trees[living != T | is.na(living), trees := 0,] # set number of trees of all non living trees to 0
summary(init_trees)
init_cohorts <- FINN::makeInitCohorts(init_trees, Nspecies = max(obs_dt$species))

m1 <- finn(
  N_species            = uniqueN(obs_dt$species),
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0, FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createProcess(~., FINN::growth,   optimizeSpecies = TRUE, optimizeEnv = TRUE, hidden = c(20L, 20L)),
  regeneration_process = createProcess(~., FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE, hidden = c(20L, 20L)),
  mortality_process    = createProcess(~., FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)

m1 |> fit(
  env        = env_dt,
  data       = obs_dt,
  init_cohort = init_cohorts,
  device     = "cpu",
  epochs     = 100,
  batchsize  = 500L,
  patch_size = 0.06,
  patches    = 4, weights = c(0.1, 10, 1.0, 10.0, 1, 1),
  lr         = 0.01#, loss = c("mse","mse","mse","mse","mse","mse","mse")
)

predictions =
  m1 |> simulateForest(env = env_dt)

plot(m1, pars = "env", env_names = env_scales_dt$variable, species_names = species_dt$species_name)
plot(m1, pars = "process", species_names = species_dt$species_name)


## Hybrid
init_cohorts <- FINN::makeInitCohorts(init_trees, Nspecies = max(obs_dt$species))

m1 <- finn(
  N_species            = uniqueN(obs_dt$species),
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0, FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createHybrid(~.,  transformer = FALSE, emb_dim = 10L),
  #regeneration_process = createHybrid(~.,  transformer = FALSE, emb_dim = 5L),
  regeneration_process = createProcess(~., FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~., FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)


m1 |> fit(
  env        = env_dt,
  data       = obs_dt,
  init_cohort = init_cohorts,
  device     = "cpu",
  epochs     = 5000L,
  batchsize  = 500L,
  patch_size = 0.06,
  patches    = 4, weights = c(0.1, 10, 1.0, 10.0, 1, 1),
  lr         = 0.01#, loss = c("mse","mse","mse","mse","mse","mse","mse")
)


model = m1


torch::torch_save(m1, path = "test_model.pt")


ALE_ce = function(X, ce, predictions = NULL, center = FALSE) {
  ales =
    lapply(1:ncol(X), function(i) {
      effs = ce[,i]
      data = X[,i]
      ord = order(data)

      x_sorted <- data[ord]
      g_sorted <- effs[ord]

      if(!is.null(predictions)) {
        intercept = mean(predictions[ord])
      } else {
        intercept = 0
      }
      browser()
      ale_vals <- cumsum( (g_sorted[-1] + g_sorted[-length(g_sorted)]) / 2 * diff(x_sorted) )
      if(!center) ale_vals <- ale_vals - mean(ale_vals)
      ale_vals = ale_vals + intercept
      data.frame(x = x_sorted[-1], ale = ale_vals, var = colnames(X)[i])
    })
  return(do.call(rbind, ales))
}


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


assignInNamespace("ALE_ce", ALE_ce, ns = "FINN")
init_cohorts <- FINN::makeInitCohorts(init_trees, Nspecies = max(obs_dt$species))


m1 = torch::torch_load("test_model.pt")
set.seed(42)
torch::torch_manual_seed(42)
ales = FINN::ALE(m1, env_dt, init_cohort = init_cohorts)
set.seed(42)
torch::torch_manual_seed(42)
ales2 = ALE2(m1, env_dt, init_cohort = init_cohorts)



ales$growthinternal[ales$growthinternal$species == 17,][order(ales$growthinternal[ales$growthinternal$species == 17,]$dbh),]

ales$mortality

ggplot(ales$growth, aes(x = x, y = ale, color = as.factor(species))) +
  geom_line() + facet_wrap(~var, scales = "free")

ggplot(ales2$growth, aes(x = x, y = ale, color = as.factor(species))) +
  geom_line() + facet_wrap(~var, scales = "free")




model = m1
env = env_dt
init_cohort = init_cohorts

ALE2 = function(model, env, init_cohort = NULL) {

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

  for(process in c("growth")) {
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

    if(!process == "regeneration") {

      for(i in unique(df$species)) {
        SPECIES = i
        df_tmp = df[df$species == SPECIES,]
        if(nrow(df_tmp) > 1) {

          predict_function = function(model, newdata, tensor = FALSE) {
            df_tmp = newdata
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
            if(!tensor) return(process_output[,1,1] |> as.numeric() )
            else return(process_output)
          }

          expl = DALEX::explain(model, data= df_tmp[,!c("species", "time")],y = runif(nrow(df_tmp)), predict_function =predict_function )
          ale <- DALEX::model_profile(expl,
                                      #variables = "light",
                                      method = "ale")

          ale$agr_profiles
          ales = data.frame(x = ale$agr_profiles$`_x_`, ale = ale$agr_profiles$`_yhat_`, var = ale$agr_profiles$`_vname_`, species = i)
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

          ales = ALE_ce(as.data.frame(df_tmp), grads)
          ales$species = J
          df_results_ale_process = rbind(df_results_ale_process, ales)
        }
      }

    }
    out[[process]] = df_results_ale_process
  }
  class(out) = c("FINNale")
  return(out)
}


#' Calculate accumulated local effect plots from gradients
#'
#' @param X data matrix
#' @param ce gradients at X with respect to the predictions from the model
ALE_ce = function(X, ce) {
  ales =
    lapply(1:ncol(X), function(i) {
      effs = ce[,i]
      data = X[,i]
      ord = order(data)

      x_sorted <- data[ord]
      g_sorted <- effs[ord]

      ale_vals <- cumsum( (g_sorted[-1] + g_sorted[-length(g_sorted)]) / 2 * diff(x_sorted) )
      #ale_vals <- ale_vals - mean(ale_vals)
      data.frame(x = x_sorted[-1], ale = ale_vals, var = colnames(X)[i])
    })
  return(do.call(rbind, ales))
}

assignInNamespace("ALE_ce", ALE_ce, ns = "FINN")
ALE_ce <- function(X, ce, nbins = 20) {
  ales <- lapply(seq_len(ncol(X)), function(i) {
    x  <- X[, i]; g <- ce[, i]
    # quantile bin conditional mean
    brk  <- unique(quantile(x, probs = seq(0,1,length.out = nbins+1), na.rm=TRUE))
    bin  <- cut(x, breaks = brk, include.lowest = TRUE, right = FALSE)
    m    <- tapply(g, bin, function(z) mean(z, na.rm=TRUE))
    mids <- (brk[-1] + brk[-length(brk)])/2
    mvec <- rep(0, length(mids)); mvec[!is.na(m)] <- m[!is.na(m)]

    ale  <- c(0, cumsum( (mvec[-1]+mvec[-length(mvec)])/2 * diff(mids) ))
    #ale  <- ale - mean(ale)  # center!
    data.frame(x = mids, ale = ale, var = colnames(X)[i])
  })
  do.call(rbind, ales)
}

library(cito)

m = dnn(as.factor(survived)~age+fare, data = titanic_imputed, loss = "binomial", lr = 0.001)
ce = m$conditional_effects$response
ceR = apply(ce[[1]]$result, 1, diag) |> t()
ALE_ce(m$data$X, ceR) |> ggplot( aes(x = x, y = ale)) +
  geom_line() + facet_wrap(~var, scales = "free")

expl = DALEX::explain(m, data= data.frame(m$data$X),y = m$data$Y, predict_function =function(model, newdata) predict(model, newdata, type = "response")[,1] )
ale <- DALEX::model_profile(expl,
                            #variables = "light",
                            method = "ale", grid_points = 50)
ales = data.frame(x = ale$agr_profiles$`_x_`, ale = ale$agr_profiles$`_yhat_`, var = ale$agr_profiles$`_vname_`)

ales |>  ggplot( aes(x = x, y = ale)) +
  geom_line() + facet_wrap(~var, scales = "free")


ale$agr_profiles
ales = data.frame(x = ale$agr_profiles$`_x_`, ale = ale$agr_profiles$`_yhat_`, var = ale$agr_profiles$`_vname_`, species = SPECIES)
df_results_ale_process = rbind(df_results_ale_process, ales)



set.seed(1)
n <- 500
x <- sort(runif(n, -1, 2))
f <- 2 * x + 1          # simple linear
g <- rep(2, n)          # derivative wrt x

df <- data.frame(x = x, f = f, g = g)

ALE_test <- ALE_ce(as.matrix(df["x"]), ce = as.matrix(df["g"]))
plot(df$x[-1], f[-1] - f[1], type = "l", col = 1, main = "f(x) vs ALE_ce", ylab = "")
lines(ALE_test$x, ALE_test$ale, col = 2)
legend("topleft", legend = c("f(x) - f(x_min)", "ALE_ce"), col = c(1,2), lty = 1)


df = data.frame(x = x, egal = runif(500),f = f )

m = lm(f~x, data = df)

expl = DALEX::explain(m, data= data.frame(x = x, egal = runif(500)),y = f )
ale <- DALEX::model_profile(expl,
                            variables = "x",
                            method = "ale", grid_points = 50)
plot(ale)
ales = data.frame(x = ale$agr_profiles$`_x_`, ale = ale$agr_profiles$`_yhat_`, var = ale$agr_profiles$`_vname_`)

ALE_test <- ALE_ce(as.matrix(df["x"]), ce = as.matrix(df["g"]))
plot(df$x[-1], f[-1] - f[1], type = "l", col = 1, main = "f(x) vs ALE_ce", ylab = "")
lines(ALE_test$x, ALE_test$ale, col = 2)
lines(ales$x, ales$ale, col = 3)
legend("topleft", legend = c("f(x) - f(x_min)", "ALE_ce"), col = c(1,2), lty = 1)
m = dnn(f~x, data = data.frame(x = x, f = f, egal = runif(500)))
ce = m$conditional_effects$response
ceR = apply(ce[[1]]$result, 1, diag) |> t()
ale_cito = ALE_ce(m$data$X, t(ceR))
lines(ale_cito$x, ale_cito$ale, col = 4, lwd = 2)

library(ale)
res = ale::ALE(m, data = df, max_num_bins = 100L)
lines(res@effect$f$ale$d1$x$x.ceil, res@effect$f$ale$d1$x$.y, )
