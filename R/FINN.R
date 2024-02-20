#' FINN model
#'
#' @param formula formula
#' @param data environmental data, must contain the variables timestep and site
#' @param response species dbh and number of trees, must contain the variables timestep, site, nTree, and Species
#' @param patches number of patches
#' @param epochs number of epochs (optimization)
#' @param learning_rate learning rate for stochastic gradient descent
#' @param batch_size batch_size for stochastic gradient descent
#' @param device cpu or cuda, if available
#' @param start_time when to record, fractional, from 0 to 1
#' @param X alternative parametrization, dims: (site, time, predictors)
#' @param Y alternative parametrization, dims: (site, time, species, 2), 2 for dbh/ba/ba*t and number of Trees (must be integers)
#' @param init_global Initial values for global parameters (optional).
#' @param init_growth Initial values for growth parameters (optional).
#' @param init_mort Initial values for mortality parameters (optional).
#' @param init_reg Initial values for regeneration parameters (optional).
#' @param init_growth_env Initial values for NN env, must be a list of matrices (each matrix == one hidden layer)
#' @param init_mort_env Initial values for NN env, must be a list of matrices (each matrix == one hidden layer)
#' @param init_reg_env Initial values for NN env, must be a list of matrices (each matrix == one hidden layer)
#' @param hidden_growth List of specifications for hidden layers in growth NN (optional).
#' @param hidden_mort List of specifications for hidden layers in mort NN (optional).
#' @param hidden_reg List of specifications for hidden layers in reg NN (optional)
#' @param simulate fit or simulation modus, if simulation, model is not fitted and will be returned
#'
#' @return An object of class `FINN` containing the fitted model, data, and functions to access the model and predictions.
#' @export
#' @author Maximilian Pichler

FINN = function(formula,
                data = NULL,
                response = NULL,
                patches = 10L,
                epochs = 1L,
                learning_rate = 1.0,
                batch_size = 70L,
                device = "cpu",
                start_time = 0.7,
                X = NULL,
                Y = NULL,
                init_global = NULL,
                init_growth = NULL,
                init_mort = NULL,
                init_reg = NULL,
                init_growth_env = NULL,
                init_mort_env = NULL,
                init_reg_env = NULL,
                hidden_growth = list(),
                hidden_mort = list(),
                hidden_reg = list(),
                simulate = FALSE) {
  require(dplyr)

  if(is.null(X)) {

    if(is.data.frame(data)) {
      mf = match.call()
      m = match("formula", names(mf))
      if(inherits(mf[3]$formula, "name")) mf[3]$formula = eval(mf[3]$formula, envir = parent.env(environment()))
      formula = stats::as.formula(mf[m]$formula)
      X = stats::model.matrix(formula, data)
    } else {
      mf = match.call()
      m = match("formula", names(mf))
      if(inherits(mf[3]$formula, "name")) mf[3]$formula = eval(mf[3]$formula, envir = parent.env(environment()))
      formula = stats::as.formula(mf[m]$formula)
      data = data.frame(data)
      X = stats::model.matrix(formula, data)
    }
    XX = as.data.frame(cbind(timestep = data$timestep, site = data$site, X))
    timesteps = length(unique(XX$timestep))
    sites = length(unique(XX$site))

    # sites, time, patch, species
    # slow as fuck...
    obs_array =
      lapply(1:sites, function(site) {
        tmp_site = response[response$site == site, ]
        tmps =
          lapply(1:timesteps, function(time) {
            tmp =
              tmp_site[tmp_site$timestep == time, ] %>%
              group_by(Species) %>%
              summarise(dbh = mean(dbh), nTree = mean(nTree)) %>%
              select(-Species) %>%
              as.data.frame()
            return(tmp)
          })
        return(abind::abind(tmps, along = 0L))
      })

    observations = abind::abind(obs_array, along = 0L)

    env_array =
      lapply(1:sites, function(site) {
        tmp_site = XX[XX$site == site, ]
        tmps =
          lapply(1:timesteps, function(time) {
            tmp =
              tmp_site[tmp_site$timestep == time, ] %>% select(-timestep, -site) %>% c()
            return(matrix(unlist(tmp), length(tmp), 1))
          })
        return(abind::abind(tmps, along = 0L))
      })
    environment = abind::abind(env_array, along = 0L)
    environment = environment[,,,1]


    ### Transform outputs
    torch = pkg.env$torch
    # basal aera * number of trees, for now the default!
    ba_N = pkg.env$FINN$BA_T_P(torch$tensor(observations[,,,1]), torch$tensor(observations[,,,2]))$unsqueeze(3L)$cpu()$data$numpy() |> force_r()

    ba_observations = observations
    ba_observations[,,,1] = ba_N
  } else {
    environment = X
    observations = Y
    ba_observations = observations
    timesteps = dim(X)[2]
    sites = dim(X)[1]
    species = dim(Y)[3]
  }
  out = list()

  out$get_model = function()  pkg.env$FINN$FINN(sp = dim(observations)[3],
                                                env = dim(environment)[3],
                                                device = device,
                                                parGlobal = init_global,
                                                parGrowth = init_growth,
                                                parMort = init_mort,
                                                parReg = init_reg,
                                                hidden_growth = hidden_growth,
                                                hidden_mort = hidden_mort,
                                                hidden_reg = hidden_reg,
                                                parGrowthEnv = init_growth_env,
                                                parMortEnv = init_mort_env,
                                                parRegEnv = init_reg_env
                                                )
  out$model = out$get_model()

  if(!simulate) {
    out$model$fit(X = environment,
                  Y = ba_observations,
                  response = "ba*T",
                  epochs = as.integer(epochs),
                  batch_size = as.integer(batch_size),
                  start_time = start_time,
                  patches = as.integer(patches),
                  learning_rate = learning_rate)
  }
  out$data = list(env = environment, ba_observations = ba_observations, observations = observations)
  class(out) = "FINN"
  return(out)
}



#' Predict from a fitted dnn model
#'
#' @param object a model created by \code{\link{FINN}}
#' @param newdata new data for predictions
#' @param device device on which network should be trained on.
#' @param ... additional arguments passed to model$predict function
#' @return prediction array
#'
#' @export
predict.FINN <- function(object, newdata = NULL,
                            type=c("link", "response", "class"),
                            device = c("cpu","cuda", "mps"),
                            reduce = c("mean", "median", "none"),...) {

  # object = checkModel(object)

  if(is.null(newdata)) newdata = object$data$env

  if(is.data.frame(newdata)) {
    data <- stats::model.matrix(object$call$formula, newdata)
    XX = as.data.frame(cbind(timestep = data$timestep, site = data$site, X))
    timesteps = length(unique(XX$timestep))
    sites = length(unique(XX$site))

    env_array =
      lapply(1:sites, function(site) {
        tmp_site = XX[XX$site == site, ]
        tmps =
          lapply(1:timesteps, function(time) {
            tmp =
              tmp_site[tmp_site$timestep == time, ] %>% select(-timestep, -site) %>% c()
            return(matrix(unlist(tmp), length(tmp), 1))
          })
        return(abind::abind(tmps, along = 0L))
      })
    environment = abind::abind(env_array, along = 0L)
    newdata = environment[,,,1]
  }

  predictions = reticulate::py_to_r(object$model$predict(env = newdata, ...))
  predictions = lapply(predictions, function(d) force_r(d$cpu()$data$numpy()))
  return(predictions)
}


#' Returns list of parameters
#'
#' @param object a model created by \code{\link{FINN}}
#' @param ... nothing implemented yet
#' @return list of weights
#'
#' @export
coef.FINN <- function(object,...){
  pars = list()
  pars$global = FINN:::force_r(m$model$parGlobal)
  pars$growth = FINN:::force_r(m$model$parGrowth)
  pars$mort = FINN:::force_r(m$model$parMort)
  pars$reg = FINN:::force_r(m$model$parReg)
  pars$ENVgrowth = FINN:::force_r(m$model$GrowthEnv)
  pars$ENVmort = FINN:::force_r(m$model$MortEnv)
  pars$ENVreg = FINN:::force_r(m$model$RegEnv)
  return(pars)
}

#' Make FINN arrays from observation data frame
#'
#' @param cohort_df a data frame with the columns Species, dbh, nTree. Optionally cohortID an be provided.
#' @param ... nothing implemented yet
#' @return list of cohort arrays dbh, nTree, Species from data frame with dendrometric data
#'
#' @export
cohort_df2arrays <- function(cohort_df, ...){
  Nspecies = length(unique(cohort_df$Species))
  Ncohorts = nrow(cohort_df)
  Species = array(data = cohort_df$Species, dim = c(1,1,Ncohorts))
  dbh = array(data = cohort_df$dbh,dim = c(1,1,Ncohorts))
  nTree = array(data = cohort_df$nTree,dim = c(1,1,Ncohorts))
  return(list(dbh = dbh, nTree = nTree, Species = Species))
}
