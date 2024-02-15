#' FINN model
#'
#' @param formula formula
#' @param data environmental data, must contain the variables timestep and site
#' @param observed species dbh and number of trees, must contain the variables timestep, site, nTree, and Species
#' @param patches number of patches
#' @param epochs number of epochs (optimization)
#' @param learning_rate learning rate for stochastic gradient descent
#' @param batch_size batch_size for stochastic gradient descent
#' @param device cpu or cuda, if available
#' @param start_time when to record, fractional, from 0 to 1
#'
#'
#' @export
#' @author Maximilian Pichler

FINN = function(formula,
                data = NULL,
                observed = NULL,
                patches = 10L,
                epochs = 1L,
                learning_rate = 1.0,
                batch_size = 70L,
                device = "cpu",
                start_time = 0.7) {
  require(dplyr)
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
      tmp_site = observed[observed$site == site, ]
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

  model = pkg.env$FINN$FINN(sp = dim(observations)[3], env = dim(environment)[3], device = device)
  model$fit(X = environment,
            Y = ba_observations,
            response = "ba*T",
            epochs = as.integer(epochs),
            batch_size = as.integer(batch_size),
            start_time = start_time,
            patches = as.integer(patches),
            learning_rate = learning_rate)
  out = list()
  out$model = model
  out$data = list(env = environment, ba_observations = ba_observations, observations = observations)
  class(out) = "FINN"
  return(out)
}
