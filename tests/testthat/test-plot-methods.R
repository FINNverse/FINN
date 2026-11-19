# plot() of a fitted model and of a simulation: both run from a bare call, draw
# what they promise, and skip a network response with a warning.

test_that("plot() of a model draws the parameters inside their bounds and the linear niches", {
  skip_if_no_torch()
  toy <- fit_toy_finn(epochs = 2L)
  m <- toy$model

  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  # the bare call draws the process parameters
  p <- plot(m)
  expect_s3_class(p, "ggplot")
  expect_true(all(c("lower", "upper", "at_bound") %in% names(p$data)))

  par_df <- p$data
  # every species-specific process parameter, on the constrained scale, inside its band
  expect_setequal(unique(par_df$process), c("growth", "mortality", "regeneration", "competition"))
  expect_true(all(par_df$value >= par_df$lower - 1e-6 & par_df$value <= par_df$upper + 1e-6))
  expect_true("recruit dispersion theta" %in% par_df$parameter)
  expect_equal(sum(par_df$parameter == "light requirement" & par_df$process == "growth"), toy$Nsp)
  g <- as.matrix(m$par_growth_r)
  expect_equal(par_df$value[par_df$process == "growth" & par_df$parameter == "size decay"], g[, 2],
               tolerance = 1e-6)

  # pars = "env", the argument and value of the 0.1.0 method, switches to the niches
  pe <- plot(m, pars = "env")
  expect_s3_class(pe, "ggplot")
  env_df <- pe$data
  expect_setequal(unique(as.character(env_df$process)), c("growth", "mortality", "regeneration"))
  expect_setequal(unique(env_df$predictor), c("temp", "prec"))
  # the x axis spans the training range of each predictor on its original scale
  expect_equal(range(env_df$x[env_df$predictor == "temp"]), range(toy$env$temp))
  # a linear response: the niche of every species is a straight line in each predictor
  one <- env_df[env_df$process == "growth" & env_df$predictor == "temp" & env_df$species == 1, ]
  expect_lt(max(abs(diff(diff(one$value)))), 1e-4)

  # subsets
  p1 <- plot(m, pars = "process", species = 1:2)
  expect_setequal(as.integer(as.character(p1$data$species)), 1:2)
  expect_error(plot(m, species = 99), "species")
  expect_error(plot(m, pars = "all"), "should be one of")
})

test_that("plot() of a model skips a network response with a warning", {
  skip_if_no_torch()
  toy <- fit_toy_finn(epochs = 1L)
  FINN.seed(1)
  m <- finn(N_species = toy$Nsp, recruits_dbh = 12.9,
            growth_process       = suppressMessages(createProcess(~ temp + prec, growth, hidden = 4L)),
            mortality_process    = createProcess(~ temp + prec, mortality),
            regeneration_process = createProcess(~ temp + prec, regeneration_capped))
  suppressMessages(fit(m, env = toy$env, data = toy$obs, init_cohort = toy$init, epochs = 1L,
                       patches = 4, patch_size = 0.06, plot_progress = FALSE))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_warning(p <- plot(m, pars = "env"), "neural network.*ALE")
  expect_false("growth" %in% as.character(p$data$process))
  # the process parameters of the network process are still drawn
  expect_true("growth" %in% plot(m)$data$process)
  # an unfitted model has no range to span
  expect_warning(plot(finn(N_species = toy$Nsp, recruits_dbh = 12.9), pars = "env"), "fit the model first")
})

test_that("a simulation has class finn_simulation and plot() draws it", {
  skip_if_no_torch()
  toy <- fit_toy_finn(epochs = 1L)
  sim <- simulateForest(toy$model, env = toy$env, init_cohort = toy$init, patches = 4, patch_size = 0.06)
  expect_s3_class(sim, "finn_simulation")
  # the tables inside are untouched
  expect_identical(class(sim$long$site), c("data.table", "data.frame"))
  expect_output(print(sim), "siteID")

  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  p <- plot(sim)
  expect_s3_class(p, "ggplot")
  expect_setequal(as.character(unique(p$data$variable)), c("ba", "dbh", "trees"))
  # the plotted mean is the mean over sites of the simulated value
  site_mean <- sim$long$site[variable == "ba" & species == 1 & year == 1, mean(value)]
  expect_equal(p$data[variable == "ba" & species == 1 & year == 1]$mean, site_mean)

  p2 <- plot(sim, response = c("growth", "mort"), sites = sort(unique(toy$env$siteID))[1])
  expect_setequal(as.character(unique(p2$data$variable)), c("growth", "mort"))
  expect_equal(unique(p2$data$n_sites), 1L)
  expect_error(plot(sim, response = "height"), "Unknown `response`")
})
