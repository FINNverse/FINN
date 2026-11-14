## The 0.2.0 interface of finn() and fit(): explicit default processes, the
## renamed and moved arguments (loss_weights, loss_recruit, finnControl(),
## createProcess(reg_floor), growth_rate), each still accepted under its old
## name with a warning and the same result, and the input checks of fit().
library(data.table)

api_setup <- function(sites = 1:3) {
  Nsp <- max(fia_obs_dt$species)
  list(obs = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites], Nsp = Nsp,
       init = makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp))
}

api_model <- function(Nsp, ...) {
  FINN.seed(1)
  finn(N_species = Nsp, recruits_dbh = 12.9,
       growth_process = createProcess(~ temp + prec, FINN::growth, optimizeSpecies = TRUE),
       mortality_process = createProcess(~ temp + prec, FINN::mortality, optimizeSpecies = TRUE),
       regeneration_process = createProcess(~ temp + prec, FINN::regeneration_capped, optimizeSpecies = TRUE),
       ...)
}

api_fit <- function(m, d, epochs = 2L, patches = 4, patch_size = 0.06, plot_progress = FALSE, ...) {
  FINN.seed(2)
  suppressMessages(fit(m, data = d$obs, env = d$env, init_cohort = d$init, epochs = epochs,
                       patches = patches, patch_size = patch_size, plot_progress = plot_progress, ...))
  invisible(m)
}

# replace tables of a setup (modifyList() would recurse into the data.tables)
api_set <- function(d, ...) {
  v <- list(...)
  d[names(v)] <- v
  d
}

api_pars <- function(m) lapply(m$parameters, function(p) as.matrix(p$cpu()))

# ---- finn() ------------------------------------------------------------------

test_that("the explicit defaults in finn()'s signature are the processes NULL used to mean", {
  skip_if_no_torch()
  Nsp <- 4
  expect_identical(formals(finn)$mortality_process, quote(createProcess(func = mortality)))
  expect_identical(formals(finn)$growth_process, quote(createProcess(func = growth)))
  expect_identical(formals(finn)$regeneration_process, quote(createProcess(func = regeneration_capped)))
  expect_identical(formals(finn)$competition_process, quote(createProcess(func = competition)))
  FINN.seed(3); a <- finn(N_species = Nsp)
  FINN.seed(3); b <- finn(N_species = Nsp, mortality_process = NULL, growth_process = NULL,
                          regeneration_process = NULL, competition_process = NULL)
  expect_identical(names(a$parameters), names(b$parameters))
  expect_equal(api_pars(a), api_pars(b), tolerance = 0)
  expect_identical(a$reg_cap, b$reg_cap)
  expect_identical(a$reg_floor, 1e-3)
  expect_identical(b$reg_floor, 1e-3)
  expect_false(a$growth_period_scale)
  for (p in c("mortality", "growth", "regeneration", "competition")) {
    expect_identical(deparse(a[[paste0(p, "_formula")]]), deparse(b[[paste0(p, "_formula")]]))
    expect_identical(a[[paste0("process_", p)]]$func, b[[paste0("process_", p)]]$func)
  }
})

test_that("reg_floor lives on the regeneration process; finn(reg_floor = ) still works with a warning", {
  skip_if_no_torch()
  Nsp <- 3
  expect_identical(formals(createProcess)$reg_floor, 1e-3)
  expect_identical(createProcess(func = FINN::regeneration_capped, reg_floor = 0.05)$reg_floor, 0.05)
  expect_error(createProcess(func = FINN::regeneration, reg_floor = -1), "reg_floor")
  m <- finn(N_species = Nsp, regeneration_process = createProcess(func = FINN::regeneration, reg_floor = 0.05))
  expect_identical(m$reg_floor, 0.05)
  m2 <- finn(N_species = Nsp, regeneration_process = createProcess(func = FINN::regeneration, reg_floor = c(0.1, 0.2, 0.3)))
  expect_identical(m2$reg_floor, c(0.1, 0.2, 0.3))
  expect_error(finn(N_species = Nsp, regeneration_process = createProcess(func = FINN::regeneration, reg_floor = c(0.1, 0.2))),
               "N_species")
  expect_warning(m3 <- finn(N_species = Nsp, reg_floor = 0.05), "createProcess\\(reg_floor = \\)")
  expect_identical(m3$reg_floor, 0.05)
  # the deprecated argument overrides the process' value, as documented
  expect_warning(m4 <- finn(N_species = Nsp, reg_floor = 0.07,
                            regeneration_process = createProcess(func = FINN::regeneration, reg_floor = 0.05)))
  expect_identical(m4$reg_floor, 0.07)
  expect_error(finn(N_species = Nsp, not_an_argument = 1), "unused argument")
})

test_that("growth_rate replaces growth_period_scale", {
  skip_if_no_torch()
  Nsp <- 3
  expect_false(finn(N_species = Nsp)$growth_period_scale)
  expect_true(finn(N_species = Nsp, growth_rate = "period")$growth_period_scale)
  expect_warning(m <- finn(N_species = Nsp, growth_period_scale = TRUE), "growth_rate")
  expect_true(m$growth_period_scale)
  expect_warning(m <- finn(N_species = Nsp, growth_period_scale = FALSE), "growth_rate")
  expect_false(m$growth_period_scale)
  expect_error(suppressWarnings(finn(N_species = Nsp, growth_rate = "period", growth_period_scale = TRUE)), "not both")
})

test_that("recruit_obs_weight is the last named argument of finn()", {
  nm <- names(formals(finn))
  expect_identical(nm[length(nm)], "...")
  expect_identical(nm[length(nm) - 1L], "recruit_obs_weight")
})

# ---- fit(): renamed arguments give the same fit with a warning ---------------

test_that("loss_weights + weight_scale = FALSE is weights + weight_scale = \"raw\"", {
  skip_if_no_torch()
  d <- api_setup()
  w <- c(1, 1, 1, 2, 1, 1)
  m_new <- api_fit(api_model(d$Nsp), d, loss_weights = w, weight_scale = FALSE)
  expect_warning(m_old <- api_fit(api_model(d$Nsp), d, weights = w, weight_scale = "raw"), "loss_weights")
  expect_equal(api_pars(m_new), api_pars(m_old), tolerance = 0)
  expect_identical(m_new$loss_weights, w)
  expect_identical(m_old$loss_weights, w)
  # TRUE is "auto", the default
  m_auto <- api_fit(api_model(d$Nsp), d, loss_weights = w)
  m_true <- api_fit(api_model(d$Nsp), d, loss_weights = w, weight_scale = TRUE)
  m_chr  <- api_fit(api_model(d$Nsp), d, loss_weights = w, weight_scale = "auto")
  expect_equal(api_pars(m_auto), api_pars(m_true), tolerance = 0)
  expect_equal(api_pars(m_auto), api_pars(m_chr), tolerance = 0)
  expect_false(isTRUE(all.equal(api_pars(m_auto), api_pars(m_new))))
  expect_error(api_fit(api_model(d$Nsp), d, weight_scale = "neither"), "weight_scale")
  expect_error(suppressWarnings(api_fit(api_model(d$Nsp), d, loss_weights = w, weights = w)), "not both")
})

test_that("loss_recruit = \"observed\" is recruits = \"observed\"", {
  skip_if_no_torch()
  d <- api_setup()
  d$obs <- copy(d$obs)[, period_length := ifelse(year == 1, 2, 1)]
  m_new <- api_fit(api_model(d$Nsp), d, loss_recruit = "observed")
  expect_warning(m_old <- api_fit(api_model(d$Nsp), d, recruits = "observed"), "loss_recruit")
  expect_equal(api_pars(m_new), api_pars(m_old), tolerance = 0)
  expect_identical(unname(m_new$loss_weights[6]), 0)
})

test_that("loss is loss_family, with a warning", {
  skip_if_no_torch()
  d <- api_setup()
  fam <- c(dbh = "mse", ba = "mse", trees = "nbinom", growth = "mse", mortality = "binomial", regeneration = "nbinom")
  m_new <- api_fit(api_model(d$Nsp), d, loss_family = fam)
  expect_warning(m_old <- api_fit(api_model(d$Nsp), d, loss = fam), "loss_family")
  expect_equal(api_pars(m_new), api_pars(m_old), tolerance = 0)
  # a partial vector keeps the other defaults
  m_part <- api_fit(api_model(d$Nsp), d, loss_family = c(trees = "nbinom"))
  expect_equal(api_pars(m_new), api_pars(m_part), tolerance = 0)
  expect_error(api_fit(api_model(d$Nsp), d, loss_family = c(stems = "poisson")), "not a response")
  expect_error(suppressWarnings(api_fit(api_model(d$Nsp), d, loss_family = fam, loss = fam)), "not both")
})

test_that("finnControl() replaces the optimizer arguments of fit()", {
  skip_if_no_torch()
  d <- api_setup()
  ctrl <- finnControl(lr_scheduler = "step", lr_scheduler_params = list(step_size = 1, gamma = 0.5), clip_norm = 1)
  expect_s3_class(ctrl, "finnControl")
  expect_output(print(ctrl), "lr_scheduler")
  m_new <- api_fit(api_model(d$Nsp), d, control = ctrl)
  expect_warning(
    expect_warning(
      expect_warning(m_old <- api_fit(api_model(d$Nsp), d, lr_scheduler = "step",
                                      lr_scheduler_params = list(step_size = 1, gamma = 0.5), clip_norm = 1),
                     "lr_scheduler"),
      "lr_scheduler_params"),
    "clip_norm")
  expect_equal(api_pars(m_new), api_pars(m_old), tolerance = 0)
  m_def <- api_fit(api_model(d$Nsp), d)
  expect_false(isTRUE(all.equal(api_pars(m_new), api_pars(m_def))))
  # the optimizer and its extra arguments, given the new and the old way
  m_a <- api_fit(api_model(d$Nsp), d, control = finnControl(optimizer = torch::optim_adam, weight_decay = 1e-3))
  expect_warning(m_b <- api_fit(api_model(d$Nsp), d, optimizer = torch::optim_adam, weight_decay = 1e-3), "optimizer")
  expect_equal(api_pars(m_a), api_pars(m_b), tolerance = 0)
  expect_error(api_fit(api_model(d$Nsp), d, control = list(lr_scheduler = "step")), "finnControl")
  expect_error(suppressWarnings(api_fit(api_model(d$Nsp), d, control = ctrl, clip_norm = 1)), "not both")
  expect_error(finnControl(lr_scheduler = "warmup"), "lr_scheduler")
})

test_that("the default fit is unchanged by the new argument layout", {
  skip_if_no_torch()
  d <- api_setup()
  m1 <- api_fit(api_model(d$Nsp), d)
  m2 <- api_fit(api_model(d$Nsp), d, loss_weights = "auto", weight_scale = TRUE,
                loss_recruit = "simulated", control = finnControl())
  expect_equal(api_pars(m1), api_pars(m2), tolerance = 0)
})

# ---- fit(): input checks -----------------------------------------------------

test_that("fit() names a missing or wrong column instead of failing deep inside", {
  skip_if_no_torch()
  d <- api_setup()
  m <- api_model(d$Nsp)
  expect_error(api_fit(m, api_set(d, obs = d$obs[, -"dbh"]), epochs = 1L), "`dbh`")
  expect_error(api_fit(m, api_set(d, obs = d$obs[, -"species"]), epochs = 1L), "`species`")
  expect_error(api_fit(m, api_set(d, obs = as.list(d$obs)), epochs = 1L), "data.frame")
  bad <- copy(d$obs)[1, species := 99L]
  expect_error(api_fit(m, api_set(d, obs = bad), epochs = 1L), "1..11")
  dup <- rbind(d$obs, d$obs[1])
  expect_error(api_fit(m, api_set(d, obs = dup), epochs = 1L), "duplicated")
  # a response that is switched off may be absent
  expect_no_error(api_fit(m, api_set(d, obs = d$obs[, -"dbh"]), epochs = 1L, loss_family = c(dbh = "none")))
  # and plain data.frames are accepted
  expect_no_error(api_fit(m, api_set(d, obs = as.data.frame(d$obs), env = as.data.frame(d$env)), epochs = 1L))
})

test_that("fit() and predict() check env and disturbance", {
  skip_if_no_torch()
  d <- api_setup()
  m <- api_model(d$Nsp)
  expect_error(api_fit(m, api_set(d, env = d$env[, .(siteID, year)]), epochs = 1L), "no predictor column")   # FINNverse/FINN#64
  expect_error(api_fit(m, api_set(d, env = d$env[, -"year"]), epochs = 1L), "`year`")
  na_env <- copy(d$env)[1, temp := NA_real_]
  expect_error(api_fit(m, api_set(d, env = na_env), epochs = 1L), "missing values")
  short <- d$env[!(siteID == 2 & year == 2)]
  expect_error(api_fit(m, api_set(d, env = short), epochs = 1L), "same years")
  late <- copy(d$obs)[year == 2, year := 5L]
  expect_error(api_fit(m, api_set(d, obs = late), epochs = 1L), "not in `env`")
  dist <- d$env[, .(siteID, year, intensity = 0)]
  expect_no_error(api_fit(m, d, epochs = 1L, disturbance = dist))
  expect_error(api_fit(m, d, epochs = 1L, disturbance = dist[, .(siteID, year)]), "`intensity`")
  expect_error(api_fit(m, d, epochs = 1L, disturbance = copy(dist)[1, intensity := 2]), "between 0 and 1")
  expect_error(api_fit(m, d, epochs = 1L, disturbance = dist[year < 2]), "sites and years of `env`")
  api_fit(m, d, epochs = 1L)
  expect_error(predict(m, env = d$env[, .(siteID, year)], init_cohort = d$init), "no predictor column")
  expect_no_error(predict(m, env = as.data.frame(d$env), init_cohort = d$init))
})

test_that("patches and patch_size are checked against init_cohort and the model (#66)", {
  skip_if_no_torch()
  d <- api_setup()
  m <- api_model(d$Nsp)
  expect_warning(api_fit(m, d, epochs = 1L, patches = 2), "`patches = 2` is ignored")
  expect_no_warning(api_fit(m, d, epochs = 1L, patches = 4))
  expect_identical(m$patch_size_ha, 0.06)
  expect_warning(api_fit(m, d, epochs = 1L, patch_size = 0.1), "0.06 ha the model was last fitted")
  expect_warning(api_fit(m, d, epochs = 1L), "0.1 ha the model was last fitted")   # back to 0.06
  # predict() takes the fitted patch size by default and warns about another one
  expect_no_warning(s <- predict(m, env = d$env, init_cohort = d$init))
  expect_identical(m$patch_size_ha, 0.06)
  expect_warning(predict(m, env = d$env, init_cohort = d$init, patch_size = 0.1), "differs from the 0.06 ha")
  expect_warning(simulateForest(m, env = d$env, init_cohort = d$init, patches = 7), "`patches = 7` is ignored")
  expect_no_warning(simulateForest(m, env = d$env, init_cohort = d$init, patch_size = 0.1))
  # a model that was never fitted simulates at 0.1 ha
  m0 <- api_model(d$Nsp)
  expect_null(m0$patch_size_ha)
})

test_that("start_time warns that it has no effect", {
  skip_if_no_torch()
  d <- api_setup()
  expect_warning(api_fit(api_model(d$Nsp), d, epochs = 1L, start_time = 2L), "start_time")
})

test_that("the progress plot draws every loss term and survives a term that is switched off", {
  skip_if_no_torch()
  d <- api_setup()
  m <- api_model(d$Nsp)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(api_fit(m, d, epochs = 3L, plot_progress = TRUE, loss_family = c(ba = "none")))
  losses <- do.call(rbind, m$history)
  expect_identical(dim(losses), c(3L, 7L))
  expect_no_error(FINN:::visualize.training(losses[, 1:6] / rep(losses[1, 1:6], each = 3), 3, 3, "Loss"))
})
