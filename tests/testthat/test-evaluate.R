# evaluate(), comparePlot(), summary(), bare-ground simulation, and the fixes
# of #77 (cohorts after readRDS), #64 (env without covariates) and #11 (data.frame
# starting values).

eval_toy <- function(epochs = 3L) {
  Nsp <- max(FINN::fia_obs_dt$species); tr <- 1:4; ho <- 5:6
  obs <- FINN::fia_obs_dt; env <- FINN::fia_env_dt; it <- FINN::fia_init_trees
  init_tr <- suppressWarnings(makeInitCohorts(it[siteID %in% tr], Nspecies = Nsp))
  init_ho <- suppressWarnings(makeInitCohorts(it[siteID %in% ho], Nspecies = Nsp))
  FINN.seed(1)
  m <- finn(N_species = Nsp, recruits_dbh = 12.9,
            growth_process       = createProcess(~ temp + prec, growth),
            mortality_process    = createProcess(~ temp + prec, mortality),
            regeneration_process = createProcess(~ temp + prec, regeneration_capped))
  holdout <- list(data = obs[siteID %in% ho], env = env[siteID %in% ho], init_cohort = init_ho)
  suppressMessages(fit(m, data = obs[siteID %in% tr], env = env[siteID %in% tr], init_cohort = init_tr,
                       epochs = epochs, patches = 4, patch_size = 0.06, plot_progress = FALSE,
                       holdout = holdout))
  list(model = m, obs = obs[siteID %in% tr], env = env[siteID %in% tr], init = init_tr, holdout = holdout, Nsp = Nsp)
}

test_that("fit(holdout = ) stores the evaluation, and fit() keeps what evaluate() needs", {
  skip_if_no_torch()
  toy <- eval_toy()
  m <- toy$model
  expect_equal(m$patches, 4L)
  expect_equal(length(m$param_recent), 3L)
  expect_equal(m$param_recent_epochs, 1:3)
  expect_named(m$evaluation, c("train", "holdout"))
  sk <- m$evaluation$holdout$skill
  expect_setequal(as.character(sk$variable), c("dbh", "ba", "trees", "growth", "mort", "reg"))
  expect_true(all(c("r2", "rmse", "r", "bias", "ratio", "auc") %in% names(sk)))
  expect_true(is.finite(sk[variable == "mort"]$auc))
  expect_true(all(is.na(sk[variable != "mort"]$auc)))
})

test_that("match_simulation() aggregates as fit() scores: states at the census, rates compounded, recruits summed", {
  skip_if_no_torch()
  toy <- eval_toy(1L)
  sim <- simulateForest(toy$model, env = toy$env, init_cohort = toy$init)
  op <- FINN:::match_simulation(sim, toy$obs, quiet = TRUE)
  s <- sim$long$site
  row <- op[variable == "growth"][1]
  yrs <- seq(row$year - row$period_length + 1, row$year)
  g <- s[siteID == row$siteID & species == row$species & variable == "growth" & year %in% yrs]$value
  expect_equal(row$pred, exp(mean(log1p(g))) - 1, tolerance = 1e-6)
  rr <- op[variable == "reg"][1]
  r <- s[siteID == rr$siteID & species == rr$species & variable == "r_mean_ha" & year %in% yrs]$value
  expect_equal(rr$pred, sum(r), tolerance = 1e-6)
  b <- op[variable == "ba"][1]
  expect_equal(b$pred, s[siteID == b$siteID & species == b$species & variable == "ba" & year == b$year]$value)
  # the weighted AUC: a perfect ranking is 1, a reversed one 0, ties one half
  expect_equal(FINN:::auc_counts(c(0.1, 0.9), pos = c(0, 3), neg = c(5, 0)), 1)
  expect_equal(FINN:::auc_counts(c(0.9, 0.1), pos = c(0, 3), neg = c(5, 0)), 0)
  expect_equal(FINN:::auc_counts(c(0.5, 0.5), pos = c(0, 3), neg = c(5, 0)), 0.5)
})

test_that("comparePlot() and plot(sim, obs): scatter and time series, pooled and by species, long input", {
  skip_if_no_torch()
  toy <- eval_toy(1L)
  sim <- simulateForest(toy$model, env = toy$env, init_cohort = toy$init)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  p <- comparePlot(sim, toy$obs)
  expect_s3_class(p, "finn_comparison")
  expect_named(unclass(p), c("scatter", "mort", "reg"))
  # mortality is a ROC curve from (0, 0) to (1, 1); recruitment a ratio per named species
  roc <- ggplot2::layer_data(p$mort, 2)
  expect_equal(range(roc$x), c(0, 1)); expect_equal(range(roc$y), c(0, 1))
  rg <- p$reg$data
  op <- attr(p, "matched")[variable == "reg"]
  sp1 <- rg$species[1]
  expect_equal(rg$ratio[1], sum(op[species == sp1]$pred) / sum(op[species == sp1]$obs))
  expect_true(all(as.character(rg$name) %in% unique(toy$obs$species_name)))
  p2 <- plot(sim, toy$obs, type = "timeseries", species = 1:2)
  expect_setequal(unique(attr(p2, "matched")$species), 1:2)
  long <- data.table::melt(toy$obs, id.vars = c("siteID", "year", "species", "species_name", "period_length",
                                                 "n_at_risk", "n_died", "growth_n"))
  p3 <- comparePlot(sim, long, response = c("ba", "mort"))
  expect_setequal(as.character(attr(p3, "skill")$variable), c("ba", "mort"))
  expect_equal(attr(p3, "skill")[variable == "ba"]$r2, attr(p, "skill")[variable == "ba"]$r2)
  expect_error(comparePlot(sim, toy$obs, species = 99), "no matching observation")
})

test_that("the largest dbh from tree_data, pooled as the site maximum; the heatmap per species", {
  skip_if_no_torch()
  toy <- eval_toy(1L)
  td <- FINN::fia_tree_dt[siteID %in% 1:4]
  sim <- simulateForest(toy$model, env = toy$env, init_cohort = toy$init, return_cohorts = c(10L, 20L))
  op <- FINN:::match_simulation(sim, toy$obs, quiet = TRUE, tree_data = td)
  dm <- op[variable == "dbh_max"]
  expect_gt(nrow(dm), 0)
  # observed: the surviving trees grown to the census
  r <- dm[1]
  t1 <- td[siteID == r$siteID & year == r$year & species == r$species & died == 0]
  expect_equal(r$obs, max(t1$dbh * (1 + t1$growth)^t1$period_length), tolerance = 1e-8)
  # simulated: the largest cohort holding trees
  co <- data.table::dcast(sim$long$cohort[siteID == r$siteID & year == r$year & species == r$species],
                          siteID + patchID + year + species + cohortID ~ variable, value.var = "value")
  expect_equal(r$pred, max(co[trees >= 0.5]$dbh))
  # pooled over species: the site maximum
  pooled <- FINN:::pool_dbh_max(op, by_species = FALSE)[variable == "dbh_max"]
  expect_equal(nrow(pooled), data.table::uniqueN(dm[, .(siteID, year)]))
  s1 <- pooled[siteID == r$siteID & year == r$year]
  expect_equal(s1$obs, max(dm[siteID == r$siteID & year == r$year]$obs))
  # without recorded cohorts the largest dbh is left out, with a message
  sim0 <- simulateForest(toy$model, env = toy$env, init_cohort = toy$init)
  expect_message(op0 <- FINN:::match_simulation(sim0, toy$obs, tree_data = td), "return_cohorts")
  expect_false(any(op0$variable == "dbh_max"))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  h <- comparePlot(sim, toy$obs, type = "heatmap", tree_data = td)
  expect_s3_class(h, "finn_heatmap")
  expect_true("dbh_max" %in% as.character(h$ratio$data$variable))
  expect_true(all(as.character(h$skill$data$species_lab) %in% unique(toy$obs$species_name)))
})

test_that("summary() prints skill, convergence and the environment; print() takes the window", {
  skip_if_no_torch()
  toy <- eval_toy(4L)
  s <- summary(toy$model)
  expect_s3_class(s, "summary.finn_class")
  expect_named(s$skill, c("train", "holdout"))
  out <- capture.output(print(s, epochs_back = 2))
  expect_true(any(grepl("Skill", out)) && any(grepl("Convergence", out)) && any(grepl("Holdout", out)))
  cv <- FINN:::convergence_table(s$convergence, 2)
  expect_equal(cv$k_loss, 2L); expect_equal(cv$k_params, 2L)
  P <- s$convergence$params
  # a network block: the relative L2 change of its weights
  nn <- cv$params[cv$params$kind == "network weights", ][1, ]
  a <- P[[4]][[nn$parameter]]; b <- P[[2]][[nn$parameter]]
  expect_equal(nn$change_pct, 100 * sqrt(sum((a - b)^2)) / sqrt(sum(b^2)))
  # a species parameter (trained only with optimizeSpecies = TRUE, as in the
  # helper's model): by name, on the scale it acts on
  tm <- fit_toy_finn(epochs = 4L)$model
  s2 <- suppressMessages(suppressWarnings(summary(tm)))
  cvt <- FINN:::convergence_table(s2$convergence, 2)
  row <- cvt$params[cvt$params$parameter == "mortality: growth", ]
  expect_equal(nrow(row), 1L)
  mk <- s2$convergence$mechanistic$par_mortality_unconstrained
  on <- function(z) mk$lower[3] + (mk$upper[3] - mk$lower[3]) * stats::plogis(z[, 3])
  Pt <- s2$convergence$params
  x1 <- on(Pt[[4]]$par_mortality_unconstrained); x0 <- on(Pt[[2]]$par_mortality_unconstrained)
  expect_equal(row$change_pct, 100 * sqrt(sum((x1 - x0)^2)) / sqrt(sum(x0^2)))
  out2 <- capture.output(print(s2, epochs_back = 2))
  expect_true(any(grepl("The species parameter that moved most", out2)))
  # the table: one row per process and part, the parameter that moved most
  g <- FINN:::parameter_groups(cvt$params)
  expect_true(all(c("growth", "mortality", "regeneration", "competition") %in% g$process))
  expect_equal(g$part[g$parameter == "mortality: growth"], "species parameters")
  expect_equal(unique(g$part[grepl("^nn_growth", g$parameter)]), "environmental response")
  expect_true(any(grepl("environmental response", out2)))
  # without the per-epoch record, the snapshots of the fit are used
  conv2 <- s$convergence; conv2$params <- list()
  conv2$snapshots <- P; conv2$snap_every <- 1L
  cv2 <- FINN:::convergence_table(conv2, 2)
  expect_equal(cv2$params_source, "snapshot")
  expect_equal(cv2$k_params, 2L)
})

test_that("simulateForest(bare_ground = TRUE) runs from an empty stand under the mean climate", {
  skip_if_no_torch()
  toy <- eval_toy(1L)
  e <- FINN:::bare_ground_env(toy$model, years = 5)
  expect_equal(nrow(e), 5L)
  expect_equal(e$temp[1], mean(toy$env$temp))
  e2 <- FINN:::bare_ground_env(toy$model, sites = c(2, 4), years = 3)
  expect_equal(unique(e2$siteID), c(2, 4))
  expect_equal(e2[siteID == 4]$prec[1], mean(toy$env[siteID == 4]$prec))
  bg <- simulateForest(toy$model, bare_ground = TRUE, years = 10)
  expect_s3_class(bg, "finn_simulation")
  expect_true(isTRUE(attr(bg, "bare_ground")))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  pb <- plot(bg)
  expect_true("site_lab" %in% names(pb$layers[[2]]$data))
  # the names fit() kept from the data label the species
  expect_true(all(stats::na.omit(unique(pb$layers[[2]]$data$colour_group)) %in% unique(toy$obs$species_name)))
  expect_equal(max(bg$long$site$year), 10L)
  expect_equal(sum(bg$long$site[year == 1 & variable == "trees"]$value > 0) >= 0, TRUE)
  expect_error(simulateForest(toy$model), "`env` is needed")
  expect_error(simulateForest(toy$model, bare_ground = TRUE, init_cohort = toy$init), "empty stand")
  evaluate(toy$model, bare_ground = TRUE, years = 10)
  expect_false(is.null(toy$model$evaluation$bare_ground))
})

test_that("initial cohorts survive saveRDS()/readRDS() (#77)", {
  skip_if_no_torch()
  toy <- eval_toy(1L)
  f <- tempfile(fileext = ".rds"); saveRDS(list(init = toy$init), f)
  back <- readRDS(f)$init
  expect_false(FINN:::cohort_tensors_valid(back))
  expect_identical(back$dbh_r, toy$init$dbh_r)          # rebuilt on access
  expect_true(FINN:::cohort_tensors_valid(back))
  back2 <- readRDS(f)$init
  sim <- simulateForest(toy$model, env = toy$env, init_cohort = back2)   # rebuilt at the entry
  expect_s3_class(sim, "finn_simulation")
})

test_that("an env without covariates builds an intercept (#64); data.frame starting values are matrices (#11)", {
  env <- data.table::data.table(siteID = rep(1:3, each = 3), year = rep(1:3, 3))
  x <- FINN:::extract_env(~ 1, env)
  expect_false(anyNA(x))
  expect_equal(dim(x), c(3L, 3L, 1L))
  p <- createProcess(~ temp, growth, initSpecies = as.data.frame(matrix(0.5, 3, 2)),
                     initEnv = as.data.frame(matrix(0, 3, 2)))
  expect_true(is.matrix(p$initSpecies))
  expect_true(is.matrix(p$initEnv[[1]]))
})

test_that("observations are placed by the position of their year in env, as fit() does (calendar years)", {
  skip_if_no_torch()
  toy <- eval_toy(1L)
  td <- FINN::fia_tree_dt[siteID %in% 1:4]
  shift <- function(x) data.table::copy(x)[, year := year + 2000L]
  FINN.seed(3); sim  <- simulateForest(toy$model, env = toy$env, init_cohort = toy$init, return_cohorts = c(10L, 20L))
  FINN.seed(3); simc <- simulateForest(toy$model, env = shift(toy$env), init_cohort = toy$init, return_cohorts = c(10L, 20L))
  expect_equal(attr(simc, "env_years"), sort(unique(toy$env$year)) + 2000L)
  a <- FINN:::match_simulation(sim, toy$obs, quiet = TRUE, tree_data = td)
  b <- FINN:::match_simulation(simc, shift(toy$obs), quiet = TRUE, tree_data = shift(td))
  expect_equal(b$year, a$year + 2000L)
  expect_equal(b$pred, a$pred); expect_equal(b$obs, a$obs)
  expect_true(any(b$variable == "dbh_max"))
  # evaluate() records the cohorts at the steps of the census years
  evaluate(toy$model, train = list(data = shift(toy$obs), env = shift(toy$env), init_cohort = toy$init,
                                   tree_data = shift(td)))
  expect_true("dbh_max" %in% as.character(toy$model$evaluation$train$skill$variable))
  # the time series is drawn on the years of the observations
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  p <- plot(simc, shift(toy$obs), type = "timeseries")
  expect_equal(range(ggplot2::layer_data(p, 1)$x), range(toy$env$year) + 2000L)
  # a year that is not in env is left out, with a message
  late <- shift(toy$obs)[year == 2020L][, year := 2050L]
  expect_message(o <- FINN:::match_simulation(simc, rbind(shift(toy$obs), late)), "not a year of the simulated")
  expect_equal(nrow(o), nrow(b[variable != "dbh_max"]))
  # plot(sim, obs, n_species = ) reaches comparePlot()
  p2 <- plot(sim, toy$obs, response = "reg", n_species = 2)
  expect_equal(nrow(p2$reg$data), 2L)
})

test_that("match_simulation(): a missing period_length is one step, the rate of a cut window keeps its length, growth_rate = 'period' compounds", {
  skip_if_no_torch()
  toy <- eval_toy(1L)
  sim <- simulateForest(toy$model, env = toy$env, init_cohort = toy$init)
  s <- sim$long$site
  o <- data.table::copy(toy$obs)[year == 10L, period_length := NA]
  op <- FINN:::match_simulation(sim, o, quiet = TRUE)
  r <- op[variable == "growth" & year == 10L][1]
  expect_equal(r$period_length, 1)
  expect_equal(r$pred, s[siteID == r$siteID & species == r$species & variable == "growth" & year == 10L]$value, tolerance = 1e-6)
  expect_true(all(is.na(o[year == 10L]$period_length)))          # the caller's table is untouched
  # a window that reaches back before the first step: spread over the full period_length
  o15 <- data.table::copy(toy$obs)[year == 10L][, period_length := 15L]
  r <- FINN:::match_simulation(sim, o15, quiet = TRUE)[variable == "mort"][1]
  mm <- s[siteID == r$siteID & species == r$species & variable == "mort" & year <= 10L]$value
  expect_equal(r$pred, 1 - exp(sum(log1p(-mm)) / 15), tolerance = 1e-6)
  # the increment over the interval for a model with growth_rate = "period"
  simp <- sim; attr(simp, "growth_rate") <- "period"
  r <- FINN:::match_simulation(simp, toy$obs, quiet = TRUE)[variable == "growth"][1]
  g <- s[siteID == r$siteID & species == r$species & variable == "growth" &
           year %in% seq(r$year - r$period_length + 1, r$year)]$value
  expect_equal(r$pred, prod(1 + g) - 1, tolerance = 1e-6)
  # species names given as a factor
  expect_equal(unname(FINN:::species_names_lookup(2, data.frame(species = 1:2, species_name = factor(c("a", "b"))))),
               c("a", "b"))
})

test_that("a failing holdout evaluation keeps the fit; a model read with torch_load() rebuilds its training cohorts", {
  skip_if_no_torch()
  Nsp <- max(FINN::fia_obs_dt$species); tr <- 1:3
  init <- suppressWarnings(makeInitCohorts(FINN::fia_init_trees[siteID %in% tr], Nspecies = Nsp))
  FINN.seed(1)
  m <- finn(N_species = Nsp, recruits_dbh = 12.9,
            growth_process       = createProcess(~ temp + prec, growth),
            mortality_process    = createProcess(~ temp + prec, mortality),
            regeneration_process = createProcess(~ temp + prec, regeneration_capped))
  # a holdout whose observations lie in no year of its env
  bad <- list(data = data.table::copy(FINN::fia_obs_dt[siteID %in% 5:6])[, year := year + 100L],
              env = FINN::fia_env_dt[siteID %in% 5:6],
              init_cohort = suppressWarnings(makeInitCohorts(FINN::fia_init_trees[siteID %in% 5:6], Nspecies = Nsp)))
  expect_warning(suppressMessages(
    fit(m, data = FINN::fia_obs_dt[siteID %in% tr], env = FINN::fia_env_dt[siteID %in% tr], init_cohort = init,
        epochs = 1L, patch_size = 0.06, plot_progress = FALSE, holdout = bad)),
    "the evaluation with `holdout` failed")
  expect_equal(length(m$history), 1L)
  f <- tempfile(fileext = ".pt"); torch::torch_save(m, f)
  m2 <- torch::torch_load(f)
  expect_false(FINN:::cohort_tensors_valid(m2$train_init_cohort$cohort))
  ic <- FINN:::.finn_train_ic(m2)
  expect_true(FINN:::cohort_tensors_valid(ic))
  expect_identical(ic$dbh_r, init$dbh_r)
})

test_that("ALE works with a process function of one's own; parameter names; 'no obs' tiles", {
  skip_if_no_torch()
  Nsp <- max(FINN::fia_obs_dt$species); tr <- 1:3
  init <- suppressWarnings(makeInitCohorts(FINN::fia_init_trees[siteID %in% tr], Nspecies = Nsp))
  # a mortality function of one's own that does not record its inputs
  my_mort <- function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL) {
    torch::torch_sigmoid(pred + light * parMort[, 1][species] + parMort[, 2][species] * dbh / 100 +
                         growth * 100 * parMort[, 3][species])
  }
  FINN.seed(1)
  m <- finn(N_species = Nsp, recruits_dbh = 12.9,
            growth_process       = createProcess(~ temp + prec, growth),
            mortality_process    = createProcess(~ temp + prec, my_mort),
            regeneration_process = createProcess(~ temp + prec, regeneration_capped))
  suppressMessages(fit(m, data = FINN::fia_obs_dt[siteID %in% tr], env = FINN::fia_env_dt[siteID %in% tr],
                       init_cohort = init, epochs = 1L, patch_size = 0.06, plot_progress = FALSE))
  a <- ALE(m, plot = FALSE)
  expect_true("mortality" %in% names(a))
  expect_gt(nrow(a$mortality), 0)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  # FINN does not know the parameters of my_mort: numbered, or named by the user
  expect_true(any(plot(m)$data$parameter == "parameter 1"))
  p <- plot(m, parameter_names = list(mortality = c("shade", "size", "growth")))
  expect_true(all(c("shade", "size", "growth") %in% p$data$parameter[p$data$process == "mortality"]))
  expect_error(plot(m, parameter_names = list(mortality = "shade")), "has 1 names")
  # a species without observed recruitment: its ratio tile says so
  sim <- simulateForest(m, env = FINN::fia_env_dt[siteID %in% tr], init_cohort = init)
  obs <- data.table::copy(FINN::fia_obs_dt[siteID %in% tr])
  sp0 <- obs[reg > 0, species][1]
  obs[species == sp0, reg := 0]
  h <- comparePlot(sim, obs, type = "heatmap")
  lab <- h$ratio$data[species == sp0 & variable == "reg"]$tile_label
  expect_equal(lab, "no obs")
})

test_that("conditional effects leave out a process without recorded inputs, growth included", {
  skip_if_no_torch()
  Nsp <- max(FINN::fia_obs_dt$species); tr <- 1:3
  init <- suppressWarnings(makeInitCohorts(FINN::fia_init_trees[siteID %in% tr], Nspecies = Nsp))
  # a model saved with an earlier FINN version keeps a forward() that does not
  # record the inputs of a process function of one's own; mimicked here by a
  # mortality function that discards the record of growth at every step
  wiping_mort <- function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL) {
    if (isTRUE(self$record_raws)) self$raw_g <- NULL
    torch::torch_sigmoid(pred + light * parMort[, 1][species] + parMort[, 2][species] * dbh / 100 +
                         growth * 100 * parMort[, 3][species])
  }
  FINN.seed(1)
  m <- finn(N_species = Nsp, recruits_dbh = 12.9,
            growth_process       = createProcess(~ temp + prec, growth),
            mortality_process    = createProcess(~ temp + prec, wiping_mort),
            regeneration_process = createProcess(~ temp + prec, regeneration_capped))
  suppressMessages(fit(m, data = FINN::fia_obs_dt[siteID %in% tr], env = FINN::fia_env_dt[siteID %in% tr],
                       init_cohort = init, epochs = 1L, patch_size = 0.06, plot_progress = FALSE))
  expect_message(a <- ALE(m, plot = FALSE), "growth process recorded no inputs")
  expect_false("growth" %in% names(a))
  expect_gt(nrow(a$mortality), 0)
  expect_gt(nrow(a$regeneration), 0)
})

test_that("a CohortMat of an earlier version (no restore()) is used when valid, and named when not", {
  skip_if_no_torch()
  old <- structure(list(dbh = torch::torch_zeros(1, 1, 1), trees = torch::torch_zeros(1, 1, 1),
                        species = torch::torch_ones(1, 1, 1)), class = "CohortMat")
  expect_identical(FINN:::restore_init_cohort(old), old)
  dead <- structure(list(dbh = NULL, trees = NULL, species = NULL), class = "CohortMat")
  expect_error(FINN:::restore_init_cohort(dead), "earlier FINN version")
})

test_that("a bare-ground run of sites with the same mean climate is drawn, the panels told apart by site", {
  skip_if_no_torch()
  toy <- eval_toy(1L)
  env <- data.table::copy(toy$env)
  for (cc in setdiff(names(env), c("siteID", "year")))
    data.table::set(env, i = which(env$siteID == 2), j = cc, value = env[siteID == 1][[cc]])
  bg <- simulateForest(toy$model, bare_ground = TRUE, years = 5, sites = 1:2, env = env)
  expect_equal(anyDuplicated(unname(attr(bg, "site_labels"))), 2L)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  pb <- plot(bg)
  expect_equal(nlevels(pb$layers[[2]]$data$site_lab), 2L)
  expect_true(all(grepl("^site [12]: ", levels(pb$layers[[2]]$data$site_lab))))
})

test_that("the convergence check from snapshots says that their spacing is an estimate", {
  skip_if_no_torch()
  toy <- eval_toy(4L)
  s <- suppressMessages(suppressWarnings(summary(toy$model)))
  s$convergence$snapshots <- s$convergence$params; s$convergence$snap_every <- 2L
  s$convergence$params <- list()
  cv <- FINN:::convergence_table(s$convergence, 4)
  expect_equal(cv$k_steps, 2L); expect_equal(cv$k_params, 4L)
  out <- capture.output(print(s, epochs_back = 4))
  expect_true(any(grepl("parameter snapshot 2 before the last to the last one \\(about 4 epochs", out)))
})
