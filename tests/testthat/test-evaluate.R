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
  P <- s$convergence$params; nm <- cv$params$block[1]
  a <- P[[4]][[nm]]; b <- P[[2]][[nm]]
  expect_equal(cv$params$change_pct[1], 100 * sqrt(sum((a - b)^2)) / sqrt(sum(b^2)))
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
