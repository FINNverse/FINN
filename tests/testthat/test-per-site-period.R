## Per-site aggregation windows: sites with different remeasurement intervals
## (period_length) in one fit. The observation tables are built from the
## model's own fit-mode predictions (lr = 0, fixed seed, no shuffling), so a
## correctly aggregated window gives an exactly zero MSE.
skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)

Nsp <- 2; Nsites <- 3; Tmax <- 10; SEED <- 11
mk_model <- function() finn(N_species = Nsp, recruits_dbh = 5,
  competition_process  = createProcess(~0, func = FINN::competition),
  growth_process       = createProcess(~1 + env1, initEnv = matrix(c(0.5, 0.5, 0.4, -0.4), Nsp, 2),
                                       initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp, 2), func = FINN::growth),
  mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(-2.5, -2.5, 0.3, -0.3), Nsp, 2),
                                       initSpecies = matrix(c(0.2, -0.2, 0.1, 0.1, 0, 0), Nsp, 3), func = FINN::mortality),
  regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(1.5, 1.0, 0.5, -0.5), Nsp, 2),
                                       initSpecies = c(0.1, 0.1), func = FINN::regeneration, sample_regeneration = FALSE))
FINN.seed(SEED); m <- mk_model()
init_trees <- data.table(siteID = rep(1:Nsites, each = 8), patchID = rep(rep(1:2, each = 4), Nsites),
                         species = rep(1:2, 12), dbh = rep(c(10, 20, 30, 40), 6) + rep(0:2, each = 8) * 5,
                         treeName = paste0("t", 1:24), trees = 3)
ic  <- makeInitCohorts(init_trees, Nspecies = Nsp)
env <- data.table(expand.grid(siteID = 1:Nsites, year = 1:Tmax)); env[, env1 := (siteID - 2) * 0.5]
LOSS <- c(dbh = "mse", ba = "mse", trees = "mse", growth = "mse", mortality = "mse", regeneration = "mse")
## fit()'s default optimizer needs torch >= 0.14 (optim_ignite_adam); fall back on older torch
OPT  <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
run_fit <- function(obs, mode = "auto") {
  FINN.seed(SEED); mm <- m$clone(deep = TRUE); mm$period_mode <- mode
  fit(mm, env = env, data = obs, init_cohort = ic, patches = 2L, patch_size = 0.1, env_autoscale = FALSE,
      loss = LOSS, epochs = 1L, lr = 0, weights = rep(1, 6), batchsize = Nsites, shuffle = FALSE,
      optimizer = OPT, plot_progress = FALSE, device = "cpu")
  mm
}
sched_rows <- function(sched) rbindlist(lapply(seq_along(sched), function(s) {
  yrs <- sched[[s]]; data.table(siteID = s, year = yrs, period_length = diff(c(0, yrs)))
}))
predictions <- function(sched) {
  d <- sched_rows(sched)[, .(species = 1:Nsp), by = .(siteID, year, period_length)]
  d[, `:=`(dbh = NA_real_, ba = NA_real_, trees = NA_real_, growth = NA_real_, mort = NA_real_, reg = NA_real_,
           species_name = paste0("sp", species))]
  p <- dcast(run_fit(d, "per_site")$pred$long$site, siteID + year + species ~ variable)
  for (cc in c("siteID", "year", "species")) p[, (cc) := as.integer(as.character(get(cc)))]
  p[]
}
mk_obs <- function(sched) {
  pred <- predictions(sched); rows <- sched_rows(sched)
  rbindlist(lapply(seq_len(nrow(rows)), function(k) {
    s <- rows$siteID[k]; y <- rows$year[k]; p <- rows$period_length[k]
    o <- pred[siteID == s & year == y, .(siteID, year, species, dbh, ba, trees)]
    w <- pred[siteID == s & year > y - p & year <= y,
              .(growth = mean(growth), mort = mean(mort), reg = sum(r_mean_ha)), by = species]
    o <- merge(o, w, by = "species"); o[, `:=`(period_length = p, species_name = paste0("sp", species))]
  }))[order(siteID, year, species)]
}

test_that("constant period: per-site path reproduces the constant path exactly", {
  obs5 <- mk_obs(list(c(5, 10), c(5, 10), c(5, 10)))
  a <- run_fit(obs5, "constant"); b <- run_fit(obs5, "per_site"); d <- run_fit(obs5, "auto")
  expect_false(a$per_site_period); expect_true(b$per_site_period); expect_false(d$per_site_period)
  expect_true(all(a$history[[1]][1:6] < 1e-8))
  expect_true(all(b$history[[1]][1:6] < 1e-8))
})

test_that("mixed periods: each site is aggregated over its own window", {
  obsm <- mk_obs(list(c(4), c(6), c(3, 7)))   # 4-yr | 6-yr | 3 then 4 yr
  r <- run_fit(obsm, "auto")
  expect_true(r$per_site_period)
  expect_true(all(r$history[[1]][1:6] < 1e-8))
  # negative control: a wrong window length on one site must score > 0
  bad <- copy(obsm); bad[siteID == 2, period_length := 5]
  expect_true(all(run_fit(bad, "auto")$history[[1]][4:6] > 1e-10))
  # the constant path cannot represent this schedule
  expect_error(run_fit(obsm, "constant"))
})

test_that("mixed periods: gradients flow through the deferred backward", {
  obsm <- mk_obs(list(c(4), c(6), c(3, 7)))
  FINN.seed(99)
  mm <- finn(N_species = Nsp, recruits_dbh = 5,
    competition_process  = createProcess(~0, func = FINN::competition),
    growth_process       = createProcess(~1 + env1, func = FINN::growth, optimizeSpecies = TRUE),
    mortality_process    = createProcess(~1 + env1, func = FINN::mortality, optimizeSpecies = TRUE),
    regeneration_process = createProcess(~1 + env1, func = FINN::regeneration, optimizeSpecies = TRUE))
  p0 <- as.numeric(torch::as_array(mm$par_growth))
  fit(mm, env = env, data = obsm, init_cohort = ic, patches = 2L, patch_size = 0.1, env_autoscale = FALSE,
      loss = LOSS, epochs = 2L, lr = 0.01, weights = rep(1, 6), batchsize = Nsites, optimizer = OPT,
      plot_progress = FALSE, device = "cpu")
  expect_true(mm$per_site_period)
  expect_gt(max(abs(as.numeric(torch::as_array(mm$par_growth)) - p0)), 1e-6)
})
