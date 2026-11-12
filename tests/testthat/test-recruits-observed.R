## fit(recruits = "observed") drives the stand with the observed recruitment of each
## interval instead of the regeneration process, freezes that process for the call
## and weights its loss term 0; simulateForest(recruits = "observed") feeds the same
## rates from a table. The rate of an interval is reg / period_length per ha and
## year, in every simulated year whose recruitment fit() sums for that observation.

ro_setup <- function(sites = 1:4) {
  Nsp  <- max(fia_obs_dt$species)
  obs  <- data.table::copy(fia_obs_dt[siteID %in% sites])
  # env runs over the years 0, 1, 2 (timesteps 1-3); the census at year 1 closes a
  # two-year interval (timesteps 1-2), the one at year 2 a one-year interval
  obs[, period_length := ifelse(year == 1, 2, 1)]
  list(obs = obs, env = fia_env_dt[siteID %in% sites], Nsp = Nsp,
       init = makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp))
}

ro_model <- function(Nsp, sample_regeneration = TRUE) {
  FINN.seed(1)
  finn(N_species = Nsp, recruits_dbh = 12.9,
       growth_process = createProcess(~ temp + prec, FINN::growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
       mortality_process = createProcess(~ temp + prec, FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
       regeneration_process = createProcess(~ temp + prec, FINN::regeneration_capped,
                                            optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                            sample_regeneration = sample_regeneration))
}

ro_fit <- function(m, d, epochs = 2L, ...) {
  FINN.seed(2)
  suppressMessages(fit(m, data = d$obs, env = d$env, init_cohort = d$init, epochs = epochs,
                       patches = 4, patch_size = 0.06, plot_progress = FALSE, ...))
}

ro_reg_names <- function(m) {
  nm <- names(m$parameters)
  nm[grepl("^par_regeneration|^nn_regeneration[.]", nm) | nm %in% c("reg_logK", "par_theta_recruits_raw")]
}

test_that("the default \"simulated\" is the fit without the argument", {
  skip_if_no_torch()
  d <- ro_setup()
  m1 <- ro_model(d$Nsp); ro_fit(m1, d)
  m2 <- ro_model(d$Nsp); ro_fit(m2, d, recruits = "simulated")
  expect_identical(m1$history, m2$history)
  expect_identical(lapply(m1$parameters, as.array), lapply(m2$parameters, as.array))
})

test_that("observed recruits freeze regeneration, train the rest and restore requires_grad", {
  skip_if_no_torch()
  d <- ro_setup()
  m <- ro_model(d$Nsp)
  # an earlier simulated fit leaves gradients and optimizer state behind
  ro_fit(m, d, epochs = 1L)
  reg_nm <- ro_reg_names(m)
  expect_true(all(c("par_regeneration_unconstrained", "nn_regeneration.0.weight",
                    "reg_logK", "par_theta_recruits_raw") %in% reg_nm))
  before <- lapply(m$parameters, function(p) as.array(p$detach()$clone()))
  rg     <- vapply(m$parameters, function(p) p$requires_grad, logical(1))

  ro_fit(m, d, epochs = 3L, recruits = "observed")
  after <- lapply(m$parameters, function(p) as.array(p$detach()))

  for (nm in reg_nm) expect_identical(after[[nm]], before[[nm]], info = nm)
  for (nm in c("par_growth_unconstrained", "par_mortality_unconstrained", "nn_growth.0.weight"))
    expect_gt(max(abs(after[[nm]] - before[[nm]])), 0)
  expect_identical(vapply(m$parameters, function(p) p$requires_grad, logical(1)), rg)
  # the regeneration term was not scored
  expect_equal(unname(m$loss_weights[6]), 0)
  expect_true(all(vapply(m$history, function(h) h[6] == 0, logical(1))))

  # nothing on the model remembers the mode: the next default fit trains regeneration
  ro_fit(m, d, epochs = 2L)
  moved <- lapply(m$parameters[reg_nm], function(p) as.array(p$detach()))
  expect_gt(max(abs(moved$par_regeneration_unconstrained - after$par_regeneration_unconstrained)), 0)
})

test_that("the fed recruitment sums to the observed reg over each interval", {
  skip_if_no_torch()
  d <- ro_setup(sites = 1:2)
  m <- ro_model(d$Nsp, sample_regeneration = FALSE)
  rec <- d$obs[, .(siteID, year, species, period_length)]
  # rates of 50-200 per ha and year: per-patch expectations (x 0.06 ha) of 3-12, well
  # above the 0.5 below which the straight-through rounding delivers no recruit at all
  set.seed(3)
  rec[, reg := round(stats::runif(.N, 50, 200)) * period_length]
  rec[species == 2, reg := 0]
  sim <- simulateForest(m, env = d$env, init_cohort = d$init, patches = 4, patch_size = 0.06,
                        recruits = "observed", obs_recruits = rec)
  out <- data.table::as.data.table(sim$wide$site)
  # timestep k = k-th year of env (0, 1, 2); census year 1 covers timesteps 1-2, year 2 timestep 3
  out[, census := ifelse(year <= 2L, 1L, 2L)]
  per_int <- out[, .(r_mean_ha = sum(r_mean_ha), reg_patch = sum(reg)), by = .(siteID, species, census)]
  chk <- merge(per_int, rec, by.x = c("siteID", "species", "census"), by.y = c("siteID", "species", "year"))
  expect_equal(nrow(chk), 2L * d$Nsp * 2L)
  # expected recruits per ha summed over the interval = the observed reg ...
  expect_equal(chk$r_mean_ha, chk$reg, tolerance = 1e-5)
  # ... and with sample_regeneration = FALSE every patch recruits its expectation, rounded
  # (checked only where that expectation is well above 1; species 2 is fed 0)
  rate <- rec[, .(siteID, species, census = year, rate = reg / period_length, period_length)]
  chk2 <- merge(out, rate, by = c("siteID", "species", "census"))
  expect_equal(chk2$reg, round(chk2$rate * 0.06), tolerance = 1e-5)
  expect_true(any(chk2$reg > 0))
})

test_that("per-site windows are fed per site, and years after a site's last census get 0", {
  skip_if_no_torch()
  d <- ro_setup(sites = 1:3)
  # site 2 has no census at year 2: per-site path, and its timestep 3 is never scored
  d$obs <- d$obs[!(siteID == 2 & year == 2)]
  m <- ro_model(d$Nsp)
  ro_fit(m, d, recruits = "observed", shuffle = FALSE)
  expect_true(m$per_site_period)
  p <- merge(data.table::as.data.table(m$pred$wide$site)[, census := ifelse(year <= 2L, 1L, 2L)],
             d$obs[, .(siteID, species, census = year, reg_obs = reg, period_length)],
             by = c("siteID", "species", "census"), all.x = TRUE)
  expect_equal(p[!is.na(reg_obs), r_mean_ha], p[!is.na(reg_obs), reg_obs / period_length], tolerance = 1e-5)
  expect_true(all(p[siteID == 2 & year == 3L, r_mean_ha] == 0))
})

test_that("observed recruits refuse what they cannot feed", {
  skip_if_no_torch()
  d <- ro_setup(sites = 1:2)
  m <- ro_model(d$Nsp)

  bad <- data.table::copy(d$obs); bad[1, reg := NA]
  expect_error(ro_fit(m, list(obs = bad, env = d$env, init = d$init), recruits = "observed"),
               "1 site x species x interval cell")
  # without period_length the census at year 1 scores timestep 2 alone: timestep 1 has no recruitment
  nopl <- data.table::copy(d$obs)[, period_length := NULL]
  expect_error(ro_fit(m, list(obs = nopl, env = d$env, init = d$init), recruits = "observed"),
               "before or between")
  ovl <- data.table::copy(d$obs)[year == 2, period_length := 2]
  expect_error(ro_fit(m, list(obs = ovl, env = d$env, init = d$init), recruits = "observed"),
               "more than one observed interval")
  expect_error(ro_fit(m, d, recruits = "observed", reanchor = TRUE), "reanchor")
  expect_error(ro_fit(m, d, recruits = "observed", teacher_forcing = TRUE), "teacher_forcing")

  expect_error(simulateForest(m, env = d$env, init_cohort = d$init, patches = 4, patch_size = 0.06,
                              recruits = "observed"), "needs `obs_recruits`")
  expect_error(simulateForest(m, env = d$env, init_cohort = d$init, patches = 4, patch_size = 0.06,
                              obs_recruits = d$obs), "only used with")
  # a simulation beyond the last observed interval
  short <- d$obs[year == 1]
  expect_error(simulateForest(m, env = d$env, init_cohort = d$init, patches = 4, patch_size = 0.06,
                              recruits = "observed", obs_recruits = short), "not an autonomous simulator")
})
