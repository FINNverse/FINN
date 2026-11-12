## The recruit dispersion is trained, and it collapses: lowering it cheapens the
## regeneration likelihood without improving the predicted mean, while the same
## parameter sets the variance of the simulated draw (FINNverse/FINN#63). The
## bound sits in the `par_theta_recruits` getter, so the draw and the likelihood
## cannot read different values.
test_that("theta_floor bounds par_theta_recruits wherever it is read", {
  skip_if_no_torch()
  m <- finn(N_species = 3, regeneration_process =
              createProcess(~1, func = FINN::regeneration_capped, theta_floor = 2))
  m$par_theta_recruits <- 0.3
  expect_equal(as.numeric(m$par_theta_recruits), rep(2, 3), tolerance = 1e-5)
  m$par_theta_recruits <- 5
  expect_equal(as.numeric(m$par_theta_recruits), rep(5, 3), tolerance = 1e-4)

  mu <- finn(N_species = 3, regeneration_process =
               createProcess(~1, func = FINN::regeneration_capped, theta_floor = 0))
  mu$par_theta_recruits <- 0.3
  expect_equal(as.numeric(mu$par_theta_recruits), rep(0.3, 3), tolerance = 1e-5)
  expect_null(mu$theta_floor)
})

test_that("the default is 0.5 and a fit cannot push theta below it", {
  skip_if_no_torch()
  library(data.table)
  expect_equal(createProcess(func = FINN::regeneration_capped)$theta_floor, 0.5)
  expect_equal(createHybrid()$theta_floor, 0.5)

  sites <- 1:3
  Nsp   <- max(fia_obs_dt$species)
  init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
  # a floor above the start value of 1, so a short fit suffices to test the bound
  FINN.seed(1)
  m <- finn(N_species = Nsp, recruits_dbh = 12.9, regeneration_process =
              createProcess(~1, func = FINN::regeneration_capped, theta_floor = 1.5))
  FINN.seed(2)
  fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
      init_cohort = init, epochs = 5L, patches = 4, patch_size = 0.06,
      plot_progress = FALSE)
  expect_true(all(as.numeric(m$par_theta_recruits) >= 1.5 - 1e-6))
})

test_that("theta_floor is validated at construction", {
  skip_if_no_torch()
  expect_error(createProcess(func = FINN::regeneration_capped, theta_floor = -1), "theta_floor")
  expect_error(createProcess(func = FINN::regeneration_capped, theta_floor = c(1, 2)), "theta_floor")
  expect_null(createProcess(func = FINN::regeneration_capped, theta_floor = NULL)$theta_floor)
})
