## optimizeTheta = FALSE keeps the recruit dispersion at dispersion_parameter for
## the whole fit: a buffer, not a parameter, so neither the optimizer nor a
## requires_grad reset can move it. The default trains it, as before.
test_that("optimizeTheta = FALSE makes theta a buffer at dispersion_parameter", {
  skip_if_no_torch()
  m <- finn(N_species = 3, regeneration_process =
              createProcess(~1, func = FINN::regeneration_capped, dispersion_parameter = 2,
                            optimizeTheta = FALSE))
  expect_false("par_theta_recruits_raw" %in% names(m$parameters))
  expect_true("par_theta_recruits_raw" %in% names(m$buffers))
  expect_equal(as.numeric(m$par_theta_recruits), rep(2, 3), tolerance = 1e-4)

  md <- finn(N_species = 3, regeneration_process = createProcess(~1, func = FINN::regeneration_capped))
  expect_true("par_theta_recruits_raw" %in% names(md$parameters))
  expect_true(createProcess(func = FINN::regeneration_capped)$optimizeTheta)
  expect_true(createHybrid()$optimizeTheta)
})

test_that("a fit leaves a fixed theta untouched and trains a free one", {
  skip_if_no_torch()
  library(data.table)
  sites <- 1:3
  Nsp   <- max(fia_obs_dt$species)
  init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
  run <- function(optimizeTheta) {
    FINN.seed(1)
    m <- finn(N_species = Nsp, recruits_dbh = 12.9, regeneration_process =
                createProcess(~1, func = FINN::regeneration_capped, optimizeSpecies = TRUE,
                              theta_floor = 0, optimizeTheta = optimizeTheta))
    FINN.seed(2)
    fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
        init_cohort = init, epochs = 3L, patches = 4, patch_size = 0.06,
        plot_progress = FALSE)
    as.numeric(m$par_theta_recruits)
  }
  expect_equal(run(FALSE), rep(1, Nsp), tolerance = 1e-5)
  expect_gt(max(abs(run(TRUE) - 1)), 1e-4)
})

test_that("optimizeTheta is validated and survives torch_save / torch_load", {
  skip_if_no_torch()
  expect_error(createProcess(func = FINN::regeneration_capped, optimizeTheta = NA), "optimizeTheta")
  expect_error(createHybrid(optimizeTheta = "no"), "optimizeTheta")
  m <- finn(N_species = 2, regeneration_process =
              createProcess(~1, func = FINN::regeneration_capped, dispersion_parameter = 3,
                            optimizeTheta = FALSE))
  f <- tempfile(fileext = ".pt")
  torch::torch_save(m, f)
  m2 <- torch::torch_load(f)
  expect_false("par_theta_recruits_raw" %in% names(m2$parameters))
  expect_equal(as.numeric(m2$par_theta_recruits), rep(3, 2), tolerance = 1e-4)
})
