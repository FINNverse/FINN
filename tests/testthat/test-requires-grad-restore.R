## A simulation turns every parameter's requires_grad off and must put back the
## state it found, not TRUE. Parameters are frozen in two places -- `create_nn()`
## for a hybrid process built with `optimize = FALSE`, and by hand wherever a
## calibration freezes one stage before fitting the next -- and a single
## simulation between two fit() calls used to make all of them trainable again.
test_that("simulating restores a frozen parameter instead of enabling everything", {
  skip_if_no_torch()
  library(data.table)
  Nsp <- 2
  m <- finn(N_species = Nsp, recruits_dbh = 5,
            growth_process       = createProcess(~1, func = FINN::growth),
            mortality_process    = createProcess(~1, func = FINN::mortality),
            regeneration_process = createProcess(~1, func = FINN::regeneration))
  m$par_theta_recruits_raw$requires_grad_(FALSE)
  frozen <- names(Filter(function(p) !p$requires_grad, m$parameters))
  expect_true("par_theta_recruits_raw" %in% frozen)

  it <- data.table(siteID = 1L, patchID = 1L, species = 1:2, dbh = c(10, 25),
                   treeName = c("t1", "t2"), trees = 4)
  env <- data.table(expand.grid(siteID = 1L, year = 1:5))[, env1 := 0]
  simulateForest(m, env = env, init_cohort = makeInitCohorts(it, Nspecies = Nsp),
                 patches = 2L, patch_size = 0.1, device = "cpu")

  expect_setequal(names(Filter(function(p) !p$requires_grad, m$parameters)), frozen)
})
