## createProcess(upper =, lower =) must reach the model: the bounds become the
## par_<process>_upper / _lower buffers that constrain the species parameters.
## (They used to be overwritten with NULL inside createProcess, so every model
## silently used the defaults, e.g. [-10, 10] for mortality.)
testthat::test_that("createProcess passes species-parameter bounds to the model", {
  p <- createProcess(~1, FINN::mortality, lower = c(-3, -3, -3), upper = c(3, 3, 3), optimizeSpecies = TRUE)
  testthat::expect_equal(p$upper, c(3, 3, 3))
  testthat::expect_equal(p$lower, c(-3, -3, -3))
  m <- finn(N_species = 3, mortality_process = p)
  testthat::expect_equal(as.numeric(m$par_mortality_upper), c(3, 3, 3))
  testthat::expect_equal(as.numeric(m$par_mortality_lower), c(-3, -3, -3))
  par <- as.matrix(m$par_mortality)
  testthat::expect_true(all(par > -3 & par < 3))
})

testthat::test_that("without bounds the process defaults still apply", {
  m <- finn(N_species = 3, mortality_process = createProcess(~1, FINN::mortality, optimizeSpecies = TRUE))
  testthat::expect_equal(as.numeric(m$par_mortality_upper), c(10, 10, 10))
  testthat::expect_equal(as.numeric(m$par_mortality_lower), c(-10, -10, -10))
})
