## createProcess(mortality_floor =): the lower bound on the annual probability with
## which the forward pass draws dead trees. The default 0.01 is the published
## behaviour; a lower floor lets a low-mortality stand keep its predicted rate.
skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(torch)
library(data.table)

test_that("the default floor is the published 0.01 and the argument reaches the model", {
  expect_equal(createProcess(func = FINN::mortality)$mortality_floor, 0.01)
  expect_equal(createHybrid()$mortality_floor, 0.01)
  m <- finn(N_species = 2,
            mortality_process = createProcess(~1, func = FINN::mortality, mortality_floor = 1e-6))
  expect_equal(m$mortality_floor, 1e-6)
  expect_equal(finn(N_species = 2)$mortality_floor, 0.01)
  expect_error(createProcess(func = FINN::mortality, mortality_floor = 0), "mortality_floor")
  expect_error(createProcess(func = FINN::mortality, mortality_floor = c(0.01, 0.02)), "mortality_floor")
  expect_error(createProcess(func = FINN::mortality, mortality_floor = NA_real_), "mortality_floor")
})

test_that("a stand that should not die loses 1 %/yr at the default floor and nothing below it", {
  sites <- 1:4
  env <- data.table::copy(FINN::fia_env_dt)[siteID %in% sites]
  env <- rbindlist(lapply(1:10, function(y) copy(env[year == 1L])[, year := y]))
  it  <- data.table::copy(FINN::fia_init_trees)[siteID %in% sites]
  Nsp <- max(FINN::fia_obs_dt$species)
  init <- makeInitCohorts(it, Nspecies = Nsp)
  stems_kept <- function(floor, draw) {
    FINN.seed(1)
    # intercepts of -20 on regeneration and mortality: no recruits, and a predicted
    # annual mortality of ~2e-9, so every stem lost is lost to the floor
    m <- finn(N_species = Nsp, recruits_dbh = 12.9,
              competition_process  = createProcess(~0, FINN::competition),
              growth_process       = createProcess(~1, FINN::growth),
              regeneration_process = createProcess(~1, FINN::regeneration, initEnv = list(matrix(-20, Nsp, 1))),
              mortality_process    = createProcess(~1, FINN::mortality, initEnv = list(matrix(-20, Nsp, 1)),
                                                   mortality_draw = draw, mortality_floor = floor))
    s <- simulateForest(m, env = env, init_cohort = init, patches = 1L, patch_size = 0.06)
    x <- s$wide$site
    x[year == max(year), sum(trees)] / x[year == min(year), sum(trees)]
  }
  expect_lt(stems_kept(0.01, "gamma"), 0.95)     # nine years at the published floor of 1 %/yr
  expect_equal(stems_kept(1e-6, "poisson"), 1)   # the predicted ~0 is kept
})
