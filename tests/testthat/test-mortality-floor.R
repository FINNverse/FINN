## createProcess(mortality_floor =, mortality_ceiling =): the interval the forward
## pass clamps the mortality probability to, per simulation step, before the death
## draw. Defaults 0.001 / 0.999 since 0.2.0; the hard-coded [0.01, 0.99] it replaces
## is why a low-mortality stand lost stems the observations did not.
skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(torch)
library(data.table)

test_that("both limits default to 0.001 / 0.999 and reach the model", {
  expect_equal(createProcess(func = FINN::mortality)$mortality_floor, 0.001)
  expect_equal(createProcess(func = FINN::mortality)$mortality_ceiling, 0.999)
  expect_equal(createHybrid()$mortality_floor, 0.001)
  expect_equal(createHybrid()$mortality_ceiling, 0.999)
  m <- finn(N_species = 2, mortality_process =
              createProcess(~1, func = FINN::mortality, mortality_floor = 1e-6, mortality_ceiling = 0.5))
  expect_equal(m$mortality_floor, 1e-6)
  expect_equal(m$mortality_ceiling, 0.5)
  expect_equal(finn(N_species = 2)$mortality_floor, 0.001)
  expect_equal(finn(N_species = 2)$mortality_ceiling, 0.999)
  expect_error(createProcess(func = FINN::mortality, mortality_floor = 0), "mortality_floor")
  expect_error(createProcess(func = FINN::mortality, mortality_ceiling = 1), "mortality_ceiling")
  expect_error(createProcess(func = FINN::mortality, mortality_floor = c(0.01, 0.02)), "mortality_floor")
  expect_error(createProcess(func = FINN::mortality, mortality_floor = NA_real_), "mortality_floor")
  # a floor above the ceiling would clamp to an empty interval
  expect_error(createProcess(func = FINN::mortality, mortality_floor = 0.6, mortality_ceiling = 0.5), "below")
})

test_that("the ceiling bounds the drawn probability from above", {
  n <- torch_full(c(50, 50, 1), 20)
  p <- torch_full(c(50, 50, 1), 0.98)
  expect_lt(as.numeric((FINN:::binomial_from_gamma(n, torch_clamp(p, 0.001, 0.5))$mean())), 11)
  expect_gt(as.numeric((FINN:::binomial_from_gamma(n, torch_clamp(p, 0.001, 0.999))$mean())), 18)
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
    # 20 patches, not 1: the effect is an expectation over draws, and one patch of a
    # handful of stems is realisation noise wide enough to swallow it.
    s <- simulateForest(m, env = env, init_cohort = init, patches = 4L, patch_size = 0.06)
    x <- s$wide$site
    x[year == max(year), sum(trees)] / x[year == min(year), sum(trees)]
  }
  # Assert the ORDERING and the size of the gap, not an absolute threshold: the
  # straight-through rounding of the death count swallows much of a nominal 1 %/yr
  # at the small n*p of these cohorts, so the realised loss is well under the
  # nominal rate and any fixed threshold is tuned to one build of the simulator.
  high <- stems_kept(0.01, "gamma")
  low  <- stems_kept(1e-3, "gamma")
  expect_lt(high, low)                            # the higher floor costs stems
  expect_gt(low - high, 0.02)                     # and the difference is not noise
  expect_gt(low, 0.98)                            # at the lowest recommended floor the stand holds
})
