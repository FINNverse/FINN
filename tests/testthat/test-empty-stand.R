skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)

## A stand that dies out before its census is still a prediction: zero basal
## area, zero trees and the window's recruits are scored against what was
## observed. Until 0.2.0 a step at which the batch had no cohort left was not
## scored at all, so a die-off cost nothing.
Nsp <- 2
# mortality near 1 per step (the ceiling 0.999 applies), no regeneration at all
m_dieoff <- function() finn(N_species = Nsp, recruits_dbh = 5,
  competition_process  = createProcess(~0, func = FINN::competition),
  growth_process       = createProcess(~1 + env1, initEnv = matrix(c(0.5, 0.5, 0, 0), Nsp, 2),
                                       initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp, 2), func = FINN::growth),
  mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(12, 12, 0, 0), Nsp, 2),
                                       initSpecies = matrix(c(0.2, -0.2, 0.1, 0.1, 0, 0), Nsp, 3), func = FINN::mortality),
  regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(-30, -30, 0, 0), Nsp, 2),
                                       initSpecies = c(0.1, 0.1), func = FINN::regeneration,
                                       sample_regeneration = FALSE, reg_floor = 0))
init <- makeInitCohorts(data.table(siteID = 1L, patchID = 1:2, species = 1L, dbh = 30, trees = 50),
                        Nspecies = Nsp)
env  <- data.table(siteID = 1L, year = 1:5, env1 = 0)
# the census at year 5 found a living stand of species 1
obs  <- data.table(siteID = 1L, year = 5L, species = 1:2, dbh = c(30, NA), ba = c(5, 0),
                   trees = c(50, 0), growth = NA_real_, mort = NA_real_, reg = 0)

test_that("the stand is empty at the census", {
  m <- m_dieoff()
  FINN.seed(1)
  s <- simulateForest(m, env = env, init_cohort = init, patch_size = 0.1, device = "cpu")$wide$site
  expect_equal(sum(s[year == 5, trees]), 0)
  expect_equal(sum(s[year == 5, ba]), 0)
})

test_that("an empty stand at a census is scored: ba and trees loss positive and finite", {
  m <- m_dieoff()
  FINN.seed(1)
  suppressWarnings(fit(m, data = obs, env = env, init_cohort = init, epochs = 1L,
                       patch_size = 0.1, plot_progress = FALSE, device = "cpu"))
  # the raw loss vector of the scored step: dbh, ba, trees, growth, mort, reg, (unused)
  expect_length(m$loss_raw, 7L)
  expect_true(all(is.finite(m$loss_raw[2:3])))
  expect_gt(m$loss_raw[2], 0)   # ba: 0 predicted against 5 observed
  expect_gt(m$loss_raw[3], 0)   # trees: 0 predicted against 50 observed
  # dbh, growth and mortality have no prediction without a tree and stay unscored
  expect_equal(m$loss_raw[c(1, 4, 5)], c(0, 0, 0))
  # and the step reaches the epoch's history
  h <- unlist(m$history[[1]])
  expect_true(all(is.finite(h)))
  expect_gt(sum(h[2:3]), 0)
})
