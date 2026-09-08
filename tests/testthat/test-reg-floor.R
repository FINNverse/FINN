skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)

## a hostile regeneration environment (exp(-20) ~ 0): the only recruitment left is the floor
Nsp <- 2
mk <- function(reg_floor) finn(N_species = Nsp, recruits_dbh = 5, reg_floor = reg_floor,
  competition_process  = createProcess(~0, func = FINN::competition),
  growth_process       = createProcess(~1 + env1, initEnv = matrix(c(0.5, 0.5, 0.4, -0.4), Nsp, 2),
                                       initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp, 2), func = FINN::growth),
  mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(-2.5, -2.5, 0.3, -0.3), Nsp, 2),
                                       initSpecies = matrix(c(0.2, -0.2, 0.1, 0.1, 0, 0), Nsp, 3), func = FINN::mortality),
  regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(-20, -20, 0, 0), Nsp, 2),
                                       initSpecies = c(0.1, 0.1), func = FINN::regeneration, sample_regeneration = FALSE))
init_trees <- data.table(siteID = rep(1:3, each = 8), patchID = rep(rep(1:2, each = 4), 3),
                         species = rep(1:2, 12), dbh = rep(c(10, 20, 30, 40), 6) + rep(0:2, each = 8) * 5,
                         treeName = paste0("t", 1:24), trees = 3)
ic  <- makeInitCohorts(init_trees, Nspecies = Nsp)
env <- data.table(expand.grid(siteID = 1:3, year = 1:6)); env[, env1 := (siteID - 2) * 0.5]
sim_rec <- function(reg_floor) {
  FINN.seed(1); m <- mk(reg_floor)
  s <- simulateForest(m, env = env, init_cohort = ic, patches = 2L, patch_size = 0.1, device = "cpu")
  s$long$site[variable == "r_mean_ha", mean(value, na.rm = TRUE)]
}

test_that("finn() exposes reg_floor with a small default", {
  expect_equal(formals(finn)$reg_floor, 1e-3)
  expect_equal(mk(0.05)$reg_floor, 0.05)
})

test_that("the recruitment floor is what reg_floor says, not a hard-coded 0.2", {
  r_hi <- sim_rec(0.2); r_lo <- sim_rec(1e-3)
  expect_true(is.finite(r_hi) && is.finite(r_lo))
  expect_equal(r_hi, 0.2, tolerance = 0.05)     # regP <= 1, exp(-20) drive: mean ~ floor
  expect_lt(r_lo, 0.01)
  expect_gt(r_hi / r_lo, 50)
})
