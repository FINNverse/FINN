skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)

## init_cohort = NULL is bare ground: no stems at the start, every tree comes from
## the regeneration process (until 0.2.0 it was one 1-cm stem of every species in
## every patch)
Nsp <- 2
mk <- function(reg_intercept) suppressWarnings(finn(N_species = Nsp, recruits_dbh = 5, reg_floor = 0,
  competition_process  = createProcess(~0, func = FINN::competition),
  growth_process       = createProcess(~1 + env1, initEnv = matrix(c(0.5, 0.5, 0.4, -0.4), Nsp, 2),
                                       initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp, 2), func = FINN::growth),
  mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(-2.5, -2.5, 0.3, -0.3), Nsp, 2),
                                       initSpecies = matrix(c(0.2, -0.2, 0.1, 0.1, 0, 0), Nsp, 3), func = FINN::mortality),
  regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(reg_intercept, 0, 0), Nsp, 2),
                                       initSpecies = c(0.1, 0.1), func = FINN::regeneration, sample_regeneration = FALSE)))
env <- data.table(expand.grid(siteID = 1:3, year = 1:8)); env[, env1 := (siteID - 2) * 0.5]
stand <- function(m) {
  FINN.seed(1)
  s <- simulateForest(m, env = env, patches = 2L, patch_size = 0.1, device = "cpu")$long$site
  s[variable %in% c("ba", "trees"), .(value = sum(value)), by = .(variable, year, species)]
}

test_that("the bare-ground start holds no stems", {
  ic <- FINN:::bare_ground_cohorts(3L, 2L, Nsp)
  expect_equal(ic$dims, c(3L, 2L, 1L))
  expect_equal(sum(torch::as_array(ic$trees)), 0)
  expect_equal(sum(torch::as_array(ic$dbh)), 0)
  expect_equal(ic$sp, Nsp)
})

test_that("bare ground stays bare when the regeneration process recruits nothing", {
  s <- stand(mk(reg_intercept = c(-30, -30)))   # exp(-30): no recruitment, no floor
  expect_lt(max(s$value), 1e-8)
})

test_that("on bare ground only the species the regeneration process recruits appear", {
  s <- stand(mk(reg_intercept = c(4, -30)))     # species 1 recruits (~55 per ha and year), species 2 does not
  expect_gt(s[variable == "trees" & species == 1 & year == max(year), value], 0)
  expect_lt(s[variable == "trees" & species == 2, max(value)], 1e-8)
})
