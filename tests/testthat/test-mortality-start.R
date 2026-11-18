skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)

## Without createProcess(initEnv = ) the mortality environment response starts
## at about 2 % per step: intercept qlogis(0.02), every other weight 0. Until
## 0.2.0 it started at torch's random weights, sigmoid(0) = 0.5 per step.
Nsp <- 2
env <- data.frame(siteID = 1L, year = 1:3, temp = 0)
nn_pars <- function(m) {
  simulateForest(m, env = env, patches = 2, patch_size = 0.1, device = "cpu")   # builds the networks
  m$nn_mortality$parameters
}

test_that("the linear mortality response starts at qlogis(0.02) with zero slopes", {
  m <- finn(N_species = Nsp)                           # default processes, formula ~ .
  w <- as.matrix(nn_pars(m)[[1]])
  expect_equal(dim(w), c(Nsp, 2L))                     # intercept, temp
  expect_equal(w[, 1], rep(qlogis(0.02), Nsp), tolerance = 1e-6)
  expect_equal(w[, 2], rep(0, Nsp))
  # the other processes keep torch's start
  expect_false(all(as.matrix(m$nn_growth$parameters[[1]]) == 0))
})

test_that("with hidden layers the output bias carries the start", {
  m <- finn(N_species = Nsp, mortality_process = createProcess(~ 1 + temp, FINN::mortality, hidden = c(5L)))
  p <- nn_pars(m)
  expect_equal(as.numeric(p[[length(p)]]), rep(qlogis(0.02), Nsp), tolerance = 1e-6)
})

test_that("initEnv and a formula without intercept are left alone", {
  m <- finn(N_species = Nsp, mortality_process = createProcess(~ 1 + temp, FINN::mortality,
                                                                  initEnv = matrix(c(-1, 0.5), Nsp, 2, byrow = TRUE)))
  expect_equal(as.matrix(nn_pars(m)[[1]]), matrix(c(-1, 0.5), Nsp, 2, byrow = TRUE), tolerance = 1e-6)
  m0 <- finn(N_species = Nsp, mortality_process = createProcess(~ 0 + temp, FINN::mortality))
  expect_false(all(as.matrix(nn_pars(m0)[[1]]) == 0))
})

test_that("an unfitted stand keeps most of its trees over ten years", {
  m <- finn(N_species = Nsp, recruits_dbh = 5)
  init <- makeInitCohorts(data.table(siteID = 1L, patchID = 1:2, species = 1L, dbh = 30, trees = 50),
                          Nspecies = Nsp)
  FINN.seed(1)
  s <- simulateForest(m, env = data.frame(siteID = 1L, year = 1:10, temp = 0), init_cohort = init,
                      patch_size = 0.1, device = "cpu")$wide$site
  trees_sp1 <- s[year == 10 & species == 1, trees] * 2   # per-patch mean x 2 patches
  expect_gt(trees_sp1, 50)                             # of 100: at most about 5 % a year died
})
