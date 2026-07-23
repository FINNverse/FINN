library(testthat)
library(FINN)

# --- pure argument-resolution logic (no torch needed, runs on CRAN) ----------
test_that("resolve_record_years maps return_cohorts to timesteps", {
  f <- FINN:::resolve_record_years
  # FALSE / NULL -> nothing
  expect_length(f(FALSE, 10), 0)
  expect_length(f(NULL, 10), 0)
  # TRUE / "all" -> every timestep
  expect_identical(f(TRUE, 10), 1:10)
  expect_identical(f("all", 10), 1:10)
  # "last" -> final timestep
  expect_identical(f("last", 10), 10L)
  # explicit vector -> sorted unique subset
  expect_identical(f(c(3, 1, 3, 7), 10), c(1L, 3L, 7L))
  # debug = TRUE forces all timesteps regardless of return_cohorts
  expect_identical(f(FALSE, 10, debug = TRUE), 1:10)
  # invalid inputs
  expect_error(f(999, 10), "between 1 and 10")
  expect_error(f(0, 10), "between 1 and 10")
  expect_error(f("nope", 10), "must be FALSE")
})

# --- end-to-end integration (needs the torch backend) ------------------------
skip_if_not(torch::torch_is_installed())
library(torch)

make_model <- function(Nsp = 3) {
  finn(N_species = Nsp,
    competition_process  = createProcess(~0, func = FINN::competition),
    growth_process       = createProcess(~1+env1,
      initEnv = matrix(c(runif(Nsp,0.1,0.5), runif(Nsp,-2,-0.5)), Nsp, 2),
      initSpecies = matrix(c(runif(Nsp,0,1), runif(Nsp,0.03,0.1)), Nsp, 2),
      func = FINN::growth),
    regeneration_process = createProcess(~1+env1,
      initEnv = matrix(c(runif(Nsp,1,5), runif(Nsp,-2,2)), Nsp, 2),
      initSpecies = runif(Nsp,0,1), func = FINN::regeneration),
    mortality_process    = createProcess(~1+env1,
      initEnv = matrix(c(runif(Nsp,-3,-2), runif(Nsp,-3,-2)), Nsp, 2),
      initSpecies = matrix(c(scale(runif(Nsp)), scale(runif(Nsp)), rep(0,Nsp)), Nsp, 3),
      func = FINN::mortality)
  )
}

env_dt <- data.table::data.table(
  siteID = rep(1:2, each = 15),
  year   = rep(1:15, times = 2),
  env1   = rnorm(30)
)

test_that("return_cohorts selects the requested timesteps", {
  set.seed(1); torch_manual_seed(1)
  m <- make_model()
  sim <- function(rc) simulateForest(m, env = env_dt, patches = 3,
                                     device = "cpu", return_cohorts = rc)

  expect_null(sim(FALSE)$long$cohort)                          # default: none
  expect_identical(sort(unique(sim(TRUE)$long$cohort$year)), 1:15)
  expect_identical(sort(unique(sim("last")$long$cohort$year)), 15L)
  expect_identical(sort(unique(sim(c(4, 9))$long$cohort$year)), c(4L, 9L))
})

test_that("recording cohorts does not change the site output", {
  set.seed(2); torch_manual_seed(2)
  m <- make_model()
  # warm up once so the env networks are built (they are created lazily on the
  # first simulate() call and consume RNG); otherwise the two seeded runs below
  # would enter the sampler from different RNG states for reasons unrelated to
  # cohort recording.
  invisible(simulateForest(m, env = env_dt, patches = 3, device = "cpu"))
  set.seed(3); torch_manual_seed(3)
  a <- simulateForest(m, env = env_dt, patches = 3, device = "cpu", return_cohorts = FALSE)
  set.seed(3); torch_manual_seed(3)
  b <- simulateForest(m, env = env_dt, patches = 3, device = "cpu", return_cohorts = c(4, 9))
  ord <- function(dt) dt[order(siteID, year, species, variable)]
  expect_equal(ord(a$long$site), ord(b$long$site))
})

test_that("cohort table carries the expected per-cohort columns", {
  set.seed(4); torch_manual_seed(4)
  m <- make_model()
  co <- simulateForest(m, env = env_dt, patches = 3, device = "cpu",
                       return_cohorts = "last")$wide$cohort
  expect_true(all(c("siteID","patchID","species","cohortID","year",
                    "trees","dbh","g","m") %in% names(co)))
  expect_gt(nrow(co), 0)
})
