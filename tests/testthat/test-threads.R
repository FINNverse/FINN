library(testthat)
library(FINN)

# fit(threads =): torch's intra-op thread count during a fit (FINN#48). torch's
# native backend (macOS) accepts one change per session, so the tests that change
# the count run only where OpenMP lets it change back.

test_that("threads is checked before anything is fitted", {
  skip_if_no_torch()
  expect_error(fit(model = NULL, env = NULL, threads = 0), "positive integer")
  expect_error(fit(model = NULL, env = NULL, threads = c(1, 2)), "positive integer")
})

test_that("threads = NULL and an unchanged count leave torch alone", {
  skip_if_no_torch()
  now = torch::torch_get_num_threads()
  expect_null(FINN:::set_torch_threads(NULL))
  expect_null(FINN:::set_torch_threads(now))
  expect_identical(torch::torch_get_num_threads(), now)
})

test_that("with OpenMP the fit's count is set and the previous one restored", {
  skip_if_no_torch()
  skip_if_not(torch::backends_openmp_is_available())
  before = torch::torch_get_num_threads()
  target = if (before == 1L) 2L else 1L
  d = list(obs = data.table::copy(FINN::fia_obs_dt)[siteID %in% 1:4],
           env = data.table::copy(FINN::fia_env_dt)[siteID %in% 1:4],
           it  = data.table::copy(FINN::fia_init_trees)[siteID %in% 1:4])
  Nsp = max(FINN::fia_obs_dt$species)
  FINN.seed(1)
  m = finn(N_species = Nsp, recruits_dbh = 12.9,
           competition_process  = createProcess(~0, FINN::competition),
           growth_process       = createProcess(~1, FINN::growth),
           regeneration_process = createProcess(~1, FINN::regeneration),
           mortality_process    = createProcess(~1, FINN::mortality))
  # the count the forward pass sees, read from inside the fit by a wrapped process
  seen = NULL
  f = m$growth_func
  m$growth_func = function(...) { seen <<- torch::torch_get_num_threads(); f(...) }
  suppressMessages(suppressWarnings(
    fit(m, env = d$env, data = d$obs, init_cohort = makeInitCohorts(d$it, Nspecies = Nsp),
        device = "cpu", epochs = 1L, patch_size = 0.06, lr = 0, plot_progress = FALSE,
        threads = target)))
  expect_identical(seen, target)
  expect_identical(torch::torch_get_num_threads(), before)
})
