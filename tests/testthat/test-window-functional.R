skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)
library(torch)

## fit() compounds the simulated yearly rates of an observation window into the
## per-year rate makeObsData(rates = "annual") reports: growth as the geometric
## mean of 1 + g minus one, mortality as one minus the geometric mean survival.
## Until 0.2.0 it compared the arithmetic mean of the yearly rates.
m_years <- c(0.5, rep(0.01, 9))                       # one disturbance year in ten
m_obs   <- 1 - (0.5 * 0.99^9)^(1 / 10)                # 0.0754, what makeObsData() derives
g_years <- c(0.5, rep(0.01, 9))                       # one year of fast growth in ten
g_obs   <- (1.5 * 1.01^9)^(1 / 10) - 1                # 0.0508; the mean would be 0.059

as_window <- function(x) torch_tensor(array(x, dim = c(1, length(x), 1)), dtype = torch_float32())

test_that("the window helpers give the compounded rate, not the mean", {
  m <- as.numeric(FINN:::window_mortality(as_window(m_years), NULL, 10))
  expect_equal(m, m_obs, tolerance = 1e-6)
  expect_gt(abs(m - mean(m_years)), 0.01)              # the mean would be 0.059
  g <- as.numeric(FINN:::window_growth(as_window(g_years), NULL, 10))
  expect_equal(g, g_obs, tolerance = 1e-6)
  expect_gt(abs(g - mean(g_years)), 1e-3)
  # with a [sites, years, 1] mask and a per-site denominator, as the per-site path calls them
  w <- as_window(rep(1, 10)); p <- torch_tensor(matrix(10, 1, 1), dtype = torch_float32())
  expect_equal(as.numeric(FINN:::window_mortality(as_window(m_years), w, p)), m_obs, tolerance = 1e-6)
  expect_equal(as.numeric(FINN:::window_growth(as_window(g_years), w, p)), g_obs, tolerance = 1e-6)
})

test_that("a masked year counts as survival 1 / growth 0 over the full period_length", {
  w <- as_window(c(0, rep(1, 9)))                     # year 1 outside the window
  m <- as.numeric(FINN:::window_mortality(as_window(rep(0.2, 10)), w, 10))
  expect_equal(m, 1 - (0.8^9)^(1 / 10), tolerance = 1e-6)
  g <- as.numeric(FINN:::window_growth(as_window(rep(0.03, 10)), w, 10))
  expect_equal(g, (1.03^9)^(1 / 10) - 1, tolerance = 1e-6)
})

test_that("the gradient stays finite where the simulated mortality is exactly 1", {
  z <- torch_tensor(array(c(30, rep(qlogis(0.01), 9)), dim = c(1, 10, 1)),
                    dtype = torch_float32(), requires_grad = TRUE)
  m <- torch_sigmoid(z)
  expect_equal(as.numeric(m[1, 1, 1]), 1)             # float32 saturates above a logit of 17.3
  out <- FINN:::window_mortality(m, NULL, 10)
  expect_true(is.finite(as.numeric(out)))
  out$backward()
  expect_true(all(is.finite(as.numeric(z$grad))))
})

## end to end: a stand whose mortality is exactly the yearly rate the environment
## dictates (no species terms), one disturbance year in a ten-year window
Nsp <- 2
m_model <- function(z_dist) finn(N_species = Nsp, recruits_dbh = 5,
  competition_process  = createProcess(~0, func = FINN::competition),
  growth_process       = createProcess(~1 + env1, initEnv = matrix(c(0.5, 0.5, 0, 0), Nsp, 2),
                                       initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp, 2), func = FINN::growth),
  # logit = qlogis(0.01) + env1 * (z_dist - qlogis(0.01)): 0.01 in a normal year, sigmoid(z_dist) in the disturbance year
  mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(qlogis(0.01), qlogis(0.01), z_dist - qlogis(0.01), z_dist - qlogis(0.01)), Nsp, 2),
                                       initSpecies = matrix(0, Nsp, 3), optimizeSpecies = FALSE, func = FINN::mortality),
  regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(-30, -30, 0, 0), Nsp, 2),
                                       initSpecies = c(0.1, 0.1), func = FINN::regeneration,
                                       sample_regeneration = FALSE, reg_floor = 0))
init <- makeInitCohorts(data.table(siteID = 1L, patchID = 1:2, species = 1L, dbh = 30, trees = 500),
                        Nspecies = Nsp)
env  <- data.table(siteID = 1L, year = 1:10, env1 = c(1, rep(0, 9)))   # the disturbance in year 1
obs  <- function(mort) data.table(siteID = 1L, year = 10L, species = 1:2, dbh = NA_real_, ba = NA_real_,
                                  trees = NA_real_, growth = NA_real_, mort = c(mort, NA), reg = NA_real_,
                                  period_length = 10)
# env_autoscale = FALSE: the 0/1 indicator must reach the logit as it is
fit_once <- function(m, d) suppressWarnings(fit(m, data = d, env = env, init_cohort = init, epochs = 1L,
  patch_size = 0.1, plot_progress = FALSE, device = "cpu", loss_family = c(mortality = "mse"),
  loss_weights = rep(1, 6), weight_scale = FALSE, env_autoscale = FALSE))

test_that("fit() scores the window mortality as the compounded rate", {
  m <- m_model(qlogis(0.5)); FINN.seed(1); fit_once(m, obs(m_obs))
  expect_lt(m$loss_raw[5], 1e-8)                     # the compounded target is matched
  m <- m_model(qlogis(0.5)); FINN.seed(1); fit_once(m, obs(mean(m_years)))
  expect_gt(m$loss_raw[5], 1e-5)                     # the arithmetic mean is not
})

test_that("a window with a year of mortality exactly 1 gives a finite loss and finite gradients", {
  m <- m_model(30); FINN.seed(1); fit_once(m, obs(0.3))
  expect_true(is.finite(m$loss_raw[5]))
  expect_true(all(vapply(m$parameters, function(p) all(is.finite(as.numeric(p$detach()))), logical(1))))
})
