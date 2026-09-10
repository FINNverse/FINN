library(testthat)
library(FINN)
library(data.table)

# FINN::regeneration_adult: recruitment limited by the density of conspecific
# adults (stems with dbh above a trainable threshold x). The unit tests call the
# process function on hand-built stands through the model it is bound to; the
# last tests fit a tiny FIA subset. Torch-gated.

build_adult_model <- function(Nsp, logx = log(20), logA = log(20)) {
  FINN.seed(1)
  finn(
    N_species            = Nsp, recruits_dbh = 12.9, reg_floor = 1e-3,
    competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    growth_process       = createProcess(~ temp + prec, FINN::growth,    optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = createProcess(~ temp + prec, FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
    regeneration_process = createProcess(
      ~ temp + prec, FINN::regeneration_adult,
      custom_parameters = list(reg_logK       = rep(log(50), Nsp),
                               reg_adult_logx = rep(logx, Nsp),
                               reg_adult_logA = rep(logA, Nsp)),
      optimizeSpecies = TRUE, optimizeEnv = TRUE)
  )
}

# one site, 2 patches, 2 species; cohort c has species sp[c], dbh d[c], n[c] stems
adult_call <- function(m, d, n, sp, debug = FALSE, stand = TRUE) {
  m$patch_size_ha <- 0.1
  shp <- c(1L, 2L, length(d))
  dbh   <- torch_tensor(array(rep(d,  each = 2), dim = shp))
  trees <- torch_tensor(array(rep(n,  each = 2), dim = shp))
  species <- torch_tensor(array(rep(sp, each = 2), dim = shp), dtype = torch_int64())
  pred  <- torch_zeros(1L, 2L)                        # exp(0) = 1 environment
  light <- torch_ones(1L, 2L, 2L)                     # full light
  args <- list(species = species, parReg = torch_tensor(c(0.1, 0.1)), pred = pred, light = light, debug = debug)
  if (stand) args <- c(args, list(dbh = dbh, trees = trees))
  do.call(m$regeneration_func, args)
}

test_that("no adults gives the floor; recruitment rises with conspecific adults only", {
  skip_if_no_torch()
  m <- build_adult_model(2L)
  floor_capped <- 50 * 1e-3 / (50 + 1e-3)

  none <- as.numeric(adult_call(m, d = c(30, 30), n = c(0, 0), sp = c(1L, 2L)))
  expect_equal(none, rep(floor_capped, 4), tolerance = 1e-6)

  small <- as.numeric(adult_call(m, d = c(14, 14), n = c(5, 0), sp = c(1L, 2L), debug = TRUE)$adults)
  expect_equal(small[1], 5 * 2 / 0.2 * plogis(-3), tolerance = 1e-4)  # 14 cm, x = 20, tau = 2: 5 % adult

  few  <- adult_call(m, d = c(30, 30), n = c(1, 0), sp = c(1L, 2L), debug = TRUE)
  many <- adult_call(m, d = c(30, 30), n = c(10, 0), sp = c(1L, 2L), debug = TRUE)
  # adults per ha of the whole site: n stems in each of 2 patches of 0.1 ha
  expect_equal(as.numeric(many$adults)[1], 10 * 2 / 0.2 * plogis(5), tolerance = 1e-4)
  expect_gt(as.numeric(many$mean[1, 1, 1]), as.numeric(few$mean[1, 1, 1]))
  expect_equal(as.numeric(many$mean[1, 1, 2]), floor_capped, tolerance = 1e-6)  # species 2 has no adults
})

test_that("without dbh/trees the seed source counts as present (f = 1)", {
  skip_if_no_torch()
  m <- build_adult_model(2L)
  out <- adult_call(m, d = c(30, 30), n = c(0, 0), sp = c(1L, 2L), debug = TRUE, stand = FALSE)
  expect_equal(as.numeric(out$seed), rep(1, 4))
})

test_that("gradients reach the adult threshold and the half-saturation density", {
  skip_if_no_torch()
  m <- build_adult_model(2L)
  out <- adult_call(m, d = c(21, 25), n = c(3, 2), sp = c(1L, 2L))
  out$sum()$backward()
  expect_true(all(abs(as.numeric(m$reg_adult_logx$grad)) > 0))
  expect_true(all(abs(as.numeric(m$reg_adult_logA$grad)) > 0))
  expect_true(all(as.numeric(m$reg_adult_logx$grad) < 0))  # higher threshold -> fewer adults -> fewer recruits
  expect_true(all(as.numeric(m$reg_adult_logA$grad) < 0))  # higher half-saturation -> fewer recruits
})

test_that("regeneration_adult fits, and its parameters move", {
  skip_if_no_torch()
  ext <- function(f) system.file("extdata", f, package = "FINN")
  obs <- fread(ext("fia_obs_dt.csv")); env <- fread(ext("fia_env_dt.csv")); it <- fread(ext("fia_init_trees.csv"))
  sites <- head(sort(unique(obs$siteID)), 4L)
  Nsp <- max(obs$species)
  m <- build_adult_model(Nsp)
  before <- list(x = as.numeric(m$reg_adult_logx), A = as.numeric(m$reg_adult_logA))
  suppressMessages(
    fit(m, env = env[siteID %in% sites], data = obs[siteID %in% sites],
        init_cohort = makeInitCohorts(it[siteID %in% sites], Nspecies = Nsp), device = "cpu",
        epochs = 3L, patches = 2, patch_size = 0.06, lr = 0.02, plot_progress = FALSE)
  )
  after <- list(x = as.numeric(m$reg_adult_logx), A = as.numeric(m$reg_adult_logA))
  expect_true(all(is.finite(unlist(after))))
  expect_false(isTRUE(all.equal(before$x, after$x)))
  expect_false(isTRUE(all.equal(before$A, after$A)))
})

test_that("regeneration_adult without its parameters errors clearly", {
  skip_if_no_torch()
  m <- finn(
    N_species = 2L, recruits_dbh = 12.9,
    competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    regeneration_process = createProcess(~ temp + prec, FINN::regeneration_adult,
                                         custom_parameters = list(reg_logK = rep(log(50), 2)),
                                         optimizeSpecies = TRUE, optimizeEnv = TRUE))
  expect_error(adult_call(m, d = c(30, 30), n = c(1, 1), sp = c(1L, 2L)), "reg_adult_logx, reg_adult_logA")
})
