library(testthat)
library(FINN)

# Allometry / process helpers. `dbh2ba` and `height` are pure R and run
# everywhere; `BA_stem` is the torch counterpart of `dbh2ba` and is gated.

test_that("dbh2ba returns stem basal area (m^2) from dbh (cm)", {
  dbh <- c(0, 10, 25, 50, 100)
  expect_equal(dbh2ba(dbh), pi * (dbh / 200)^2, tolerance = 1e-10)
  expect_equal(dbh2ba(0), 0)
  expect_equal(dbh2ba(50), pi * 0.25^2, tolerance = 1e-8)   # ~0.196 m^2
})

test_that("BA_stand sums stem basal area per hectare", {
  # one 50 cm stem: pi*(0.25)^2 m^2 spread over 0.1 ha
  expect_equal(BA_stand(50, 1, 0.1), pi * 0.25^2 / 0.1, tolerance = 1e-8)
  # linear in tree count, inverse in patch area
  expect_equal(BA_stand(30, 4, 0.1), 4 * BA_stand(30, 1, 0.1), tolerance = 1e-10)
  expect_equal(BA_stand(30, 1, 0.05), 2 * BA_stand(30, 1, 0.1), tolerance = 1e-10)
})

test_that("height is non-negative, increasing in dbh, and ~0 at dbh 0", {
  ph <- 0.5
  d  <- seq(0, 200, by = 1)
  h  <- height(d, ph)
  expect_true(all(h >= 0))
  expect_true(all(diff(h) > 0))            # strictly increasing in dbh
  expect_lt(height(0, ph), 0.01)           # a seedling has ~no height
  expect_gt(height(50, 0.8), height(50, 0.3))  # taller allometry -> taller tree
})

test_that("BA_stem (torch) matches the numeric dbh2ba", {
  skip_if_no_torch()
  dbh <- c(10, 25, 50)
  expect_equal(as.numeric(BA_stem(torch::torch_tensor(dbh))),
               dbh2ba(dbh), tolerance = 1e-6)
})

test_that("light_mask takes the bin edges over the live cohorts only", {
  skip_if_no_torch()
  # Padding slots have height 0.001 and occupy the lower deciles, so with enough of
  # them every live cohort lands in the top bin and shades nothing.
  mk <- function(n_pad) {
    K <- 3L + n_pad
    d <- array(c(rep(c(15, 30, 55), 1), rep(0, n_pad)), dim = c(1, 1, K))
    n <- array(c(rep(1, 3), rep(0, n_pad)), dim = c(1, 1, K))
    list(dbh = torch::torch_tensor(d, dtype = torch::torch_float32()),
         trees = torch::torch_tensor(n, dtype = torch::torch_float32()),
         species = torch::torch_tensor(array(1L, dim = c(1, 1, K)), dtype = torch::torch_long()))
  }
  pc <- torch::torch_tensor(cbind(0.5, 1), dtype = torch::torch_float32())
  lt <- function(x, ...) as.array(competition(dbh = x$dbh, species = x$species, trees = x$trees,
                                              parComp = pc, patch_size_ha = 0.1, ...))[1, 1, 1:3]
  dense <- mk(0L); sparse <- mk(60L)
  # unmasked: the padding pushes the live cohorts into the top bin
  expect_false(isTRUE(all.equal(lt(dense), lt(sparse), tolerance = 1e-5)))
  # masked: the padding is invisible, so a sparse patch gets the dense patch's light
  expect_equal(lt(dense, light_mask = TRUE), lt(sparse, light_mask = TRUE), tolerance = 1e-6)
})

test_that("mortality scales growth by growth_scale (default 0.01, i.e. %/yr)", {
  skip_if_no_torch()
  # FINN binds `self` when it installs the process; here a stub does it
  mort <- FINN::mortality
  environment(mort) <- list2env(list(self = list(record_raws = FALSE)), parent = environment(FINN::mortality))
  sp  <- torch::torch_tensor(matrix(c(1L, 2L), 1, 2), dtype = torch::torch_long())
  dbh <- torch::torch_tensor(matrix(c(20, 40), 1, 2))
  g   <- torch::torch_tensor(matrix(c(0.005, 0.02), 1, 2))
  li  <- torch::torch_tensor(matrix(c(0.3, 0.8), 1, 2))
  pr  <- torch::torch_tensor(matrix(c(-4, -3), 1, 2))
  par <- torch::torch_tensor(matrix(c(-0.5, 0.2, 0.3, -0.1, -2, 1), 2, 3))
  m_default <- as.numeric(mort(dbh, sp, NULL, par, pr, li, growth = g))
  # the same model on the raw rate needs the growth coefficient multiplied by 100
  par_unscaled <- as.matrix(par); par_unscaled[, 3] <- par_unscaled[, 3] / 0.01
  m_raw <- as.numeric(mort(dbh, sp, NULL, torch::torch_tensor(par_unscaled), pr, li, growth = g, growth_scale = 1))
  expect_equal(m_default, m_raw, tolerance = 1e-6)
  logit <- as.numeric(pr) + c(0.005, 0.02) / 0.01 * c(-2, 1) + c(0.3, 0.8) * c(-0.5, 0.2) + c(0.2, 0.4) * c(0.3, -0.1)
  expect_equal(m_default, plogis(logit), tolerance = 1e-6)
})


test_that("createProcess(growth_scale =) is a reachable opt-out of the 0.01 default", {
  skip_if_no_torch()
  # The documented wrapper is not: FINN::mortality reads self$record_raws from its own
  # namespace environment and set_environment() rebinds the outer function only. The
  # value therefore has to arrive as the installed kernel's formal.
  mk <- function(gs) finn(N_species = 2L,
    competition_process  = createProcess(~0, func = FINN::competition),
    growth_process       = createProcess(~1, func = FINN::growth),
    mortality_process    = createProcess(~1, func = FINN::mortality, growth_scale = gs),
    regeneration_process = createProcess(~1, func = FINN::regeneration))
  expect_equal(as.numeric(formals(mk(NULL)$mortality_func)$growth_scale), 0.01)
  expect_equal(as.numeric(formals(mk(1)$mortality_func)$growth_scale), 1)
  # a custom kernel without the argument is left alone
  own <- function(dbh, species, trees, parMort, pred, light, growth = NULL) torch::torch_sigmoid(pred)
  m <- finn(N_species = 2L,
    competition_process  = createProcess(~0, func = FINN::competition),
    growth_process       = createProcess(~1, func = FINN::growth),
    mortality_process    = createProcess(~1, func = own, growth_scale = 1),
    regeneration_process = createProcess(~1, func = FINN::regeneration))
  expect_false("growth_scale" %in% names(formals(m$mortality_func)))
})


test_that("createProcess rejects a growth_scale that is not a positive finite number", {
  # The value is used as a bare divisor, so every wrong value is silent rather than loud:
  # a bare NA is logical and coerces to TRUE, so `growth / NA` divides by 1 and restores the
  # pre-0.2.0 objective; NA_integer_ divides by -2^31; 0 saturates mortality to {0, 1}.
  cp <- function(gs) createProcess(~1, func = FINN::mortality, growth_scale = gs)
  expect_error(cp(NA), "positive finite")
  expect_error(cp(NA_real_), "positive finite")
  expect_error(cp(0), "positive finite")
  expect_error(cp(c(0.01, 1)), "positive finite")
  expect_error(cp("0.01"), "positive finite")
  expect_silent(cp(NULL))
  expect_silent(cp(1))
})
