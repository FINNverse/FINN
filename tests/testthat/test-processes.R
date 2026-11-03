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
