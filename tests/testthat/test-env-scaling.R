library(testthat)
library(FINN)
library(data.table)

# compute_env_scaling() learns z-standardization constants (center/scale) from
# an env table; apply_env_scaling() re-applies the STORED constants so that
# calibration and prediction share one transformation. Both are pure R and run
# on CRAN.

test_that("compute_env_scaling learns per-variable center and scale", {
  env <- data.table(
    siteID = 1:5, year = 0L,
    temp   = c(1, 2, 3, 4, 5),
    prec   = c(10, 20, 30, 40, 50)
  )
  s <- compute_env_scaling(env)

  expect_setequal(s$variable, c("temp", "prec"))          # siteID/year excluded
  expect_equal(s$center[s$variable == "temp"], mean(env$temp))
  expect_equal(s$scale[s$variable == "prec"],  sd(env$prec))
})

test_that("apply_env_scaling z-standardizes and inverts back to the raw values", {
  env <- data.table(
    siteID = 1:5, year = 0L,
    temp   = c(1, 2, 3, 4, 5),
    prec   = c(10, 20, 30, 40, 50)
  )
  s <- compute_env_scaling(env)
  z <- apply_env_scaling(env, s)

  # scaled predictors are mean 0, sd 1; keys are untouched
  expect_equal(mean(z$temp), 0, tolerance = 1e-12)
  expect_equal(sd(z$prec),   1, tolerance = 1e-12)
  expect_equal(z$siteID, env$siteID)

  # un-scaling with the stored constants recovers the raw column exactly
  recovered <- z$temp * s$scale[s$variable == "temp"] + s$center[s$variable == "temp"]
  expect_equal(recovered, env$temp, tolerance = 1e-12)
})

test_that("apply_env_scaling uses stored constants, not the new data's own", {
  train <- data.table(siteID = 1:3, year = 0L, temp = c(1, 2, 3))
  s <- compute_env_scaling(train)

  # a new site with a different mean must be scaled by the TRAINING center/sd,
  # so its mean is NOT forced to 0 (the classic leak of scale()-in-formula)
  newdat <- data.table(siteID = 4:6, year = 0L, temp = c(11, 12, 13))
  z <- apply_env_scaling(newdat, s)
  expect_equal(z$temp, (newdat$temp - s$center) / s$scale, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(mean(z$temp), 0)))
})

test_that("constant predictor is centred only (scale = 1, no divide-by-zero)", {
  env <- data.table(siteID = 1:4, year = 0L, flat = rep(7, 4))
  s <- compute_env_scaling(env)
  expect_equal(s$scale, 1)
  z <- apply_env_scaling(env, s)
  expect_true(all(is.finite(z$flat)))
  expect_true(all(z$flat == 0))          # (7 - 7) / 1
})

test_that("apply_env_scaling is a no-op for NULL and errors on a missing var", {
  env <- data.table(siteID = 1:3, year = 0L, temp = c(1, 2, 3))
  expect_identical(apply_env_scaling(env, NULL), env)

  s <- compute_env_scaling(env)
  expect_error(apply_env_scaling(data.table(siteID = 1:3, year = 0L), s),
               "missing environmental variable")
})

test_that("compute_env_scaling returns NULL when there are no numeric predictors", {
  expect_null(compute_env_scaling(data.table(siteID = 1:3, year = 0L)))
})
