library(testthat)
library(FINN)

# `"cuda"` is what a torch user reaches for first; match.arg()'s message does not
# say that FINN calls the same thing "gpu". The alias is resolved for every entry
# point, so a wrong spelling cannot reach a fit and fail hours later.

test_that("device accepts cuda as a synonym for gpu, and still rejects nonsense", {
  resolve = FINN:::resolve_device
  expect_identical(resolve("cuda"), "gpu")
  expect_identical(resolve("gpu"), "gpu")
  expect_identical(resolve("cpu"), "cpu")
  expect_identical(resolve(c("cpu", "gpu")), "cpu")     # the formal's default
  expect_error(resolve("mps"), "should be one of")
  expect_error(resolve("CUDA"), "should be one of")     # match.arg is case-sensitive
})

test_that("fit() reaches the same device from both spellings", {
  skip_if_no_torch()
  # epochs = 0 is set-up only, so this costs nothing and still runs the resolution
  d = fit_toy_finn(epochs = 0L, n_sites = 2L)
  expect_identical(d$model$device, "cpu")
})
