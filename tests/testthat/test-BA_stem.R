# tests/testthat/test-BA_P.R
library(testthat)
library(FINN)
library(torch)

test_that("BA_stem calculates basal area correctly", {
  dbh <- torch_tensor(50)
  expected_result <- pi * (50 / 100 / 2)^2
  expect_equal(as.numeric(BA_stem(dbh)), expected_result, tolerance = 0.00001)
})


