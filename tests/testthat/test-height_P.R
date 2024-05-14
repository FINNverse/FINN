# tests/testthat/test-height_P.R
library(testthat)
library(data.table)
library(ggplot2)
library(FINN)

test_that("height_P calculates heights within a specified range", {

  test_dt =
    data.table(
      expand.grid(
          list(
            dbh = seq(0, 500,1),
            parHeight = seq(0,1,0.1)
          )
        )
      )
  test_dt$rowID = 1:nrow(test_dt)

  min_height = 0   # lower limit of height
  max_height = 135 # higher limit of height

  test_dt[, height := height_P(dbh, parHeight), by = rowID]

  expect_true(min(test_dt$height) >= min_height, info = "All heights should be greater than or equal to the minimum height")
  expect_true(max(test_dt$height) <= max_height, info = "All heights should be less than or equal to the maximum height")
})
