# tests/testthat/test-height.R
library(testthat)
library(FINN)

test_that("height calculates heights within a specified range", {
  # Generate the test data using expand.grid and store as data.frame
  test_df <- as.data.frame(expand.grid(
    dbh = seq(0, 500, 1),
    parHeight = seq(0, 1, 0.1)
  ))

  # Calculate the height and add it to the data frame
  test_df$height <- with(test_df, height(dbh, parHeight))

  # Define the minimum and maximum height
  min_height <- 0   # Lower limit of height
  max_height <- 135 # Upper limit of height

  dbh_height_positive = !any(unlist(by(test_df$height,test_df$parHeight,function(x) diff(x)<0)))

  # Perform the tests
  expect_true(dbh_height_positive, info = "All heights increase with dbh")
  expect_true(min(test_df$height) >= min_height, info = "All heights should be greater than or equal to the minimum height")
  expect_true(max(test_df$height) <= max_height, info = "All heights should be less than or equal to the maximum height")
})
