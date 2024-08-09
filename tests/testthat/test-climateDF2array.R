library(testthat)
library(data.table)
library(FINN)

test_that("climateDF2array function works correctly", {
  # Define the number of days in each month
  days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

  # Define the years and site IDs
  years <- 1970:2019
  site_ids <- 1:4

  # Create an empty list to store data
  climate_data_list <- list()

  # Loop over each year, month, and site to generate data
  for (site_id in site_ids) {
    for (year in years) {
      for (month in 1:12) {
        days <- 1:days_in_month[month]

        # Handle leap years
        if (month == 2 && (year %% 4 == 0 && (year %% 100 != 0 || year %% 400 == 0))) {
          days <- 1:29
        }

        # Create a data frame for the current site, year, and month
        climate_data <- data.frame(
          siteID = site_id,
          year = year,
          month = month,
          day = days,
          tmp = runif(length(days), -10, 30),  # Random temperature values
          pre = runif(length(days), 0, 200)    # Random precipitation values
        )

        # Add to the list
        climate_data_list <- append(climate_data_list, list(climate_data))
      }
    }
  }

  # Combine the list into a single data frame
  climate_dt <- rbindlist(climate_data_list)

  # Use the function to create the arrays
  resultDay <- climateDF2array(climate_dt, env_vars = c("tmp", "pre"))
  climate_dt_month <- climate_dt[, .(tmp = mean(tmp), pre = sum(pre)), by = .(siteID, year, month)]
  resultMonth <- climateDF2array(climate_dt_month, env_vars = c("tmp", "pre"))
  climate_dt_year <- climate_dt[, .(tmp = mean(tmp), pre = sum(pre)), by = .(siteID, year)]
  resultYear <- climateDF2array(climate_dt_year, env_vars = c("tmp", "pre"))

  # Start checking the arrays
  min_year <- min(climate_dt$year)
  check_vec <- c()

  for (site_i in 1:4) {
    # Check Yearly Data
    array1 <- as.matrix(resultYear[site_i,,])
    array2 <- as.matrix(climate_dt_year[siteID == site_i, .(tmp, pre)])
    array1 <- array1[!is.na(array1[,1]),]
    array2 <- array2[!is.na(array2[,1]),]
    dimnames(array1) <- NULL
    dimnames(array2) <- NULL
    check <- identical(round(array1, 4), round(array2, 4))
    if (!check) stop("Error: arrays are not equal for site ", site_i)
    check_vec <- c(check_vec, check)

    for (year_i in 1:uniqueN(climate_dt$year)) {
      # Check Monthly Data
      array1 <- as.matrix(resultMonth[site_i, year_i,,])
      array2 <- as.matrix(climate_dt_month[siteID == site_i & year == (year_i + min_year - 1), .(tmp, pre)])
      array1 <- array1[!is.na(array1[,1]),]
      array2 <- array2[!is.na(array2[,1]),]
      dimnames(array1) <- NULL
      dimnames(array2) <- NULL
      check <- identical(round(array1, 4), round(array2, 4))
      if (!check) stop("Error: arrays are not equal for site ", site_i, ", year ", year_i)
      check_vec <- c(check_vec, check)

      for (month_i in 1:12) {
        # Check Daily Data
        array1 <- as.matrix(resultDay[site_i, year_i, month_i,,])
        array2 <- as.matrix(climate_dt[siteID == site_i & year == (year_i + min_year - 1) & month == month_i, .(tmp, pre)])
        array1 <- array1[!is.na(array1[,1]),]
        array2 <- array2[!is.na(array2[,1]),]
        dimnames(array1) <- NULL
        dimnames(array2) <- NULL
        check <- identical(round(array1, 4), round(array2, 4))
        if (!check) stop("Error: arrays are not equal for site ", site_i, ", year ", year_i, ", month ", month_i)
        check_vec <- c(check_vec, check)
      }
    }
  }

  expect_true(all(check_vec == TRUE))
})

