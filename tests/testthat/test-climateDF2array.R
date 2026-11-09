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

## The site/year lookup of climateDF2arrayCpp() used to be a linear scan per row
## (O(rows x sites), ~90 s per call at 32,686 sites). This is that algorithm,
## written out in R, as the reference the hash lookup must reproduce exactly.
climate_cpp_reference <- function(m, site_ids, year_ids, include_month, include_day,
                                  Nsites, Nyears, Nmonths, Ndays, Nenv) {
  n <- if (include_day) Nsites * Nyears * Nmonths * Ndays * Nenv
       else if (include_month) Nsites * Nyears * Nmonths * Nenv else Nsites * Nyears * Nenv
  out <- rep(NA_real_, n)
  for (i in seq_len(nrow(m))) {
    site <- which(site_ids == m[i, 1])[1] - 1L    # first match; NA never matches
    year <- which(year_ids == m[i, 2])[1] - 1L
    if (is.na(site) || is.na(year)) stop("Site or year index not found. Check your data.")
    for (e in seq_len(Nenv) - 1L) {
      if (include_day) {
        mo <- trunc(m[i, 3] - 1); d <- trunc(m[i, 4] - 1)
        out[1 + site + Nsites * (year + Nyears * (mo + Nmonths * (d + Ndays * e)))] <- m[i, 5 + e]
      } else if (include_month) {
        mo <- trunc(m[i, 3] - 1)
        out[1 + site + Nsites * (year + Nyears * (mo + Nmonths * e))] <- m[i, 4 + e]
      } else {
        out[1 + site + Nsites * (year + Nyears * e)] <- m[i, 3 + e]
      }
    }
  }
  out
}

test_that("climateDF2arrayCpp's lookup reproduces the linear scan exactly", {
  set.seed(42)
  for (res in c("year", "month", "day")) {
    sites <- sample(c(3L, 17L, 101L, 5000L, 9L, 42L))   # unordered, non-contiguous ids
    years <- c(1990L, 1995L, 2001L, 2002L)
    g <- expand.grid(siteID = sites, year = years, month = if (res != "year") c(1L, 6L, 12L) else NA,
                     day = if (res == "day") c(1L, 15L, 31L) else NA)
    g <- g[sample(nrow(g)), , drop = FALSE]
    g <- g[-sample(nrow(g), 10), , drop = FALSE]       # missing site-years -> NA cells
    g <- rbind(g, g[1:3, ])                            # duplicated rows: the last one wins
    g$a <- rnorm(nrow(g)); g$b <- runif(nrow(g))
    g$a[c(2, 7)] <- NA                                 # NA values are copied as NA
    keys <- c("siteID", "year", if (res != "year") "month", if (res == "day") "day")
    m <- as.matrix(g[, c(keys, "a", "b")])
    s_ids <- sort(unique(g$siteID)); y_ids <- sort(unique(g$year))
    Nm <- if (res == "year") 1L else 12L; Nd <- if (res == "day") 31L else 1L
    args <- list(m, s_ids, y_ids, res == "month", res == "day", length(s_ids), length(y_ids), Nm, Nd, 2L)
    got <- do.call(FINN:::climateDF2arrayCpp, c(args, list(c("a", "b"))))
    want <- do.call(climate_cpp_reference, args)
    expect_identical(got, want)
    # an id that is not in the lookup stops, as before
    bad <- m; bad[1, 1] <- 0.5
    expect_error(do.call(FINN:::climateDF2arrayCpp, c(list(bad), args[-1], list(c("a", "b")))),
                 "Site or year index not found")
  }
})

test_that("climateDF2array places unordered rows and leaves missing site-years NA", {
  env <- data.table(expand.grid(siteID = c(30L, 10L, 20L), year = c(2002L, 2001L, 2003L)))
  env[, x := siteID * 1000 + year]
  env <- env[sample(.N)][!(siteID == 20L & year == 2002L)]
  a <- climateDF2array(env, env_vars = "x")
  expect_identical(dim(a), c(3L, 3L, 1L))
  want <- outer(c(10, 20, 30) * 1000, 2001:2003, "+"); want[2, 2] <- NA
  expect_identical(a[, , 1], want)
})


test_that("climateDF2array() stops on a duplicated (siteID, year) and names it (#35)", {
  dt <- data.table(siteID = c(1L, 1L, 2L, 2L), year = c(1L, 2L, 1L, 2L), temp = c(10, 11, 20, 21))
  expect_identical(drop(climateDF2array(dt, "temp")), matrix(c(10, 20, 11, 21), 2))
  dup <- rbind(dt, data.table(siteID = 1L, year = 1L, temp = -999))
  expect_error(climateDF2array(dup, "temp"), regexp = "appear more than once, e.g. \\(siteID 1, year 1\\)")
  # an exact copy is a duplicate too
  expect_error(climateDF2array(rbind(dt, dt[1]), "temp"), regexp = "appear more than once, e.g. \\(siteID 1, year 1\\)")
  # the bundled climate tables hold one row per site x year
  for (e in list(FINN::fia_env_dt, FINN::fia_env_test))
    expect_equal(anyDuplicated(data.table::as.data.table(e), by = c("siteID", "year")), 0L)
  # a missing combination stays NA: the response tables of fit() rely on it
  a <- drop(climateDF2array(dt[!(siteID == 2L & year == 2L)], "temp"))
  expect_true(is.na(a[2, 2]))
})
