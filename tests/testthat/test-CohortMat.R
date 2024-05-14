# tests/testthat/test-CohortMat.R
library(testthat)
library(data.table)
library(R6)
library(torch)

# Define the CohortMat class as shown above if not sourced
# source("path/to/CohortMat.R")  # Make sure to source the class if defined elsewhere

test_that("CohortMat class correctly transforms data table to arrays and back", {
  # Initialize parameters
  Nsites <- 3
  Npatches <- 5
  maxNcohorts <- 20
  Nspecies <- 5

  # Create the data table
  obs_dt <- data.table()
  for (i_site in 1:Nsites) {
    for (i_patch in 1:Npatches) {
      Ncohorts <- sample(1:maxNcohorts, 1)
      for (i_cohort in 1:Ncohorts) {
        obs_dt <- rbind(
          obs_dt,
          data.table(
            siteID = i_site,
            patchID = i_patch,
            cohortID = i_cohort,
            species = sample(1:Nspecies, 1),
            dbh = round(runif(1, 1, 200), 2),
            nTree = sample(1:20, 1)
          )
        )
      }
    }
  }

  # Create an instance of CohortMat
  cohort <- CohortMat$new(obs_df = obs_dt)

  # Convert back to data.table
  obs_dt_restored <- cohort$asDF()

  # Test if the original and restored data.tables are equal
  expect_equal(obs_dt_restored, obs_dt, tolerance = 0.00001)
})
