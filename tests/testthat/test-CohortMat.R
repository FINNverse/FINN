# tests/testthat/test-CohortMat.R
library(testthat)
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
  obs_dt <- data.frame()
  for (i_site in 1:Nsites) {
    for (i_patch in 1:Npatches) {
      Ncohorts <- sample(1:maxNcohorts, 1)
      for (i_cohort in 1:Ncohorts) {
        obs_dt <- rbind(
          obs_dt,
          data.frame(
            siteID = i_site,
            patchID = i_patch,
            cohortID = i_cohort,
            species = sample(1:Nspecies, 1),
            dbh = round(runif(1, 1, 200), 2),
            trees = sample(1:20, 1)
          )
        )
      }
    }
  }

  # Create an instance of CohortMat
  cohort <- CohortMat$new(obs_df = obs_dt)

  # Convert back to data.frame
  obs_dt_restored <- cohort$asDF()

  # Test if the original and restored data.frames are equal
  expect_equal(obs_dt_restored, obs_dt, tolerance = 0.00001)
})

test_that("rweibull_cohorts creates valid data frame that can be use to create a cohort tensor", {
  #examples
  obs_df =
    rweibull_cohorts(
      trees = c(300,10),
      dbh_shape = c(3,3),
      dbh_scale = c(1,50),
      dbh_class_range = 0.1,
      patchID = c(1,1),
      species = c(3,4)
    )
  cohortMat1 = CohortMat$new(obs_df = obs_df)
  testthat::expect_equal(cohortMat1$asDF(),obs_df, tolerance = 0.00001)
})
