library(testthat)
library(FINN)
library(data.table)

# makeInitCohorts() turns an initial tree list into a CohortMat (the torch-backed
# starting state), so it is torch-gated.

test_that("makeInitCohorts builds a CohortMat that preserves the input cohorts", {
  skip_if_no_torch()

  init <- data.table(
    siteID   = 1L,
    patchID  = 1L,
    species  = c(1L, 2L),
    dbh      = c(20, 30),
    trees    = c(5L, 3L),
    treeName = c("a", "b")
  )

  ic <- makeInitCohorts(init, Nspecies = 2)
  df <- as.data.table(ic$asDF())

  expect_equal(nrow(df), 2L)                       # two cohorts in
  expect_setequal(round(df$dbh, 4), c(20, 30))     # ... two out, dbh preserved
  expect_setequal(df$species, c(1, 2))
  expect_equal(sum(df$trees), 8)                   # 5 + 3 trees
})

test_that("makeInitCohorts bins dbh when dbh_binsize is given", {
  skip_if_no_torch()

  # two stems 2 cm apart collapse into one cohort under a 5 cm bin
  init <- data.table(
    siteID   = 1L, patchID = 1L, species = c(1L, 1L),
    dbh      = c(20, 22), trees = c(4L, 6L), treeName = c("a", "b")
  )
  ic <- makeInitCohorts(init, Nspecies = 1, dbh_binsize = 5)
  df <- as.data.table(ic$asDF())

  expect_equal(nrow(df), 1L)          # merged into a single cohort
  expect_equal(sum(df$trees), 10)     # trees summed across the merged stems
})
