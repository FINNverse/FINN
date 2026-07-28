library(testthat)
library(FINN)

# obsDF2arrays() <-> array2obsDF() convert a cohort table to the [site, patch,
# cohort] arrays FINN works on, and back. Rcpp-backed but torch-free, so these
# run on CRAN.

test_that("obsDF2arrays -> array2obsDF round-trips a full cohort grid", {
  # a rectangular grid (every site x patch has the same cohorts) round-trips
  # cleanly; padding only appears for ragged inputs.
  obs <- data.frame(
    siteID   = c(1L, 1L, 1L, 1L),
    patchID  = c(1L, 1L, 2L, 2L),
    cohortID = c(1L, 2L, 1L, 2L),
    species  = c(1L, 2L, 1L, 2L),
    dbh      = c(20, 30, 25, 35),
    trees    = c(5L, 3L, 4L, 2L)
  )

  arr  <- obsDF2arrays(obs)
  expect_named(arr, c("species", "dbh", "trees"), ignore.order = TRUE)
  expect_equal(dim(arr$dbh), c(1, 2, 2))   # 1 site, 2 patches, 2 cohorts

  back <- as.data.frame(array2obsDF(arr))

  key <- function(d) d[order(d$siteID, d$patchID, d$dbh),
                       c("siteID", "patchID", "species", "dbh", "trees")]
  a <- key(obs); b <- key(back)
  rownames(a) <- rownames(b) <- NULL
  expect_equal(b$dbh,     a$dbh,     tolerance = 1e-6)
  expect_equal(b$trees,   a$trees)
  expect_equal(b$species, a$species)
})
