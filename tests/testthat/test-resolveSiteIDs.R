library(testthat)
library(FINN)
library(data.table)

# resolveSiteIDs() turns name-keyed tables (siteName/patchName/species_name) into
# the integer siteID/patchID/species/0-indexed-year keys FINN simulates on. No
# torch needed.

test_that("contiguous integer site IDs and 0-indexed time", {
  tree <- toy_tree_dt()
  obs  <- suppressMessages(
    makeObsData(tree, plotsize = 1, aggregate_by_site = FALSE)
  )$obs_dt
  env  <- data.table(siteName = "s1", year = c(2000, 2010),
                     temp = c(8.0, 8.2), prec = c(900, 950))

  # createInitCohorts = FALSE keeps this torch-free (we only check the ID tables,
  # not the torch CohortMat), so it runs on CRAN where libtorch is absent.
  res <- suppressMessages(resolveSiteIDs(tree, env, obs, createInitCohorts = FALSE))

  expect_true(all(c("siteID_dt", "tree_dt", "env_dt", "obs_dt",
                    "obs_dt_patches", "species_dt") %in% names(res)))
  # one site -> siteID 1; species coded as contiguous integers 1..N
  expect_equal(sort(unique(res$siteID_dt$siteID)), 1L)
  sp <- res$species_dt[order(species)]
  expect_equal(sp$species, seq_len(nrow(sp)))
  expect_setequal(sp$species_name, c("oak", "pine"))
  # time is remapped to 0-indexed model steps (period - 1)
  expect_equal(min(res$obs_dt$year), 0L)
  expect_true(is.integer(res$obs_dt$siteID))
})

test_that("'other' is always coded as the last species", {
  # capping to 1 species lumps the rest into 'other'; it must sort last so the
  # integer code is stable regardless of alphabetical order.
  mo  <- suppressMessages(
    makeObsData(toy_tree_dt(), plotsize = 1, aggregate_by_site = FALSE, Nspecies = 1)
  )
  env <- data.table(siteName = "s1", year = c(2000, 2010),
                    temp = c(8.0, 8.2), prec = c(900, 950))

  # feed the RECODED tree list (with 'other') back in, as the pipeline does
  res <- suppressMessages(resolveSiteIDs(mo$tree_dt, env, mo$obs_dt, createInitCohorts = FALSE))
  sp  <- res$species_dt[order(species)]
  expect_identical(sp$species_name[nrow(sp)], "other")
})

test_that("env_dt comes back with one row per site and year, however many patches", {
  # Two patches on one site. The site/patch/year map has one row per patch, and
  # the environment used to be merged onto it row by row, so a site with four
  # patches got every climate row four times -- which climateDF2array() now
  # refuses, so fit() failed on the output of the documented pipeline.
  tree <- rbind(toy_tree_dt(),
                toy_tree_dt()[, `:=`(patchName = "p2", treeName = paste0(treeName, "b"))])
  obs  <- suppressMessages(makeObsData(tree, plotsize = 1, aggregate_by_site = FALSE))$obs_dt
  env  <- data.table(siteName = "s1", year = c(2000, 2010),
                     temp = c(8.0, 8.2), prec = c(900, 950))

  res <- suppressMessages(resolveSiteIDs(tree, env, obs, createInitCohorts = FALSE))
  expect_equal(uniqueN(res$obs_dt_patches$patchID), 2L)
  expect_equal(nrow(res$env_dt), 2L)
  expect_false(anyDuplicated(res$env_dt, by = c("siteID", "year")) > 0)
  # and the array builder behind fit() accepts it
  expect_equal(dim(climateDF2array(res$env_dt, c("temp", "prec"))), c(1L, 2L, 2L))
})
