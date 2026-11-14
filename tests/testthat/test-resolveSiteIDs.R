library(testthat)
library(FINN)
library(data.table)

# resolveSiteIDs() turns name-keyed tables (siteName/patchName/species_name) into
# the integer siteID/patchID/species/year keys FINN simulates on. No torch
# needed.

# one row per calendar year 2000..2010 for the toy site(s): the annual axis needs
# the years after the first inventory
toy_env <- function(sites = "s1") CJ(siteName = sites, year = 2000:2010)[, `:=`(temp = 8, prec = 900)][]

test_that("annual axis: years since the first inventory, env 1..Tmax, no year-0 rows", {
  tree <- toy_tree_dt()
  mo   <- suppressMessages(makeObsData(tree, plotsize = 1, aggregate_by_site = FALSE))
  res  <- suppressMessages(resolveSiteIDs(mo$tree_dt, toy_env(), mo$obs_dt, createInitCohorts = FALSE))

  expect_equal(unique(res$obs_dt$year), 10)               # the 2010 inventory, 10 years on
  expect_equal(unique(res$obs_dt$period_length), 10)
  expect_equal(res$env_dt$year, 1:10)                     # no row 0
  # fit()'s env: siteID, year and the predictors only
  expect_identical(names(res$env_dt), c("siteID", "year", "temp", "prec"))
  expect_setequal(res$tree_dt$year, c(0, 10))
  # the yearly mortality rate survives the aggregation to the site
  expect_equal(res$obs_dt[species_name == "oak", mort], 1 - 0.5^(1/10), tolerance = 1e-10)
  # and the env is what climateDF2array() / fit() take
  expect_equal(dim(climateDF2array(res$env_dt, c("temp", "prec"))), c(1L, 10L, 2L))

  # a missing year of env is refused with a message that says which
  expect_error(resolveSiteIDs(mo$tree_dt, toy_env()[year != 2005], mo$obs_dt, createInitCohorts = FALSE),
               "no row for 1 site-year")
})

test_that("time_axis = \"census\" keeps the inventory index of FINN 0.1.0", {
  tree <- toy_tree_dt()
  mo   <- suppressMessages(makeObsData(tree, plotsize = 1, aggregate_by_site = FALSE, rates = "interval"))
  env  <- data.table(siteName = "s1", year = c(2000, 2010), temp = c(8.0, 8.2), prec = c(900, 950))
  res  <- suppressMessages(resolveSiteIDs(mo$tree_dt, env, mo$obs_dt, createInitCohorts = FALSE,
                                          time_axis = "census"))
  expect_setequal(res$obs_dt$year, c(0, 1))
  expect_equal(res$env_dt$year, c(0, 1))
  expect_equal(res$obs_dt[species_name == "oak" & year == 1, mort], 0.5)
  expect_equal(unique(res$obs_dt[year == 1, period_length]), 1)
  # yearly rates on the census axis are a mismatch, and said so
  ma <- suppressMessages(makeObsData(tree, plotsize = 1, aggregate_by_site = FALSE))
  expect_warning(resolveSiteIDs(ma$tree_dt, env, ma$obs_dt, createInitCohorts = FALSE, time_axis = "census"),
                 "rates = \"interval\"")
})

test_that("contiguous integer site IDs and 0-indexed time", {
  tree <- toy_tree_dt()
  obs  <- suppressMessages(
    makeObsData(tree, plotsize = 1, aggregate_by_site = FALSE)
  )$obs_dt
  env  <- toy_env()

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
  # time is remapped to years since the first inventory, which is the initial
  # state (year 0 of the tree list) and not an observation
  expect_equal(min(res$tree_dt$year), 0L)
  expect_equal(min(res$obs_dt$year), 10L)
  expect_true(is.integer(res$obs_dt$siteID))
})

test_that("'other' is always coded as the last species", {
  # capping to 1 species lumps the rest into 'other'; it must sort last so the
  # integer code is stable regardless of alphabetical order.
  mo  <- suppressMessages(
    makeObsData(toy_tree_dt(), plotsize = 1, aggregate_by_site = FALSE, Nspecies = 1)
  )
  env <- toy_env()

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
  env  <- toy_env()

  res <- suppressMessages(resolveSiteIDs(tree, env, obs, createInitCohorts = FALSE))
  expect_equal(uniqueN(res$obs_dt_patches$patchID), 2L)
  expect_equal(nrow(res$env_dt), 10L)
  expect_false(anyDuplicated(res$env_dt, by = c("siteID", "year")) > 0)
  # and the array builder behind fit() accepts it
  expect_equal(dim(climateDF2array(res$env_dt, c("temp", "prec"))), c(1L, 10L, 2L))

  # the same on the census axis, with the inventory years only
  env2 <- data.table(siteName = "s1", year = c(2000, 2010), temp = c(8.0, 8.2), prec = c(900, 950))
  obs2 <- suppressMessages(makeObsData(tree, plotsize = 1, aggregate_by_site = FALSE, rates = "interval"))$obs_dt
  res2 <- suppressMessages(resolveSiteIDs(tree, env2, obs2, createInitCohorts = FALSE, time_axis = "census"))
  expect_equal(nrow(res2$env_dt), 2L)
  expect_equal(dim(climateDF2array(res2$env_dt, c("temp", "prec"))), c(1L, 2L, 2L))
})
