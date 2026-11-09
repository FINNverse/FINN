## Cohort identity: an initial cohort keeps the id makeInitCohorts() gave it,
## and an id is never re-issued to a different cohort.
##
## This is what lets a simulated cohort be traced back to the inventory tree it
## came from: makeInitCohorts() without `dbh_binsize` does not bin, so one row
## of the initial-tree table becomes one cohort, and the join key back to that
## row is (siteID, patchID, cohortID).
skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)

Nsp <- 3L; Nsites <- 2L; Npatch <- 3L; Tmax <- 12L

mk_model <- function() finn(N_species = Nsp, recruits_dbh = 5,
  competition_process  = createProcess(~0, func = FINN::competition),
  growth_process       = createProcess(~1 + env1, initEnv = matrix(c(rep(0.5, Nsp), rep(0.2, Nsp)), Nsp, 2),
                                       initSpecies = matrix(c(rep(0.1, Nsp), rep(0.05, Nsp)), Nsp, 2), func = FINN::growth),
  mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(rep(-2.5, Nsp), rep(0.2, Nsp)), Nsp, 2),
                                       initSpecies = matrix(c(rep(0.2, Nsp), rep(0.1, Nsp), rep(0, Nsp)), Nsp, 3), func = FINN::mortality),
  regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(rep(1.0, Nsp), rep(0.2, Nsp)), Nsp, 2),
                                       initSpecies = rep(0.1, Nsp), func = FINN::regeneration,
                                       sample_regeneration = FALSE))

## one row per tree: the case the identity claim is about
init_trees <- CJ(siteID = 1:Nsites, patchID = 1:Npatch, k = 1:4)[
  , .(siteID, patchID, species = rep_len(1:Nsp, .N), dbh = rep_len(c(12, 20, 28, 36), .N),
      treeName = paste0("t", .I), trees = 1)]
env <- data.table(expand.grid(siteID = 1:Nsites, year = 1:Tmax))[, env1 := (siteID - 1.5)]

test_that("an initial cohort keeps the id makeInitCohorts() gave it", {
  ic  <- makeInitCohorts(copy(init_trees), Nspecies = Nsp, treeID_table = TRUE)
  key <- as.data.table(ic$init_trees)[, .(siteID, patchID, cohortID, treeName, sp0 = species, dbh0 = dbh)]
  expect_equal(nrow(key), nrow(init_trees))          # no binning: one cohort per tree

  FINN.seed(3)
  sim <- simulateForest(mk_model(), env = env, init_cohort = ic$initCohort, patches = Npatch,
                        patch_size = 0.1, device = "cpu", return_cohorts = 1L)
  C <- as.data.table(sim$wide$cohort)
  for (cc in c("siteID", "patchID", "species", "cohortID")) C[, (cc) := as.integer(as.character(get(cc)))]

  ## every initial tree is there, under its own id, with its own species
  j <- merge(C, key, by = c("siteID", "patchID", "cohortID"))
  expect_equal(nrow(j), nrow(key))
  expect_true(all(j$species == j$sp0))
  ## after one year of growth a LIVING cohort's diameter is its input diameter
  ## plus a year of growth. Dead cohorts are recorded with dbh 0 (the recorded
  ## state is masked by trees > 0.5), so they are excluded here - worth knowing
  ## before reading a traced cohort's diameter as a prediction.
  ## (one year of growth is species- and light-dependent, so this is a per-cohort
  ## bound rather than an exact value: the point is that the row we joined to is
  ## THAT tree and not another one.)
  alive <- j[trees > 0.5]
  expect_gt(nrow(alive), 0.5 * nrow(j))
  expect_true(all(alive$dbh >= alive$dbh0 - 1e-6))
  expect_gt(cor(alive$dbh, alive$dbh0), 0.95)
})

test_that("a cohort id is never re-issued to a different cohort", {
  ic <- makeInitCohorts(copy(init_trees), Nspecies = Nsp)
  FINN.seed(4)
  sim <- simulateForest(mk_model(), env = env, init_cohort = ic, patches = Npatch,
                        patch_size = 0.1, device = "cpu", return_cohorts = TRUE)
  C <- as.data.table(sim$wide$cohort)
  for (cc in c("siteID", "patchID", "species", "cohortID", "year")) C[, (cc) := as.integer(as.character(get(cc)))]
  ## an id identifies ONE cohort within a patch, so its species can never change
  spp <- C[, .(n_species = uniqueN(species)), by = .(siteID, patchID, cohortID)]
  expect_true(all(spp$n_species == 1L))
  ## and once an id has left the table it must not come back
  gaps <- C[, .(reappears = any(diff(sort(unique(year))) > 1L)), by = .(siteID, patchID, cohortID)]
  expect_false(any(gaps$reappears))
})
