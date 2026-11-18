library(testthat)
library(FINN)
library(data.table)

# fit() and predict() pair data, env and init_cohort by site POSITION. The ids
# need not be 1:N, but the three tables must hold the same sites (#33), and the
# output carries the caller's own siteIDs, mapped back through the lookup the
# model keeps (#34).

fia_subset <- function(sites) {
  obs <- as.data.table(FINN::fia_obs_dt); env <- as.data.table(FINN::fia_env_dt)
  it  <- as.data.table(FINN::fia_init_trees)
  keep <- sort(unique(obs$siteID))[sites]
  list(obs = obs[siteID %in% keep], env = env[siteID %in% keep], it = it[siteID %in% keep],
       keep = keep, nsp = max(obs$species))
}

fia_model <- function(nsp) {
  finn(N_species = nsp, recruits_dbh = 12.9,
       competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
       growth_process       = createProcess(~ temp + prec, FINN::growth),
       regeneration_process = createProcess(~ temp + prec, FINN::regeneration),
       mortality_process    = createProcess(~ temp + prec, FINN::mortality))
}

test_that("fit() stops with the offending ids when data, env and init_cohort hold different sites", {
  skip_if_no_torch()
  d <- fia_subset(1:5)
  four <- d$keep[1:4]
  expect_error(
    fit(fia_model(d$nsp), env = d$env, data = d$obs[siteID %in% four],
        init_cohort = makeInitCohorts(d$it[siteID %in% four], Nspecies = d$nsp),
        epochs = 0L, patches = 4L, patch_size = 0.06, plot_progress = FALSE),
    regexp = sprintf("in `env` but not in `data`: %s.*in `env` but not in `init_cohort`: %s", d$keep[5], d$keep[5]))
  expect_error(
    fit(fia_model(d$nsp), env = d$env[siteID %in% four], data = d$obs[siteID %in% four],
        init_cohort = makeInitCohorts(d$it, Nspecies = d$nsp),
        epochs = 0L, patches = 4L, patch_size = 0.06, plot_progress = FALSE),
    regexp = sprintf("in `init_cohort` but not in `env`: %s", d$keep[5]))
  expect_error(
    fit(fia_model(d$nsp), env = d$env[, siteID := as.character(siteID)], data = d$obs,
        epochs = 0L, patches = 4L, patch_size = 0.06, plot_progress = FALSE),
    regexp = "must be integer-valued.*resolveSiteIDs\\(\\)")
})

test_that("non-contiguous site ids fit, are kept as a lookup, and come back in the output", {
  skip_if_no_torch()
  d <- fia_subset(1:5)
  s <- d$keep[c(1, 2, 4, 5)]
  init <- function() makeInitCohorts(d$it[siteID %in% s], Nspecies = d$nsp)
  expect_identical(init()$site_ids, s)
  m <- fia_model(d$nsp)
  expect_no_warning(fit(m, env = d$env[siteID %in% s], data = d$obs[siteID %in% s], init_cohort = init(),
                        epochs = 0L, patches = 4L, patch_size = 0.06, plot_progress = FALSE))
  expect_equal(m$site_lookup$siteID, s)
  expect_equal(m$site_lookup$index, seq_along(s))
  p <- predict(m, env = d$env[siteID %in% s], init_cohort = init(), patches = 4L, patch_size = 0.06)
  expect_setequal(unique(p$long$site$siteID), s)
  expect_setequal(unique(p$wide$site$siteID), s)
  expect_error(predict(m, env = d$env[siteID %in% s[1:3]], init_cohort = init(), patches = 4L, patch_size = 0.06),
               regexp = sprintf("in `init_cohort` but not in `env`: %s", s[4]))
})

test_that("obsDF2arrays() numbers sites in radix order and records the user ids", {
  obs <- data.frame(siteID = c(10, 3, 3, 7), patchID = 1, cohortID = c(1, 1, 2, 1),
                    species = 1, dbh = c(10, 20, 30, 40), trees = 1)
  expect_message(a <- obsDF2arrays(obs), regexp = "siteIDs are not 1:3")
  expect_identical(attr(a, "site_ids"), c(3, 7, 10))
  expect_no_message(obsDF2arrays(transform(obs, siteID = match(siteID, c(3, 7, 10)))))
  expect_equal(dim(a$dbh)[1], 3L)
  expect_equal(a$dbh[3, 1, 1], 10)     # site 10 is the third position
  expect_equal(a$dbh[1, 1, 1:2], c(20, 30))
})
