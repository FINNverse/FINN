library(testthat)
library(FINN)
library(data.table)

# Multiple shooting (reanchor = TRUE): at every timestep that has an anchor, the
# simulated cohort state is discarded after that timestep's loss has been
# backpropagated and rebuilt from the observed tree list, instead of being carried
# forward. Tests run on the same 4-site FIA subset as test-fit.R: env covers 3
# timesteps and observations exist at timesteps 2 and 3, so year_sequence = c(2, 3)
# and an anchor at timestep 2 is visible in the year-3 state.

reanchor_data = function(n_sites = 4L) {
  ext = function(f) system.file("extdata", f, package = "FINN")
  obs = data.table::fread(ext("fia_obs_dt.csv"))
  env = data.table::fread(ext("fia_env_dt.csv"))
  it  = data.table::fread(ext("fia_init_trees.csv"))
  sites = utils::head(sort(unique(obs$siteID)), n_sites)
  obs = obs[siteID %in% sites]
  list(obs = obs, env = env[siteID %in% sites], it = it[siteID %in% sites],
       Nsp = max(obs$species))
}

# a fresh, identically seeded model + init cohorts, so that two fits differ only
# in the arguments under test
reanchor_model = function(d) {
  FINN.seed(1)
  init = makeInitCohorts(d$it, Nspecies = d$Nsp)
  m = finn(
    N_species            = d$Nsp, recruits_dbh = 12.9,
    competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    growth_process       = createProcess(~ temp + prec, FINN::growth,       optimizeSpecies = TRUE, optimizeEnv = TRUE),
    regeneration_process = createProcess(~ temp + prec, FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = createProcess(~ temp + prec, FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
  )
  list(model = m, init = init)
}

reanchor_fit = function(f, d, epochs = 2L, lr = 0.01, ...) {
  suppressMessages(
    fit(f$model, env = d$env, data = d$obs, init_cohort = f$init, device = "cpu",
        epochs = epochs, patch_size = 0.06, lr = lr, plot_progress = FALSE, ...)
  )
}

# `species` is recycled along the site dimension, so a vector of length `sites`
# gives every site its own species - the marker the batch-subsetting test needs
anchor_mat = function(species, sp, sites = 4L, patches = 4L, ncohorts = 2L,
                      dbh = 25, trees = 50) {
  dims = c(sites, patches, ncohorts)
  CohortMat(dbh = array(dbh, dim = dims), trees = array(trees, dim = dims),
            species = array(as.integer(species), dim = dims), dims = dims, sp = sp)
}

pars = function(m) unlist(lapply(m$parameters, function(p) as.numeric(p$cpu())))

# trees per species at one timestep, summed over the sites of the last batch
species_trees = function(m, timestep) {
  p = m$pred$wide$site[year == timestep]
  p[, .(trees = sum(trees)), by = species][order(species)]
}


test_that("reanchor = FALSE leaves fitting unchanged", {
  skip_if_no_torch()
  d = reanchor_data()
  # The reference check against origin/main cannot live in the suite (it needs a
  # second package tree), so it was run separately: 5 epochs over 2 batches,
  # 198 parameter values, max |diff| = 0, loss history identical. What is checked
  # here is that the new arguments leave the default path alone.
  a = reanchor_model(d)
  reanchor_fit(a, d)
  b = reanchor_model(d)
  reanchor_fit(b, d, reanchor = FALSE)

  expect_identical(pars(a$model), pars(b$model))
  expect_identical(a$model$history, b$model$history)
})


test_that("reanchor = TRUE changes the fitted parameters", {
  skip_if_no_torch()
  d = reanchor_data()
  a = reanchor_model(d)
  reanchor_fit(a, d)
  b = reanchor_model(d)
  reanchor_fit(b, d, reanchor = TRUE,
               anchor_cohorts = list("2" = anchor_mat(2L, sp = d$Nsp)))

  expect_false(isTRUE(all.equal(pars(a$model), pars(b$model))))
  expect_true(all(is.finite(pars(b$model))))
})


test_that("the cohort state is replaced by the anchor, not carried forward", {
  skip_if_no_torch()
  d = reanchor_data()
  free = reanchor_model(d)
  reanchor_fit(free, d, epochs = 1L)
  anch = reanchor_model(d)
  reanchor_fit(anch, d, epochs = 1L, reanchor = TRUE,
               anchor_cohorts = list("2" = anchor_mat(2L, sp = d$Nsp)))

  # only timestep 2 is anchored, not every observation year (timestep 3 is the last
  # one, nothing is simulated after it), so the year-3 state is the anchor evolved by
  # one year: anchored on species 2 only, it must be dominated by species 2
  a3 = species_trees(anch$model, 3)
  expect_gt(a3[species == 2, trees] / sum(a3$trees), 0.9)
  # ... and that is not what the free-running state does
  f3 = species_trees(free$model, 3)
  expect_gt(a3[species == 2, trees], 10 * f3[species == 2, trees])
  # dbh at year 3 must descend from the anchor's 25 cm, one year of growth on
  expect_gt(anch$model$pred$wide$site[year == 3 & species == 2, max(dbh)], 25)
})


test_that("the anchor year is scored before the state is replaced", {
  skip_if_no_torch()
  d = reanchor_data()
  # Observations at timestep 2 only, anchored at timestep 2: the loss is measured
  # from the free-running state, so lr = 0 must give exactly the same loss with and
  # without the anchor. Re-anchoring anywhere before that year's aggregation would
  # score the observation against the state just rebuilt from it.
  d$obs = d$obs[year == 1]
  a = reanchor_model(d)
  reanchor_fit(a, d, epochs = 1L, lr = 0)
  b = reanchor_model(d)
  reanchor_fit(b, d, epochs = 1L, lr = 0, reanchor = TRUE,
               anchor_cohorts = list("2" = anchor_mat(2L, sp = d$Nsp)))

  expect_identical(a$model$history, b$model$history)
})


test_that("fitting does not modify the anchors it was given", {
  skip_if_no_torch()
  d = reanchor_data()
  f = reanchor_model(d)
  A = anchor_mat(2L, sp = d$Nsp)
  before = list(dbh = torch::as_array(A$dbh), trees = torch::as_array(A$trees),
                species = torch::as_array(A$species))
  reanchor_fit(f, d, epochs = 1L, reanchor = TRUE, anchor_cohorts = list("2" = A))

  expect_identical(torch::as_array(A$dbh), before$dbh)
  # the model keeps its own copy: writing to it must not reach the caller's object
  f$model$anchor_cohorts[["2"]]$dbh$mul_(0)
  f$model$anchor_cohorts[["2"]]$trees$mul_(0)
  expect_identical(torch::as_array(A$dbh), before$dbh)
  expect_identical(torch::as_array(A$trees), before$trees)
  expect_identical(torch::as_array(A$species), before$species)
})


test_that("anchors are subset by the batch's sites", {
  skip_if_no_torch()
  d = reanchor_data()
  f = reanchor_model(d)
  # site s is anchored on species s; with batchsize 2 and no shuffling the last
  # batch is sites 3 and 4, and model$pred holds that batch (labelled 1 and 2).
  # Subsetting the anchors with anything other than the batch index shows up as
  # the wrong species (verified by mutating the slicing: it reports 1 and 2).
  reanchor_fit(f, d, epochs = 1L, batchsize = 2, shuffle = FALSE, reanchor = TRUE,
               anchor_cohorts = list("2" = anchor_mat(1:4, sp = d$Nsp)))

  p = f$model$pred$wide$site[year == 3]
  expect_identical(p[siteID == 1][which.max(trees), species], 3L)
  expect_identical(p[siteID == 2][which.max(trees), species], 4L)
})


test_that("invalid anchors are refused", {
  skip_if_no_torch()
  d = reanchor_data()
  fit_with = function(...) reanchor_fit(reanchor_model(d), d, epochs = 1L, reanchor = TRUE, ...)

  expect_error(fit_with(anchor_cohorts = NULL), "requires")
  # an unnamed list would never match a timestep either
  expect_error(fit_with(anchor_cohorts = list(anchor_mat(2L, sp = d$Nsp))), "named")
  # timestep 1 carries no observation, so an anchor there would never fire
  expect_error(fit_with(anchor_cohorts = list("1" = anchor_mat(2L, sp = d$Nsp))), "year_sequence")
  # nor would one that is not on an update_step boundary
  expect_error(fit_with(update_step = 2L, anchor_cohorts = list("3" = anchor_mat(2L, sp = d$Nsp))), "is not a multiple of update_step")
  expect_error(fit_with(anchor_cohorts = list("2" = anchor_mat(2L, sp = d$Nsp - 1L))), "sp")
  expect_error(fit_with(anchor_cohorts = list("2" = anchor_mat(2L, sp = d$Nsp, sites = 3L))), "sites")
  expect_error(fit_with(anchor_cohorts = list("2" = anchor_mat(2L, sp = d$Nsp, patches = 3L))), "patches")

  # an anchor inside another observation's rate-aggregation window: with
  # period_length 2 the observation at timestep 3 aggregates over timesteps 2:3
  dp = d
  dp$obs = data.table::copy(d$obs)[, period_length := 2L]
  expect_error(reanchor_fit(reanchor_model(dp), dp, epochs = 1L, reanchor = TRUE,
                            anchor_cohorts = list("2" = anchor_mat(2L, sp = dp$Nsp))),
               "aggregation window")
})


test_that("an update_step that scores no observation is refused", {
  skip_if_no_torch()
  d = reanchor_data()
  # observations at timesteps 2 and 3; with update_step = 4 no loss ever fires
  expect_error(reanchor_fit(reanchor_model(d), d, epochs = 1L, update_step = 4L),
               "scores no observation")
})


test_that("simulation ignores the anchors a model was fitted with", {
  skip_if_no_torch()
  d = reanchor_data()
  f = reanchor_model(d)
  reanchor_fit(f, d, epochs = 1L, reanchor = TRUE,
               anchor_cohorts = list("2" = anchor_mat(2L, sp = d$Nsp)))

  sim = function(m) {
    FINN.seed(42)
    predict(m, env = d$env, init_cohort = f$init, patch_size = 0.06, device = "cpu")$long$site
  }
  with_anchors = sim(f$model)
  f$model$reanchor = FALSE
  f$model$anchor_cohorts = NULL
  expect_identical(with_anchors, sim(f$model))
  expect_true(any(is.finite(with_anchors$value)))

  FINN.seed(42)
  expect_identical(with_anchors,
                   simulateForest(f$model, env = d$env, init_cohort = f$init,
                                  patch_size = 0.06, device = "cpu")$long$site)
})


# anchor_mode = "merge": an inventory sees nothing below its dbh threshold (FIA: 12.7 cm),
# so an anchor built from it holds no small trees and "replace" deletes every simulated
# cohort below the threshold - FINN's whole regeneration pool - at every census. "merge"
# keeps those and concatenates the full anchor, which double counts nothing because the
# observation is complete at and above the threshold.

# The anchor only fires in fit mode (y is not NULL), and only forward() exposes the cohort
# arrays it acts on, so the merge is driven through forward() directly. Growth and
# mortality are switched off through the intercept column of the env networks, which makes
# the state recorded at the step AFTER the anchor the merged array itself: dbh is carried
# unchanged and `trees_before` holds the trees that entered the step.
merge_setup = function(d, thr = 12.7) {
  FINN.seed(1)
  init = makeInitCohorts(d$it, Nspecies = d$Nsp)
  m = finn(
    N_species            = d$Nsp, recruits_dbh = 1.0,
    competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    growth_process       = createProcess(~ temp + prec, FINN::growth,       optimizeSpecies = TRUE, optimizeEnv = TRUE),
    regeneration_process = createProcess(~ temp + prec, FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = createProcess(~ temp + prec, FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
  )
  suppressMessages(fit(m, env = d$env, data = d$obs, init_cohort = init, device = "cpu",
                       epochs = 1L, patch_size = 0.06, lr = 0, plot_progress = FALSE))
  # column 1 of the design matrix is the intercept, so a constant there is a constant
  # environmental effect on every site and in every year
  const_env = function(nn, value) torch::with_no_grad({
    w = nn$parameters[[1]]
    w$fill_(0)
    w[, 1] = value
  })
  const_env(m$nn_growth, -1e5)   # exp(pred) underflows, so dbh is carried unchanged
  const_env(m$nn_mortality, -1e5)
  const_env(m$nn_regeneration, 1)
  envs = lapply(c("mortality_formula", "growth_formula", "regeneration_formula"),
                function(k) FINN:::extract_env(m[[k]], FINN:::apply_env_scaling(d$env, m$env_scaling)))
  sites = length(unique(d$env$siteID))
  list(model = m, thr = thr, anchor = anchor_mat(2L, sp = d$Nsp, sites = sites),
       run = function(anchor = NULL, ...) {
         FINN.seed(42)
         m$forward(dbh = init$dbh, trees = init$trees, species = init$species, env = envs,
                   y = torch::torch_zeros(c(sites, 1L, d$Nsp, 9L)), year_sequence = integer(0),
                   verbose = FALSE, reanchor = !is.null(anchor),
                   anchor_cohorts = if(is.null(anchor)) NULL else list("2" = anchor), ...)
       })
}


test_that("anchor_mode defaults to \"replace\"", {
  skip_if_no_torch()
  d = reanchor_data()
  A = list("2" = anchor_mat(2L, sp = d$Nsp))
  a = reanchor_model(d)
  reanchor_fit(a, d, reanchor = TRUE, anchor_cohorts = A)
  b = reanchor_model(d)
  reanchor_fit(b, d, reanchor = TRUE, anchor_cohorts = A, anchor_mode = "replace")

  expect_identical(pars(a$model), pars(b$model))
  expect_identical(a$model$history, b$model$history)

  # and fit() does hand the mode down: merging carries cohorts that replacing deletes
  # (this model's recruits enter at 12.9 cm), which has to move the parameters
  cc = reanchor_model(d)
  reanchor_fit(cc, d, reanchor = TRUE, anchor_cohorts = A, anchor_mode = "merge", anchor_min_dbh = 13)
  expect_false(isTRUE(all.equal(pars(a$model), pars(cc$model))))
})


test_that("merging keeps the cohorts below the observation threshold, and nothing else", {
  skip_if_no_torch()
  d = reanchor_data()
  s = merge_setup(d)
  out = s$run(s$anchor, anchor_mode = "merge", anchor_min_dbh = s$thr, return_cohorts = TRUE)
  before = out$Predictions$Cohort$cohortStates[[2]]   # recorded before the anchor fires
  after  = out$Predictions$Cohort$cohortStates[[3]]
  id_before = out$Predictions$Cohort$cohortID[[2]]
  id_after  = out$Predictions$Cohort$cohortID[[3]]
  # trees_before is NA for the cohorts recruited during the step, i.e. everything the
  # anchor did not hand over
  merged = !is.na(after$trees_before)
  sub = before$trees > 0.5 & before$dbh < s$thr
  expect_gt(sum(sub), 10)                             # there is a regeneration pool to lose

  carried = match(id_before[sub], id_after)
  expect_false(anyNA(carried))
  expect_identical(after$trees_before[carried], before$trees[sub])
  expect_identical(after$dbh[carried], before$dbh[sub])
  # every simulated cohort the anchor can see is gone ...
  expect_false(any(id_before[before$trees > 0.5 & before$dbh >= s$thr] %in% id_after[merged]))
  # ... so the basal area above the threshold is the anchor's, exactly
  above = merged & after$dbh >= s$thr
  expect_identical(sum(BA_stem(after$dbh[above])*after$trees_before[above]),
                   sum(BA_stem(torch::as_array(s$anchor$dbh))*torch::as_array(s$anchor$trees)))

  # replacing drops the pool instead
  rep = s$run(s$anchor, return_cohorts = TRUE)$Predictions$Cohort$cohortStates[[3]]
  expect_identical(sum(rep$trees_before[!is.na(rep$trees_before) & rep$dbh < s$thr]), 0)
})


test_that("merging reuses no cohort id", {
  skip_if_no_torch()
  d = reanchor_data()
  s = merge_setup(d)
  out = s$run(s$anchor, anchor_mode = "merge", anchor_min_dbh = s$thr, return_cohorts = TRUE)
  after = out$Predictions$Cohort$cohortStates[[3]]
  id_before = out$Predictions$Cohort$cohortID[[2]]
  id_after  = out$Predictions$Cohort$cohortID[[3]]

  live = after$trees > 0.5
  expect_identical(anyDuplicated(id_after[live]), 0L)
  # the anchor's ids continue above the carried state's, they do not restart
  fresh = id_after[live & !(id_after %in% id_before)]
  expect_gt(length(fresh), 0)
  expect_gt(min(fresh), max(id_before))
})


test_that("the merged cohorts carry no gradient past the anchor", {
  skip_if_no_torch()
  d = reanchor_data()
  s = merge_setup(d)
  # with only the mortality parameter trainable, the growth rate of the step after the
  # anchor can carry a graph only through the state the anchor handed over
  probe = function(...) {
    lapply(s$model$parameters, function(p) p$requires_grad_(FALSE))
    s$model$parameters[["par_mortality_unconstrained"]]$requires_grad_(TRUE)
    s$run(...)
    s$model$g$requires_grad
  }

  expect_true(probe(update_step = 3L))   # free-running and untruncated: the probe sees one
  expect_false(probe(s$anchor, anchor_mode = "merge", anchor_min_dbh = s$thr))
})


test_that("merging without a threshold is refused", {
  skip_if_no_torch()
  d = reanchor_data()
  fit_with = function(...) reanchor_fit(reanchor_model(d), d, epochs = 1L, reanchor = TRUE,
                                        anchor_cohorts = list("2" = anchor_mat(2L, sp = d$Nsp)), ...)

  expect_error(fit_with(anchor_mode = "merge"), "anchor_min_dbh")
  expect_error(fit_with(anchor_mode = "keep"), "should be one of")
})
