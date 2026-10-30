library(testthat)
library(FINN)
library(data.table)

# Teacher forcing (teacher_forcing = TRUE): no trajectory is ever simulated. The cohort
# state is reset at EVERY timestep to the last observed stand - init_cohort for the first
# interval, then anchor_cohorts - the demographic processes are evaluated there, and only
# the three RATE responses are scored. Tests share test-reanchor.R's 4-site FIA subset:
# env covers 3 timesteps, observations at timesteps 2 and 3, so year_sequence = c(2, 3),
# the last observation is timestep 3 and the only anchor teacher forcing needs is at 2.

tf_data = function(n_sites = 4L) {
  ext = function(f) system.file("extdata", f, package = "FINN")
  obs = data.table::fread(ext("fia_obs_dt.csv"))
  env = data.table::fread(ext("fia_env_dt.csv"))
  it  = data.table::fread(ext("fia_init_trees.csv"))
  sites = utils::head(sort(unique(obs$siteID)), n_sites)
  obs = obs[siteID %in% sites]
  list(obs = obs, env = env[siteID %in% sites], it = it[siteID %in% sites],
       Nsp = max(obs$species))
}

tf_model = function(d) {
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

tf_fit = function(f, d, epochs = 2L, lr = 0.01, ...) {
  suppressMessages(
    fit(f$model, env = d$env, data = d$obs, init_cohort = f$init, device = "cpu",
        epochs = epochs, patch_size = 0.06, lr = lr, plot_progress = FALSE, ...)
  )
}

tf_anchor = function(species, sp, sites = 4L, patches = 4L, ncohorts = 2L,
                     dbh = 25, trees = 50) {
  dims = c(sites, patches, ncohorts)
  CohortMat(dbh = array(dbh, dim = dims), trees = array(trees, dim = dims),
            species = array(as.integer(species), dim = dims), dims = dims, sp = sp)
}

tf_pars = function(m) unlist(lapply(m$parameters, function(p) as.numeric(p$cpu())))


test_that("teacher_forcing = FALSE leaves fitting unchanged", {
  skip_if_no_torch()
  d = tf_data()
  a = tf_model(d)
  tf_fit(a, d)
  b = tf_model(d)
  tf_fit(b, d, teacher_forcing = FALSE)

  expect_identical(tf_pars(a$model), tf_pars(b$model))
  expect_identical(a$model$history, b$model$history)
})


test_that("teacher_forcing = TRUE changes the fitted parameters", {
  skip_if_no_torch()
  d = tf_data()
  a = tf_model(d)
  tf_fit(a, d)
  b = tf_model(d)
  tf_fit(b, d, teacher_forcing = TRUE,
         anchor_cohorts = list("2" = tf_anchor(2L, sp = d$Nsp)))

  expect_false(isTRUE(all.equal(tf_pars(a$model), tf_pars(b$model))))
  expect_true(all(is.finite(tf_pars(b$model))))
})


# The reset only fires in fit mode (y is not NULL), and only forward() exposes the cohort
# arrays it acts on, so it is driven through forward() directly. Growth is switched off
# through the intercept column of the env networks, so the state entering a timestep is
# exactly the state that was forced into it: `trees_before` holds the trees that entered
# the step and is NA for the cohorts recruited during it.
tf_setup = function(d) {
  f = tf_model(d)
  suppressMessages(fit(f$model, env = d$env, data = d$obs, init_cohort = f$init, device = "cpu",
                       epochs = 1L, patch_size = 0.06, lr = 0, plot_progress = FALSE))
  m = f$model
  const_env = function(nn, value) torch::with_no_grad({
    w = nn$parameters[[1]]
    w$fill_(0)
    w[, 1] = value
  })
  const_env(m$nn_growth, -1e5)      # exp(pred) underflows, so dbh is carried unchanged
  const_env(m$nn_regeneration, 3)   # ... but recruits keep arriving, so a free run accumulates
  envs = lapply(c("mortality_formula", "growth_formula", "regeneration_formula"),
                function(k) FINN:::extract_env(m[[k]], FINN:::apply_env_scaling(d$env, m$env_scaling)))
  sites = length(unique(d$env$siteID))
  list(model = m, init = f$init, sites = sites,
       run = function(...) {
         FINN.seed(42)
         m$forward(dbh = f$init$dbh, trees = f$init$trees, species = f$init$species, env = envs,
                   y = torch::torch_zeros(c(sites, 1L, d$Nsp, 9L)), year_sequence = integer(0),
                   verbose = FALSE, return_cohorts = TRUE, ...)
       })
}

# Cohort ids are issued above the running maximum, so the highest id alive at a timestep
# only rises while the state is CARRIED. Resetting the state puts back the ids it was built
# with, which makes the peak id an exact, deterministic signature of what each step started
# from - unlike the stand itself, which mortality and recruitment resample every year.
peak_id = function(out, step) max(out$Predictions$Cohort$cohortID[[step]])
tot_trees = function(out, step) sum(out$Predictions$Site$trees[, step, ])


test_that("the state entering every timestep is the last observed stand", {
  skip_if_no_torch()
  d = tf_data()
  s = tf_setup(d)
  A = tf_anchor(2L, sp = d$Nsp, sites = s$sites)

  # no anchors: every timestep restarts from init_cohort, so no id is ever issued twice
  # over and the peak stands still
  tf = s$run(teacher_forcing = TRUE)
  expect_identical(peak_id(tf, 2), peak_id(tf, 1))
  expect_identical(peak_id(tf, 3), peak_id(tf, 1))
  # ... and the stand does not accumulate three years of recruits either
  expect_lt(tot_trees(tf, 3) / tot_trees(tf, 1), 1.35)  # only the year env differs

  # a free run carries its state instead, so both rise
  fr = s$run()
  expect_gt(peak_id(fr, 2), peak_id(fr, 1))
  expect_gt(peak_id(fr, 3), peak_id(fr, 2))
  expect_gt(tot_trees(fr, 3) / tot_trees(fr, 1), 1.5)

  # the anchor at timestep 2 takes over from timestep 3 on, and not before
  tfa = s$run(teacher_forcing = TRUE, anchor_cohorts = list("2" = A))
  expect_identical(peak_id(tfa, 2), peak_id(tfa, 1))
  expect_gt(peak_id(tfa, 3), peak_id(tfa, 2))
  # the anchor is 50 stems per slot of species 2 only, so it must dominate timestep 3 and
  # be absent from timestep 2
  sp3 = colSums(tfa$Predictions$Site$trees[, 3, ])
  sp2 = colSums(tfa$Predictions$Site$trees[, 2, ])
  expect_gt(sp3[2] / sum(sp3), 0.85)
  expect_lt(sp2[2] / sum(sp2), 0.3)

  # multiple shooting free-runs inside the segment and only rebuilds AT the census
  ms = s$run(reanchor = TRUE, anchor_cohorts = list("2" = A))
  expect_identical(tot_trees(ms, 2), tot_trees(fr, 2))
  expect_false(isTRUE(all.equal(tot_trees(ms, 2), tot_trees(tfa, 2))))
})


test_that("the state responses leave the loss and the rate responses stay", {
  skip_if_no_torch()
  d = tf_data()
  f = tf_model(d)
  tf_fit(f, d, epochs = 1L, lr = 0, teacher_forcing = TRUE,
         anchor_cohorts = list("2" = tf_anchor(2L, sp = d$Nsp)))
  h = as.numeric(f$model$history[[1]])

  expect_identical(h[1:3], c(0, 0, 0))     # dbh, ba, trees
  expect_true(all(h[4:6] > 0))             # growth, mortality, regeneration
  expect_true(all(is.finite(h)))
  # the default path scores all six, so the zeros above are the treatment and not the data
  g = tf_model(d)
  tf_fit(g, d, epochs = 1L, lr = 0)
  expect_true(all(as.numeric(g$model$history[[1]])[1:3] > 0))
})


test_that("update_step is inert under teacher forcing", {
  skip_if_no_torch()
  d = tf_data()
  # There is no multi-year trajectory to truncate, so where the gradient is cut cannot
  # matter. Restricted to the timestep-2 observation, both update_step 1 and 2 score it,
  # and everything is then conditioned on init_cohort - so the two fits must be identical.
  ds = d
  ds$obs = d$obs[year == 1]
  a = tf_model(ds)
  tf_fit(a, ds, teacher_forcing = TRUE, update_step = 1L)
  b = tf_model(ds)
  tf_fit(b, ds, teacher_forcing = TRUE, update_step = 2L)

  expect_identical(tf_pars(a$model), tf_pars(b$model))
  expect_identical(a$model$history, b$model$history)
  # ... and the same two values are NOT interchangeable when a trajectory is simulated
  fa = tf_model(ds)
  tf_fit(fa, ds, update_step = 1L)
  fb = tf_model(ds)
  tf_fit(fb, ds, update_step = 2L)
  expect_false(isTRUE(all.equal(tf_pars(fa$model), tf_pars(fb$model))))

  # an update_step that would silently drop an observation is refused instead
  expect_error(tf_fit(tf_model(d), d, epochs = 1L, teacher_forcing = TRUE, update_step = 2L,
                      anchor_cohorts = list("2" = tf_anchor(2L, sp = d$Nsp))),
               "drops the observations")
})


test_that("teacher forcing and re-anchoring are mutually exclusive", {
  skip_if_no_torch()
  d = tf_data()
  expect_error(tf_fit(tf_model(d), d, epochs = 1L, reanchor = TRUE, teacher_forcing = TRUE,
                      anchor_cohorts = list("2" = tf_anchor(2L, sp = d$Nsp))),
               "cannot be combined")
})


test_that("a missing anchor is refused, and the last observation needs none", {
  skip_if_no_torch()
  d = tf_data()
  # year_sequence is c(2, 3): timestep 2 must be anchored, timestep 3 must not be, and
  # the first interval is conditioned on init_cohort.
  expect_error(tf_fit(tf_model(d), d, epochs = 1L, teacher_forcing = TRUE,
                      anchor_cohorts = list("3" = tf_anchor(2L, sp = d$Nsp))),
               "missing 2")
  expect_error(tf_fit(tf_model(d), d, epochs = 1L, teacher_forcing = TRUE), "missing 2")
  f = tf_model(d)
  tf_fit(f, d, epochs = 1L, teacher_forcing = TRUE,
         anchor_cohorts = list("2" = tf_anchor(2L, sp = d$Nsp)))
  expect_true(all(is.finite(tf_pars(f$model))))
})


test_that("anchors are subset by the batch's sites under teacher forcing", {
  skip_if_no_torch()
  d = tf_data()
  f = tf_model(d)
  # site s anchored on species s; batchsize 2 without shuffling leaves sites 3 and 4 in
  # model$pred, labelled 1 and 2
  tf_fit(f, d, epochs = 1L, batchsize = 2, shuffle = FALSE, teacher_forcing = TRUE,
         anchor_cohorts = list("2" = tf_anchor(1:4, sp = d$Nsp)))

  p = f$model$pred$wide$site[year == 3]
  expect_identical(p[siteID == 1][which.max(trees), species], 3L)
  expect_identical(p[siteID == 2][which.max(trees), species], 4L)
})


test_that("simulation ignores the teacher forcing a model was fitted with", {
  skip_if_no_torch()
  d = tf_data()
  f = tf_model(d)
  tf_fit(f, d, epochs = 1L, teacher_forcing = TRUE,
         anchor_cohorts = list("2" = tf_anchor(2L, sp = d$Nsp)))

  sim = function(m) {
    FINN.seed(42)
    predict(m, env = d$env, init_cohort = f$init, patch_size = 0.06, device = "cpu")$long$site
  }
  forced = sim(f$model)
  f$model$teacher_forcing = FALSE
  f$model$anchor_cohorts = NULL
  expect_identical(forced, sim(f$model))
  expect_true(any(is.finite(forced$value)))
})
