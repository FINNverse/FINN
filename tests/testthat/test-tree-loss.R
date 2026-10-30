library(testthat)
library(FINN)
library(data.table)

# `loss_resolution`: growth and mortality scored per individual tree INSTEAD of per
# (site, year, species) cell, keyed by the cohort ids. Same 4-site FIA subset as
# test-reanchor.R; `period_length` is explicit so the window is known and
# `shuffle = FALSE` keeps the batch order equal to the site order the targets are keyed in.

tl_data = function(years_obs = c(1, 2), period = 1L, n_sites = 4L) {
  ext = function(f) system.file("extdata", f, package = "FINN")
  obs = data.table::fread(ext("fia_obs_dt.csv"))
  env = data.table::fread(ext("fia_env_dt.csv"))
  it  = data.table::fread(ext("fia_init_trees.csv"))
  sites = utils::head(sort(unique(obs$siteID)), n_sites)
  obs = obs[siteID %in% sites & year %in% years_obs]
  obs[, period_length := period]
  list(obs = obs, env = env[siteID %in% sites], it = it[siteID %in% sites],
       Nsp = max(obs$species))
}

tl_model = function(d) {
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

tl_fit = function(f, d, epochs = 1L, lr = 0, ...) {
  suppressMessages(
    fit(f$model, env = d$env, data = d$obs, init_cohort = f$init, device = "cpu",
        epochs = epochs, patch_size = 0.06, lr = lr, plot_progress = FALSE,
        shuffle = FALSE, ...)
  )
}

# a per-tree target of the right shape; the values only have to be observable, the
# identity tests read the predictions and not the residual
tl_target = function(dims) list(growth = array(0.01, dim = dims),
                                mort = array(rep(c(0, 1), length.out = prod(dims)), dim = dims))

tl_anchor = function(sp, dims = c(4L, 4L, 2L)) {
  CohortMat(dbh = array(25, dim = dims), trees = array(1, dim = dims),
            species = array(2L, dim = dims), dims = dims, sp = sp)
}

tl_pars = function(m) unlist(lapply(m$parameters, function(p) as.numeric(p$cpu())))


test_that("loss_resolution = \"site\" leaves fitting unchanged", {
  skip_if_no_torch()
  # The reference check against the branch head needs a second package tree and was run
  # separately: 5 epochs x 5 configurations (auto and numeric weights, update_step 1 and 3,
  # and all three aggregation branches), 198 parameters, identical() on parameters, history,
  # weights and baselines with the arguments absent and with loss_resolution = "site".
  d = tl_data()
  a = tl_model(d)
  tl_fit(a, d, epochs = 2L, lr = 0.01)
  b = tl_model(d)
  tl_fit(b, d, epochs = 2L, lr = 0.01, loss_resolution = "site")

  expect_identical(tl_pars(a$model), tl_pars(b$model))
  expect_identical(a$model$history, b$model$history)
  expect_length(a$model$history[[1]], 7L)
})


test_that("loss_resolution = \"tree\" moves two terms and leaves the other four alone", {
  skip_if_no_torch()
  d = tl_data()
  dims = dim(tl_model(d)$init$dbh_r)
  a = tl_model(d)
  tl_fit(a, d, epochs = 2L, lr = 0.01)
  b = tl_model(d)
  tl_fit(b, d, epochs = 2L, lr = 0.01, loss_resolution = "tree",
         tree_data = list("2" = tl_target(dims), "3" = tl_target(dims)))

  # the loss vector keeps its width and its positions: the same response, another resolution
  expect_length(b$model$history[[1]], 7L)
  expect_named(b$model$loss_weights,
               c("dbh", "ba", "trees", "growth", "mortality", "regeneration"))
  expect_true(all(is.finite(tl_pars(b$model))))
  expect_false(isTRUE(all.equal(tl_pars(a$model), tl_pars(b$model))))
  expect_true(all(b$model$history[[1]][4:5] > 0))
  # growth and mortality are scaled by a baseline computed over the TREES, the other four
  # by the same cell baselines as at site resolution
  expect_identical(a$model$loss_baseline[c("dbh", "ba", "trees", "regeneration")],
                   b$model$loss_baseline[c("dbh", "ba", "trees", "regeneration")])
  expect_false(isTRUE(all.equal(unname(a$model$loss_baseline[c("growth", "mortality")]),
                                unname(b$model$loss_baseline[c("growth", "mortality")]))))
})


test_that("the per-tree and the aggregate predictions agree over a one-year window", {
  skip_if_no_torch()
  # One census, a one-year window, and the cohorts in it are exactly the initial ones. The
  # aggregate growth/mortality of a (site, species) cell is then exactly the tree-weighted
  # mean of that cell's per-tree predictions - the invariant that says a residual belongs
  # to the tree the model thinks it does.
  d = tl_data(years_obs = 1, period = 1L)
  d$env = d$env[year >= 1]          # the observed year becomes the first simulated year
  f = tl_model(d)
  dims = dim(f$init$dbh_r)
  tl_fit(f, d, loss_resolution = "tree", tree_data = list("1" = tl_target(dims)))

  tp = f$model$tree_pred
  expect_length(tp$growth, prod(dims))
  g_slot = array(tp$growth, dim = dims)
  m_slot = array(tp$mort, dim = dims)
  site = f$model$pred$wide$site[year == 1]
  for (s in seq_len(dims[1])) for (k in seq_len(d$Nsp)) {
    w = (f$init$species_r[s,,] == k) * f$init$trees_r[s,,]
    if (sum(w) < 0.5) next
    cell = site[siteID == s & species == k]
    expect_equal(sum(g_slot[s,,]*w)/sum(w), cell$growth, tolerance = 1e-5)
    expect_equal(sum(m_slot[s,,]*w)/sum(w), cell$mort, tolerance = 1e-5)
  }
  # liveness is `trees > 0.5`, never a dbh test: a dead cohort is a zeroed slot that keeps
  # a diameter, so a dead-at-year-0 tree must contribute no year to any accumulator
  expect_true(all(array(tp$n, dim = dims)[f$init$trees_r < 0.5] == 0))
  expect_true(all(array(tp$n, dim = dims)[f$init$trees_r > 0.5] == 1))
})


test_that("no recruit cohort enters the per-tree residual", {
  skip_if_no_torch()
  # Free running with regeneration on: by the census the cohort array carries recruits the
  # initial tree list knows nothing about. Their ids sit above the generation's block, so
  # the residual keeps exactly one slot per initial tree and cannot reach them.
  d = tl_data(years_obs = 2, period = 1L)
  f = tl_model(d)
  dims = dim(f$init$dbh_r)
  tl_fit(f, d, loss_resolution = "tree", tree_data = list("3" = tl_target(dims)))

  cohorts = predict(f$model, env = d$env, init_cohort = f$init, patch_size = 0.06,
                    return_cohorts = "last")$wide$cohort
  expect_gt(nrow(cohorts[dbh < 13]), 0)          # recruits exist in the simulated state
  expect_length(f$model$tree_pred$growth, prod(dims))
  expect_true(all(f$model$tree_pred$n <= 1))     # one year per slot, never two cohorts

  # the model's own per-tree predictions fed back as the target drive growth to zero and
  # mortality to the predictions' entropy over exactly the tracked slots; re-keying the
  # same values to other slots does not. Both need every residual paired with its slot,
  # and a growth term of 0 is also what says the CELL term no longer sits in position 4.
  own = list(growth = array(f$model$tree_pred$growth, dim = dims),
             mort   = array(f$model$tree_pred$mort, dim = dims))
  score = function(td) {
    g = tl_model(d)
    tl_fit(g, d, loss_resolution = "tree", weights = rep(1, 6), tree_data = list("3" = td))
    as.numeric(g$model$history[[1]])[4:5]
  }
  right = score(own)
  off = score(lapply(own, function(a) aperm(a, c(2, 1, 3))))
  expect_lt(right[1], 1e-8)
  expect_gt(off[1], 1e-4)

  y = f$model$tree_pred$mort[f$model$tree_pred$n > 0.5]
  p = pmin(pmax(y, 1e-6), 1 - 1e-6)
  expect_equal(right[2], mean(-(y*log(p) + (1 - y)*log(1 - p))), tolerance = 1e-5)
})


test_that("the interval prediction compounds the annual rates", {
  skip_if_no_torch()
  # Mortality is the compounded per-step probability 1 - prod(1 - m), not the sampled
  # fate; growth follows the aggregate convention. Three lr = 0 fits share one trajectory,
  # so the two-year window must be the exact composition of the two one-year ones.
  runs = lapply(list(list(1, 1L, "2"), list(2, 1L, "3"), list(2, 2L, "3")), function(a) {
    d = tl_data(years_obs = a[[1]], period = a[[2]])
    f = tl_model(d)
    tl_fit(f, d, loss_resolution = "tree", tree_data = stats::setNames(list(tl_target(dim(f$init$dbh_r))), a[[3]]))
    f$model
  })
  y3 = lapply(runs[2:3], function(m) m$pred$wide$site[year == 3][order(siteID, species), mort])
  expect_equal(y3[[1]], y3[[2]], tolerance = 1e-6)   # same trajectory in both

  p = lapply(runs, function(m) m$tree_pred)
  ok = p[[1]]$n > 0.5 & p[[2]]$n > 0.5 & p[[3]]$n > 1.5
  expect_gt(sum(ok), 20)
  expect_equal(p[[3]]$mort[ok], 1 - (1 - p[[1]]$mort[ok])*(1 - p[[2]]$mort[ok]), tolerance = 1e-5)
  expect_equal(p[[3]]$growth[ok], (p[[1]]$growth[ok] + p[[2]]$growth[ok])/2, tolerance = 1e-5)
})


test_that("a new generation of tree identities starts at every anchor", {
  skip_if_no_torch()
  # The state rebuilt from an observed tree list starts a new generation, so the interval
  # after an anchor is scored against THAT list. Under "merge" the carried sub-threshold
  # cohorts keep their old ids and stay out of the residual - the anchor never saw them.
  d = tl_data(years_obs = c(1, 2), period = 1L)
  dims = dim(tl_model(d)$init$dbh_r)
  A = tl_anchor(d$Nsp)
  td = list("2" = tl_target(dims), "3" = tl_target(c(4, 4, 2)))

  for (arg in list(list(reanchor = TRUE),
                   list(reanchor = TRUE, anchor_mode = "merge", anchor_min_dbh = 12.7),
                   list(teacher_forcing = TRUE))) {
    f = tl_model(d)
    do.call(tl_fit, c(list(f, d, loss_resolution = "tree", tree_data = td,
                           anchor_cohorts = list("2" = A)), arg))
    expect_length(f$model$tree_pred$growth, 4*4*2)
    expect_true(all(f$model$tree_pred$n == 1))   # every anchored tree is tracked
    expect_true(all(is.finite(as.numeric(f$model$history[[1]])[4:5])))
    # the interval after the anchor starts from a stand in which every tree is alive, so
    # no window is truncated there (the year-0 layout carries empty padding slots, which
    # the toy target covers and which therefore count as truncated at the first census)
    expect_equal(f$model$tree_coverage$frac_full[f$model$tree_coverage$census == 3], 1)
  }

  # single shooting keeps one generation for the whole horizon, so the tracked set shrinks
  s = tl_model(d)
  tl_fit(s, d, loss_resolution = "tree", tree_data = list("2" = tl_target(dims), "3" = tl_target(dims)))
  expect_length(s$model$tree_pred$growth, prod(dims))
  expect_lt(sum(s$model$tree_pred$n > 0.5), sum(s$init$trees_r > 0.5))
  # ... and the erosion is reported per census rather than assumed
  expect_identical(s$model$tree_coverage$census, c(2L, 3L))
  expect_gt(s$model$tree_coverage$frac_full[1], s$model$tree_coverage$frac_full[2])
})


test_that("tree_data is checked against the layout the interval starts from", {
  skip_if_no_torch()
  d = tl_data()
  dims = dim(tl_model(d)$init$dbh_r)
  bad = function(re, ...) expect_error(tl_fit(tl_model(d), d, loss_resolution = "tree", ...), re)
  bad("requires `tree_data`")
  bad("year_sequence", tree_data = list("9" = tl_target(dims)))
  bad("starts from", tree_data = list("2" = tl_target(dims + 1L)))
  bad("`mort`", tree_data = list("2" = list(growth = array(0, dims))))
  # the weight vector is six long in both modes - one weight per response, not per resolution
  bad("length 6", weights = rep(1, 5), tree_data = list("2" = tl_target(dims)))
  # the anchor, not init_cohort, is the layout of the interval that follows it
  bad("starts from", reanchor = TRUE, anchor_cohorts = list("2" = tl_anchor(d$Nsp)),
      tree_data = list("3" = tl_target(dims)))
  # one accumulator per cohort holds one window, so a window may not reach past the census
  # before it - refused rather than silently double counting a year
  d$obs[year == 2, period_length := 2L]
  bad("overlap the one before it", tree_data = list("2" = tl_target(dims), "3" = tl_target(dims)))
})


# ---- windows longer than one year -------------------------------------------------------
# The FIA extdata carries three env years, so `period > 1` needs its own panel: 3 sites,
# 2 patches, one tree per cohort, 10 years. Mortality is set high enough that cohorts die
# INSIDE a window, which is the only condition under which the whole-window rule is
# testable at all. `lr = 0` freezes the trajectory, the six aggregate responses are NA, so
# only the tree terms are scored, and `weights = rep(1, 6)` keeps the loss hand-computable.

tll = local({
  Nsp = 2L
  Nsites = 3L
  Tmax = 10L
  SEED = 11
  init_trees = data.table(siteID = rep(1:Nsites, each = 8), patchID = rep(rep(1:2, each = 4), Nsites),
                          species = rep(1:2, 12), dbh = rep(c(10, 20, 30, 40), 6) + rep(0:2, each = 8)*5,
                          treeName = paste0("t", 1:24), trees = 1)
  env = data.table(expand.grid(siteID = 1:Nsites, year = 1:Tmax))[, env1 := (siteID - 2)*0.5]
  list(
    Nsp = Nsp, dims = c(Nsites, 2L, 4L),
    obs = function(years, period) data.table::CJ(siteID = 1:Nsites, year = years, species = 1:Nsp)[
      , `:=`(species_name = paste0("sp", species), period_length = period, dbh = NA_real_,
             ba = NA_real_, trees = NA_real_, growth = NA_real_, mort = NA_real_, reg = NA_real_)][],
    target = function() list(growth = array(0.02, dim = c(Nsites, 2L, 4L)),
                            mort = array(rep(c(0, 1), length.out = 24), dim = c(Nsites, 2L, 4L))),
    fit = function(obs, td, period_scale = FALSE, weights = rep(1, 6), ...) {
      FINN.seed(SEED)
      ic = makeInitCohorts(copy(init_trees), Nspecies = Nsp)
      m = finn(N_species = Nsp, recruits_dbh = 5,
        competition_process  = createProcess(~0, func = FINN::competition),
        growth_process       = createProcess(~1 + env1, initEnv = matrix(c(0.5, 0.5, 0.4, -0.4), Nsp, 2),
                                             initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp, 2), func = FINN::growth),
        mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(-0.5, -0.5, 0.3, -0.3), Nsp, 2),
                                             initSpecies = matrix(c(0.2, -0.2, 0.1, 0.1, 0, 0), Nsp, 3), func = FINN::mortality),
        regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(1.5, 1.0, 0.5, -0.5), Nsp, 2),
                                             initSpecies = c(0.1, 0.1), func = FINN::regeneration, sample_regeneration = FALSE))
      m$growth_period_scale = period_scale
      suppressMessages(fit(m, env = env, data = obs, init_cohort = ic, patches = 2L, patch_size = 0.1,
          env_autoscale = FALSE, epochs = 1L, lr = 0, weights = weights, batchsize = Nsites,
          shuffle = FALSE, plot_progress = FALSE, device = "cpu",
          loss = c(dbh = "mse", ba = "mse", trees = "mse", growth = "mse",
                   mortality = "binomial", regeneration = "mse"),
          loss_resolution = "tree", tree_data = td, ...))
      m
    })
})


test_that("growth is scored only where the cohort lived through the WHOLE window", {
  skip_if_no_torch()
  # A cohort killed in year 2 of a 4-year window carries a 2-year increment, which is not a
  # prediction of the interval's growth; the mortality term keeps it (a 2-year hazard against
  # a 4-year fate is an attenuated prediction, not a different quantity - see tree_coverage).
  m = tll$fit(tll$obs(4L, 4L), list("4" = tll$target()))
  tp = m$tree_pred
  y = as.numeric(aperm(tll$target()$growth, c(3, 2, 1)))
  whole = tp$n >= 3.5
  scored = tp$n > 0.5
  expect_gt(sum(scored & !whole), 4)                 # truncated slots exist to be excluded
  expect_equal(as.numeric(m$history[[1]])[4], mean((tp$growth[whole] - y[whole])^2), tolerance = 1e-6)
  # scoring truncated increments as interval growth is the mutation the one-year suite cannot
  # see; here it is a different number by a factor of two
  expect_gt(abs(mean((tp$growth[scored] - y[scored])^2) - mean((tp$growth[whole] - y[whole])^2)), 1e-3)
  expect_equal(m$tree_coverage$frac_full, mean(whole), tolerance = 1e-12)
  expect_equal(m$tree_coverage$mean_n_over_L, mean(tp$n[scored])/4, tolerance = 1e-6)
})


test_that("the accumulators cover the window, not everything since the last census", {
  skip_if_no_torch()
  # period 3 with censuses 5 years apart: zeroing at the census (or only at year 1) would
  # accumulate 5 years into a 3-year window, which is exactly what `n` reports.
  one = tll$fit(tll$obs(5L, 3L), list("5" = tll$target()))
  expect_equal(max(one$tree_pred$n), 3)
  expect_gt(sum(one$tree_pred$n == 3), 4)

  two = tll$fit(tll$obs(c(5L, 10L), 3L), list("5" = tll$target(), "10" = tll$target()))
  expect_equal(max(two$tree_pred$n), 3)              # the second window, read after the fit
  expect_identical(two$tree_coverage$census, c(5L, 10L))
  expect_true(all(two$tree_coverage$n_target == 24))
  # single shooting conditions the tracked set on ever more simulated survival
  expect_gt(two$tree_coverage$frac_full[1], two$tree_coverage$frac_full[2])
})


test_that("update_step and growth_period_scale reach the per-tree terms", {
  skip_if_no_torch()
  a = tll$fit(tll$obs(4L, 4L), list("4" = tll$target()))
  b = tll$fit(tll$obs(4L, 4L), list("4" = tll$target()), update_step = 4L)
  # where the gradient is truncated cannot change the estimator: at lr = 0 the same
  # trajectory must give the same loss and the same per-slot prediction
  expect_identical(a$tree_pred, b$tree_pred)
  expect_equal(as.numeric(a$history[[1]]), as.numeric(b$history[[1]]), tolerance = 1e-12)

  p = tll$fit(tll$obs(4L, 4L), list("4" = tll$target()), period_scale = TRUE)
  expect_identical(p$tree_pred$n, a$tree_pred$n)
  expect_equal(p$tree_pred$mort, a$tree_pred$mort, tolerance = 1e-12)   # mortality is unscaled
  whole = a$tree_pred$n >= 3.5
  # compounded over the window, and never above the same rate compounded (AM >= GM)
  expect_gt(mean(p$tree_pred$growth[whole]), 3*mean(a$tree_pred$growth[whole]))
  expect_true(all(1 + p$tree_pred$growth[whole] <= (1 + a$tree_pred$growth[whole])^4 + 1e-8))
})


test_that("the per-tree terms carry a gradient to the growth and mortality kernels", {
  skip_if_no_torch()
  # only the two tree terms are weighted, so any gradient at all comes through them
  gr = function(m, p) { g = m$parameters[[p]]$grad
                        if (is.null(g) || prod(g$shape) == 0) 0 else as.numeric(g$abs()$sum()) }
  for (us in c(1L, 4L)) {
    m = tll$fit(tll$obs(4L, 4L), list("4" = tll$target()), update_step = us,
                weights = c(0, 0, 0, 1, 1, 0))
    expect_gt(gr(m, "nn_growth.0.weight"), 0)
    expect_gt(gr(m, "nn_mortality.0.weight"), 0)
    expect_equal(gr(m, "nn_regeneration.0.weight"), 0)   # regeneration stays a count
  }
})
