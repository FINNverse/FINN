## shooting = "multiple": the stand is reset to the observed tree list at every
## census, the recurrent state is not truncated inside an interval, and the
## anchored trees are tracked through the simulated interval.
##
## Checked here: the default reproduces the single-shooting objective exactly;
## the state after an anchor IS the census tree list (species, stem numbers,
## cohort ids 1..n); anchored trees the simulation kills stay in the state until
## their interval ends; the objective differs from single shooting because the
## gradient window differs; input validation.
skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)

Nsp <- 2; Nsites <- 3; Tmax <- 6; SEED <- 11
LOSS <- c(dbh = "mse", ba = "mse", trees = "mse", growth = "mse",
          mortality = "binomial", regeneration = "mse")

mk_model <- function(mort_intercept = -2.5) finn(N_species = Nsp, recruits_dbh = 5,
  competition_process  = createProcess(~0, func = FINN::competition),
  growth_process       = createProcess(~1 + env1, initEnv = matrix(c(0.5, 0.5, 0.4, -0.4), Nsp, 2),
                                       initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp, 2), func = FINN::growth),
  mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(mort_intercept, mort_intercept, 0.3, -0.3), Nsp, 2),
                                       initSpecies = matrix(c(0.2, -0.2, 0.1, 0.1, 0, 0), Nsp, 3), func = FINN::mortality),
  regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(1.5, 1.0, 0.5, -0.5), Nsp, 2),
                                       initSpecies = c(0.1, 0.1), func = FINN::regeneration,
                                       sample_regeneration = FALSE))

init_trees <- data.table(siteID = rep(1:Nsites, each = 8), patchID = rep(rep(1:2, each = 4), Nsites),
                         species = rep(1:2, 12), dbh = rep(c(10, 20, 30, 40), 6) + rep(0:2, each = 8) * 5,
                         treeName = paste0("t", 1:24), trees = 3)
ic  <- makeInitCohorts(init_trees, Nspecies = Nsp)
env <- data.table(expand.grid(siteID = 1:Nsites, year = 1:Tmax)); env[, env1 := (siteID - 2) * 0.5]
## censuses at years 3 and 6; site 1 has a 3-year first interval, the others 2
obs <- data.table(expand.grid(siteID = 1:Nsites, year = c(3L, 6L), species = 1:Nsp))
obs[, `:=`(dbh = 20, ba = 10, trees = 300, growth = 0.02, mort = 0.02, reg = 5,
           period_length = fifelse(year == 3L & siteID == 1L, 3L, 2L),
           species_name = paste0("sp", species))]
setorder(obs, siteID, year, species)
## the tree list that opens each interval: 4 trees per patch, distinct per census
## so an anchor is recognisable (species pattern flips between the censuses)
tree_dt <- data.table(expand.grid(siteID = 1:Nsites, year = c(3L, 6L), patchID = 1:2, slot = 1:4))
tree_dt[, `:=`(species = fifelse(year == 3L, (slot %% 2L) + 1L, 2L - (slot %% 2L)),   # 2,1,2,1 at census 3; 1,2,1,2 at census 6
               dbh = 8 * slot + fifelse(year == 6L, 50, 0),
               trees = fifelse(year == 3L, 3, 7), growth = 0.02, died = as.numeric(slot == 3L))]
setorder(tree_dt, siteID, year, patchID, slot)
tree_dt[, slot := NULL]

## the recorded state of one cohort id (slots are re-sorted every year, the id
## travels with the slot)
by_id <- function(out, step, site, patch, id, what) {
  ids <- out$Predictions$Cohort$cohortID[[step]][site, patch, ]
  w <- which(ids == id)
  if (!length(w)) return(NA)
  out$Predictions$Cohort$cohortStates[[step]][[what]][site, patch, w[1]]
}

OPT <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
run_fit <- function(..., data = obs, model = mk_model()) {
  FINN.seed(SEED); mm <- model
  fit(mm, env = env, data = data, init_cohort = ic, patches = 2L, patch_size = 0.1,
      env_autoscale = FALSE, epochs = 1L, lr = 0, weights = rep(1, 6), batchsize = Nsites,
      shuffle = FALSE, optimizer = OPT, plot_progress = FALSE, device = "cpu",
      loss_family = LOSS, ...)
  mm
}

## forward() in fit mode with the cohort record switched on, so the state after
## every step can be inspected. Mirrors what fit() hands to forward().
fwd_record <- function(mm, shooting) {
  ys <- which(levels(as.factor(env$year)) %in% levels(as.factor(obs$year)))
  tob <- lapply(mm$tree_obs[c("dbh", "trees", "species", "growth", "died")], function(t) t)
  ## the response tensor, rebuilt as fit() does (period_length in slice 7)
  resp <- lapply(c("dbh", "ba", "trees", "growth", "mort", "reg", "period_length"), function(v)
    abind::abind(lapply(1:Nsp, function(i) FINN:::extract_env(stats::as.formula(paste0("~0+", v)), obs[species == i])), along = 3L))
  resp <- c(resp, list(array(1, dim = dim(resp[[1]])), array(1, dim = dim(resp[[1]]))))
  y <- torch::torch_cat(lapply(resp, function(a) torch::torch_tensor(a, dtype = torch::torch_float32())$unsqueeze(4)), 4)
  ## the design matrices the processes expect (~1 + env1: intercept + env1)
  e <- mm$.__enclos_env__$private$extract_env_method(env)
  FINN.seed(SEED)
  mm$forward(dbh = ic$dbh, trees = ic$trees, species = ic$species,
             env = list(torch::torch_tensor(e$mortality_env), torch::torch_tensor(e$growth_env),
                        torch::torch_tensor(e$regeneration_env)),
             y = y, tree_obs = tob, year_sequence = ys, shooting = shooting,
             return_cohorts = TRUE, verbose = FALSE)
}

test_that("the default reproduces the single-shooting objective", {
  a <- run_fit(loss_aggregation = list(growth = "tree", mortality = "tree"), tree_data = tree_dt)
  b <- run_fit(loss_aggregation = list(growth = "tree", mortality = "tree"), tree_data = tree_dt,
               shooting = "single")
  expect_equal(as.numeric(a$history[[1]]), as.numeric(b$history[[1]]), tolerance = 1e-10)
  expect_identical(b$shooting, "single")
})

test_that("multiple shooting needs the tree list and forces the per-site path", {
  expect_error(run_fit(shooting = "multiple"), "needs `tree_data`")
  m <- run_fit(shooting = "multiple", tree_data = tree_dt)
  expect_identical(m$shooting, "multiple")
  expect_true(m$per_site_period)
  h <- as.numeric(m$history[[1]])
  expect_true(all(is.finite(h)))
})

test_that("the state after an anchor is the census tree list, with cohort ids 1..n", {
  m <- run_fit(shooting = "multiple", tree_data = tree_dt)
  out <- fwd_record(m, "multiple")
  ## Intervals: site 1 has 1..3 and 5..6 (opens at steps 1 and 5, year 4 belongs
  ## to no interval); sites 2 and 3 have 2..3 and 5..6 (open at steps 2 and 5).
  ## The state recorded after step 5
  ## for site 2 must therefore carry the census-6 list (species 1,2,1,2; dbh
  ## 58..82 before the year's growth; 7 stems per slot) under ids 1..4.
  sp5 <- vapply(1:4, function(j) by_id(out, 5, 2, 1, j, "species"), numeric(1))
  expect_equal(as.integer(sp5), c(1L, 2L, 1L, 2L))
  d5 <- vapply(1:4, function(j) by_id(out, 5, 2, 1, j, "dbh"), numeric(1))
  expect_true(all(d5 >= 8 * (1:4) + 50 - 1e-4))          # the census diameter, grown one year
  expect_true(all(d5 < 8 * (1:4) + 50 + 10))
  n5 <- vapply(1:4, function(j) by_id(out, 5, 2, 1, j, "trees"), numeric(1))
  expect_true(all(n5 <= 7 & n5 >= 5))
  ## Step 4 opens nothing for site 2: its state is the simulated continuation of
  ## the census-3 list (species 2,1,2,1), not the census-6 one.
  sp4 <- vapply(1:4, function(j) by_id(out, 4, 2, 1, j, "species"), numeric(1))
  expect_equal(as.integer(sp4), c(2L, 1L, 2L, 1L))
  d4 <- vapply(1:4, function(j) by_id(out, 4, 2, 1, j, "dbh"), numeric(1))
  expect_true(all(d4 > 8 * (1:4) & d4 < 58))
  ## site 1 runs on from its census-3 list through the gap year 4 and is reset
  ## at step 5 like the others
  sp4_1 <- vapply(1:4, function(j) by_id(out, 4, 1, 1, j, "species"), numeric(1))
  expect_equal(as.integer(sp4_1), c(2L, 1L, 2L, 1L))
  sp5_1 <- vapply(1:4, function(j) by_id(out, 5, 1, 1, j, "species"), numeric(1))
  expect_equal(as.integer(sp5_1), c(1L, 2L, 1L, 2L))
})

test_that("anchored trees the simulation kills stay in the state until the interval ends", {
  ## mortality intercept 5 -> annual death probability ~ 0.99: every anchored
  ## tree dies in the first year of its interval
  m <- run_fit(shooting = "multiple", tree_data = tree_dt, model = mk_model(mort_intercept = 5))
  out <- fwd_record(m, "multiple")
  ## step 3 lies inside every site's first interval, two or three years after
  ## its anchor: every anchored tree is dead, yet still recorded under its id,
  ## with a diameter that kept growing
  for (s in 1:Nsites) for (j in 1:4) {
    expect_lt(by_id(out, 3, s, 1, j, "trees"), 1e-4)   # float residue of the straight-through rounding
    expect_true(by_id(out, 3, s, 1, j, "dbh") > 8 * j)
  }
  ## the same model under single shooting prunes the dead cohorts
  m1 <- run_fit(shooting = "single", tree_data = tree_dt, model = mk_model(mort_intercept = 5),
                loss_aggregation = list(growth = "tree", mortality = "tree"))
  out1 <- fwd_record(m1, "single")
  expect_true(all(is.na(vapply(1:4, function(j) by_id(out1, 3, 1, 1, j, "trees"), numeric(1)))))
})

test_that("tracked per-tree terms score and differ from the teacher-forced ones", {
  agg <- list(growth = "tree", mortality = "tree")
  s <- run_fit(shooting = "single",   tree_data = tree_dt, loss_aggregation = agg)
  m <- run_fit(shooting = "multiple", tree_data = tree_dt, loss_aggregation = agg)
  hs <- as.numeric(s$history[[1]]); hm <- as.numeric(m$history[[1]])
  expect_true(all(is.finite(hm)) && all(hm[4:5] > 0))
  ## tracked trees grow and see a changing stand, the held stand does not
  expect_false(isTRUE(all.equal(hs[4:5], hm[4:5], tolerance = 1e-6)))
})

test_that("parameters receive gradients under multiple shooting", {
  agg <- list(growth = "tree", mortality = "tree")
  FINN.seed(SEED); mm <- mk_model()
  fit(mm, env = env, data = obs, init_cohort = ic, patches = 2L, patch_size = 0.1,
      env_autoscale = FALSE, epochs = 2L, lr = 0.05, weights = rep(1, 6), batchsize = Nsites,
      shuffle = FALSE, optimizer = OPT, plot_progress = FALSE, device = "cpu",
      loss_family = LOSS, loss_aggregation = agg, tree_data = tree_dt, shooting = "multiple",
      record_gradients = TRUE, checkpoints = 1L)
  g <- mm$gradients[[1]]
  expect_true(length(g) > 0)
  expect_true(any(vapply(g, function(x) !is.null(x) && any(is.finite(as.numeric(x)) & as.numeric(x) != 0), logical(1))))
  expect_true(all(is.finite(as.numeric(mm$history[[2]]))))
})
