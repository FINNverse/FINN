## loss_family / loss_aggregation: the declarative loss API.
##
## Checked here: the old `loss` argument still works and gives the identical
## objective (so nothing that exists today changes), "none" switches a response
## off cleanly, the size-class operator converges on hard diameter binning as
## tau -> 0 while staying differentiable at a break, the tree-level operator
## runs, is scored on the observed trees and averages over the interval window,
## weights = "auto" gives the new aggregations their own baseline, and a single
## non-finite response no longer discards the whole timestep.
skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)

Nsp <- 2; Nsites <- 3; Tmax <- 6; SEED <- 11
LOSS <- c(dbh = "mse", ba = "mse", trees = "mse", growth = "mse",
          mortality = "mse", regeneration = "mse")

mk_model <- function() finn(N_species = Nsp, recruits_dbh = 5,
  competition_process  = createProcess(~0, func = FINN::competition),
  growth_process       = createProcess(~1 + env1, initEnv = matrix(c(0.5, 0.5, 0.4, -0.4), Nsp, 2),
                                       initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp, 2), func = FINN::growth),
  mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(-2.5, -2.5, 0.3, -0.3), Nsp, 2),
                                       initSpecies = matrix(c(0.2, -0.2, 0.1, 0.1, 0, 0), Nsp, 3), func = FINN::mortality),
  regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(1.5, 1.0, 0.5, -0.5), Nsp, 2),
                                       initSpecies = c(0.1, 0.1), func = FINN::regeneration,
                                       sample_regeneration = FALSE))

init_trees <- data.table(siteID = rep(1:Nsites, each = 8), patchID = rep(rep(1:2, each = 4), Nsites),
                         species = rep(1:2, 12), dbh = rep(c(10, 20, 30, 40), 6) + rep(0:2, each = 8) * 5,
                         treeName = paste0("t", 1:24), trees = 3)
ic  <- makeInitCohorts(init_trees, Nspecies = Nsp)
env <- data.table(expand.grid(siteID = 1:Nsites, year = 1:Tmax)); env[, env1 := (siteID - 2) * 0.5]
obs <- data.table(expand.grid(siteID = 1:Nsites, year = c(3L, 6L), species = 1:Nsp))
## period_length differs between sites, which puts fit() on the per-site
## aggregation-window path -- the branch the tree-level overrides live in. With a
## constant period they would never be exercised at all.
obs[, `:=`(dbh = 20, ba = 10, trees = 300, growth = 0.02, mort = 0.02, reg = 5,
           period_length = fifelse(siteID == 1L, 3L, 2L),
           species_name = paste0("sp", species))]
setorder(obs, siteID, year, species)

OPT <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
run_fit <- function(..., data = obs) {
  FINN.seed(SEED); mm <- mk_model()
  fit(mm, env = env, data = data, init_cohort = ic, patches = 2L, patch_size = 0.1,
      env_autoscale = FALSE, epochs = 1L, lr = 0, weights = rep(1, 6), batchsize = Nsites,
      shuffle = FALSE, optimizer = OPT, plot_progress = FALSE, device = "cpu", ...)
  mm
}

test_that("`loss` is accepted with a warning and gives the identical objective", {
  old <- NULL
  expect_warning(old <- run_fit(loss = LOSS), "deprecated")
  new <- run_fit(loss_family = LOSS)
  expect_equal(as.numeric(old$history[[1]]), as.numeric(new$history[[1]]), tolerance = 1e-10)
})

test_that("loss_aggregation validates its input", {
  expect_error(run_fit(loss_family = LOSS, loss_aggregation = c(nonsense = "mean")),
               "not a response")
  expect_error(run_fit(loss_family = LOSS, loss_aggregation = c(trees = "size_classes")),
               "needs `breaks`")
  expect_error(run_fit(loss_family = LOSS, loss_aggregation = c(dbh = "tree")),
               "'growth' and 'mortality' only")
  expect_error(run_fit(loss_family = LOSS, loss_aggregation = c(trees = "median")),
               "Unsupported loss_aggregation")
  expect_error(run_fit(loss_family = LOSS, loss_aggregation = c(growth = "tree")),
               "needs `tree_data`")
})

test_that("`none` switches a response off without changing the others", {
  full <- run_fit(loss_family = LOSS)
  off  <- run_fit(loss_family = replace(LOSS, "dbh", "none"))
  h_full <- as.numeric(full$history[[1]]); h_off <- as.numeric(off$history[[1]])
  expect_equal(h_off[1], 0)
  expect_equal(h_off[2:6], h_full[2:6], tolerance = 1e-10)
})

test_that("size-class counts converge on hard binning as tau -> 0", {
  m <- mk_model()
  cc <- m$.__enclos_env__$private$class_counts
  dbh     <- torch::torch_tensor(array(c(5, 15, 25, 35), dim = c(1, 1, 4)))
  trees   <- torch::torch_tensor(array(c(2, 3, 4, 5),    dim = c(1, 1, 4)))
  species <- torch::torch_tensor(array(c(1L, 1L, 2L, 2L), dim = c(1, 1, 4)), dtype = torch::torch_int64())
  breaks  <- c(10, 20, 30)
  got  <- as.array(cc(dbh, trees, species, breaks, tau = 1e-3))[1, , ]
  # hard binning: sp1 has 2 stems in (-Inf,10] and 3 in (10,20]; sp2 4 in (20,30] and 5 in (30,Inf)
  want <- matrix(c(2, 3, 0, 0,
                   0, 0, 4, 5), nrow = 2, byrow = TRUE)
  expect_equal(dim(got), c(2L, 4L))
  expect_equal(got, want, tolerance = 1e-3)
})

test_that("size-class counts are soft (and differentiable) at a break", {
  m <- mk_model()
  cc <- m$.__enclos_env__$private$class_counts
  dbh <- torch::torch_tensor(array(20, dim = c(1, 1, 1)), requires_grad = TRUE)
  out <- cc(dbh, torch::torch_tensor(array(1, dim = c(1, 1, 1))),
            torch::torch_tensor(array(1L, dim = c(1, 1, 1)), dtype = torch::torch_int64()),
            breaks = 20, tau = 1)
  # a tree exactly on the break splits evenly between the two classes
  expect_equal(as.array(out)[1, 1, ], c(0.5, 0.5), tolerance = 1e-5)
  out[1, 1, 1]$backward()
  expect_true(is.finite(as.numeric(dbh$grad)) && as.numeric(dbh$grad) != 0)
})

test_that("tree-level growth and mortality run and are scored on the observed trees", {
  tree_dt <- data.table(expand.grid(siteID = 1:Nsites, year = c(3L, 6L), patchID = 1:2, slot = 1:4))
  tree_dt[, `:=`(species = rep_len(1:Nsp, .N), dbh = rep_len(c(8, 16, 24, 32), .N),
                 trees = 3, growth = 0.02, died = rep_len(c(0, 0, 1, 0), .N))]
  tree_dt[, slot := NULL]
  agg <- list(trees = list(type = "size_classes", breaks = c(10, 20, 30)),
              growth = "tree", mortality = "tree")
  obs_cls <- copy(obs)
  for (k in 1:4) obs_cls[[paste0("n_class_", k)]] <- 50
  m <- run_fit(loss_family = replace(LOSS, c("dbh", "ba"), c("none", "none")),
               loss_aggregation = agg, tree_data = tree_dt, data = obs_cls)
  h <- as.numeric(m$history[[1]])
  expect_equal(h[1:2], c(0, 0))                       # dbh and ba switched off
  expect_true(all(is.finite(h)))
  expect_true(all(h[3:5] > 0))                        # counts, growth, mortality all score
  expect_equal(dim(as.array(m$class_obs)), c(Nsites, 2L, Nsp, 4L))
  expect_equal(dim(as.array(m$tree_obs$dbh)), c(Nsites, 2L, 2L, 4L))
})

test_that("weights = 'auto' gives the new aggregations their own baseline", {
  tree_dt <- data.table(expand.grid(siteID = 1:Nsites, year = c(3L, 6L), patchID = 1:2, slot = 1:4))
  tree_dt[, `:=`(species = rep_len(1:Nsp, .N), dbh = rep_len(c(8, 16, 24, 32), .N),
                 trees = 3, growth = rep_len(c(0.01, 0.03), .N), died = rep_len(c(0, 0, 1, 0), .N))]
  tree_dt[, slot := NULL]
  obs_cls <- copy(obs)
  for (k in 1:4) obs_cls[[paste0("n_class_", k)]] <- c(80, 60, 40, 20)[k]
  FINN.seed(SEED); mm <- mk_model()
  fit(mm, env = env, data = obs_cls, init_cohort = ic, patches = 2L, patch_size = 0.1,
      env_autoscale = FALSE, epochs = 1L, lr = 0, weights = "auto", batchsize = Nsites,
      shuffle = FALSE, optimizer = OPT, plot_progress = FALSE, device = "cpu",
      loss_family = c(dbh = "none", ba = "none", trees = "poisson", growth = "mse",
                      mortality = "binomial", regeneration = "mse"),
      loss_aggregation = list(trees = list(type = "size_classes", breaks = c(10, 20, 30)),
                              growth = "tree", mortality = "tree"),
      tree_data = tree_dt)
  b <- mm$loss_baseline
  # trees/growth/mortality no longer live in Y; without their own baseline these
  # would be NA and silently fall back to a weight of 1.
  expect_true(all(is.finite(b[c("trees", "growth", "mortality")])))
  expect_true(all(b[c("trees", "growth", "mortality")] > 0))
})

test_that("the tree operator averages over the interval window", {
  tree_dt <- data.table(expand.grid(siteID = 1:Nsites, year = c(3L, 6L), patchID = 1:2, slot = 1:4))
  tree_dt[, `:=`(species = rep_len(1:Nsp, .N), dbh = rep_len(c(8, 16, 24, 32), .N),
                 trees = 3, growth = 0.02, died = rep_len(c(0, 0, 1, 0), .N))]
  tree_dt[, slot := NULL]
  obs_cls <- copy(obs); for (k in 1:4) obs_cls[[paste0("n_class_", k)]] <- 50
  m <- run_fit(loss_family = replace(LOSS, c("dbh", "ba"), c("none", "none")),
               loss_aggregation = list(trees = list(type = "size_classes", breaks = c(10, 20, 30)),
                                       growth = "tree", mortality = "tree"),
               tree_data = tree_dt, data = obs_cls)
  tp   <- m$.__enclos_env__$private$tree_predictions
  envt <- list(growth = torch::torch_randn(c(Nsites, 6, 2)), mort = torch::torch_randn(c(Nsites, 6, 2)))
  # make year 5 and 6 identical, so a 2-year window must return the same growth
  # rate as a 1-year window ending at 6 -- the mean of two identical values
  envt$growth[, 5, ] <- envt$growth[, 6, ]; envt$mort[, 5, ] <- envt$mort[, 6, ]
  mk_M <- function(p) {
    M <- torch::torch_zeros(c(Nsites, 6)); for (j in (6 - p + 1):6) M[, j] <- 1; M
  }
  torch::with_no_grad({
    a <- tp(m$tree_obs, 2L, envt, 6L, mk_M(1), torch::torch_ones(Nsites))
    b <- tp(m$tree_obs, 2L, envt, 6L, mk_M(2), torch::torch_ones(Nsites) * 2)
  })
  expect_equal(as.array(a$growth), as.array(b$growth), tolerance = 1e-5)
  # mortality is a survival PRODUCT over the window, so two identical years give
  # a strictly higher interval death probability than one
  expect_true(all(as.array(b$mortality) >= as.array(a$mortality) - 1e-6))
  expect_true(mean(as.array(b$mortality)) > mean(as.array(a$mortality)))
})

test_that("one non-finite response no longer discards the other five", {
  bad <- copy(obs)[, trees := -1]      # a negative count makes the Poisson term NaN
  expect_warning(m <- run_fit(loss_family = replace(LOSS, "trees", "poisson"), data = bad),
                 "Non-finite loss skipped for trees")
  h <- as.numeric(m$history[[1]])
  expect_true(m$nonfinite_counts[3] > 0)
  expect_true(all(is.finite(h[c(1, 2, 4, 5, 6)])))
  expect_true(all(h[c(1, 2, 4, 5, 6)] > 0))   # the other responses still scored
})

test_that("predictTrees returns one prediction per observed tree", {
  tree_dt <- data.table(expand.grid(siteID = 1:Nsites, year = c(3L, 6L), patchID = 1:2, slot = 1:4))
  tree_dt[, `:=`(species = rep_len(1:Nsp, .N), dbh = rep_len(c(8, 16, 24, 32), .N),
                 trees = 3, growth = 0.02, died = rep_len(c(0, 0, 1, 0), .N),
                 period_length = 3L)]
  tree_dt[, slot := NULL]
  obs_cls <- copy(obs); for (k in 1:4) obs_cls[[paste0("n_class_", k)]] <- 50
  m <- run_fit(loss_family = replace(LOSS, c("dbh", "ba"), c("none", "none")),
               loss_aggregation = list(trees = list(type = "size_classes", breaks = c(10, 20, 30)),
                                       growth = "tree", mortality = "tree"),
               tree_data = tree_dt, data = obs_cls)
  P <- predictTrees(m, tree_dt, env, patches = 2L, patch_size = 0.1)
  expect_equal(nrow(P), nrow(tree_dt))
  expect_true(all(is.finite(P$growth_pred)) && all(is.finite(P$mort_pred)))
  expect_true(all(P$mort_pred >= 0 & P$mort_pred <= 1))
  expect_true(all(c("siteID", "year", "patchID") %in% names(P)))
  # larger trees grow relatively less: the kernel's exp(-b * dbh) term
  D <- merge(P, tree_dt[, .(siteID, year, patchID, dbh)], by = c("siteID", "year", "patchID"), allow.cartesian = TRUE)
  expect_lt(cor(D$dbh, D$growth_pred), 0)
})

test_that("the quantile operator reproduces weighted diameter quantiles", {
  m  <- mk_model()
  qv <- m$.__enclos_env__$private$quantile_values
  ## one site, one patch, five cohorts of species 1 with known abundances
  d  <- c(10, 20, 30, 40, 50); w <- c(1, 1, 1, 1, 1)
  dbh     <- torch::torch_tensor(array(d, dim = c(1, 1, 5)))
  trees   <- torch::torch_tensor(array(w, dim = c(1, 1, 5)))
  species <- torch::torch_tensor(array(rep(1L, 5), dim = c(1, 1, 5)), dtype = torch::torch_int64())
  probs <- c(0.1, 0.5, 0.9)
  got <- as.array(qv(dbh, trees, species, probs, tau = 1e-3))[1, 1, ]
  ## with equal weights the exact weighted quantiles are the 1st, 3rd and 5th stem
  expect_equal(got, c(10, 30, 50), tolerance = 1e-2)
})

test_that("the quantile operator follows the abundance weights, not the values", {
  m  <- mk_model()
  qv <- m$.__enclos_env__$private$quantile_values
  ## the same diameters, but almost all stems are thin: the median must move down
  dbh     <- torch::torch_tensor(array(c(10, 20, 30, 40, 50), dim = c(1, 1, 5)))
  species <- torch::torch_tensor(array(rep(1L, 5), dim = c(1, 1, 5)), dtype = torch::torch_int64())
  even <- as.array(qv(dbh, torch::torch_tensor(array(rep(1, 5), dim = c(1, 1, 5))),
                      species, 0.5, tau = 1e-3))[1, 1, 1]
  thin <- as.array(qv(dbh, torch::torch_tensor(array(c(50, 1, 1, 1, 1), dim = c(1, 1, 5))),
                      species, 0.5, tau = 1e-3))[1, 1, 1]
  expect_equal(even, 30, tolerance = 1e-2)
  expect_equal(thin, 10, tolerance = 1e-2)
})

test_that("the quantile operator is differentiable and species-separated", {
  m  <- mk_model()
  qv <- m$.__enclos_env__$private$quantile_values
  dbh <- torch::torch_tensor(array(c(10, 20, 30, 40), dim = c(1, 1, 4)), requires_grad = TRUE)
  out <- qv(dbh, torch::torch_tensor(array(rep(1, 4), dim = c(1, 1, 4))),
            torch::torch_tensor(array(c(1L, 1L, 2L, 2L), dim = c(1, 1, 4)), dtype = torch::torch_int64()),
            probs = 0.5, tau = 1e-2)
  v <- as.array(out)[1, , 1]
  expect_true(v[1] < 25 && v[2] > 25)     # species 1 is the thin pair, species 2 the thick one
  out[1, 1, 1]$backward()
  expect_true(any(as.array(dbh$grad) != 0))
})

test_that("dbh quantiles run end to end and validate their input", {
  expect_error(run_fit(loss_family = LOSS, loss_aggregation = list(dbh = list(type = "quantiles"))),
               "needs `probs`")
  expect_error(run_fit(loss_family = LOSS,
                       loss_aggregation = list(dbh = list(type = "quantiles", probs = c(0, 0.5)))),
               "strictly between 0 and 1")
  expect_error(run_fit(loss_family = LOSS,
                       loss_aggregation = list(ba = list(type = "quantiles", probs = 0.5))),
               "'dbh' only")
  obs_q <- copy(obs)
  for (p in c(10, 50, 90)) obs_q[[sprintf("dbh_q%g", p)]] <- 15 + p / 10
  mm <- run_fit(loss_family = LOSS, data = obs_q,
                loss_aggregation = list(dbh = list(type = "quantiles", probs = c(0.1, 0.5, 0.9))))
  h <- as.numeric(mm$history[[1]])
  expect_true(all(is.finite(h)))
  expect_gt(h[1], 0)
  expect_equal(dim(as.array(mm$quant_obs)), c(Nsites, 2L, Nsp, 3L))
})
