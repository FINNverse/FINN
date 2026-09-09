## loss_family / loss_aggregation: the declarative loss API.
##
## Three things are checked: the old `loss` argument still works and gives the
## identical objective (so nothing that exists today changes), the size-class
## operator converges on hard diameter binning as tau -> 0, and the tree-level
## operator runs and scores the observed trees.
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
obs[, `:=`(dbh = 20, ba = 10, trees = 300, growth = 0.02, mort = 0.02, reg = 5,
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
