# rb_gradient keeps the forward pass identical and only changes the gradient path.
# Without this check the straight-through substitution could silently alter the
# simulated stand, which is the one thing it must not do (FINNverse/FINN#53).

test_that("rb_gradient leaves the forward pass and the default behaviour untouched", {
  skip_if_not(torch::torch_is_installed())
  library(data.table)
  sites <- 1:3
  Nsp   <- max(fia_obs_dt$species)
  init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
  mk <- function(draw = "gamma") {
    FINN.seed(1)
    finn(N_species = Nsp, recruits_dbh = 12.9,
         growth_process       = createProcess(~ temp + prec, FINN::growth),
         mortality_process    = createProcess(~ temp + prec, FINN::mortality, mortality_draw = draw),
         regeneration_process = createProcess(~ temp + prec, FINN::regeneration))
  }
  run <- function(rb, draw = "gamma") {
    m <- mk(draw)
    FINN.seed(2)
    fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
        init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
        rb_gradient = rb, plot_progress = FALSE)
    as.data.table(do.call(rbind, m$history))
  }
  # With the sampler held fixed, rb_gradient changes only the gradient: the simulated
  # stand, and therefore epoch 1's losses, are identical.
  expect_equal(run(TRUE)[1, ], run(FALSE)[1, ], tolerance = 1e-6)

  # `mortality_draw = "binomial"` is the exception, and deliberately so: it resolves to
  # torch's exact (non-differentiable) binomial under rb_gradient and to the Gamma draw
  # without it, so there the two do NOT simulate the same stand.
  expect_false(isTRUE(all.equal(run(TRUE, "binomial")[1, ], run(FALSE, "binomial")[1, ])))
})

test_that("trees_offset reaches both count losses and NULL keeps the historical constants", {
  skip_if_not(torch::torch_is_installed())
  library(data.table)
  sites <- 1:3
  Nsp   <- max(fia_obs_dt$species)
  init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
  loss1 <- function(off, family = "poisson") {
    FINN.seed(1)
    m <- finn(N_species = Nsp, recruits_dbh = 12.9,
              growth_process       = createProcess(~ temp + prec, FINN::growth),
              mortality_process    = createProcess(~ temp + prec, FINN::mortality),
              regeneration_process = createProcess(~ temp + prec, FINN::regeneration))
    FINN.seed(2)
    fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
        init_cohort = init, epochs = 1L, patches = 4, patch_size = 0.06, lr = 0,
        loss_family = c(dbh = "mse", ba = "mse", trees = family, growth = "mse",
                        mortality = "binomial", regeneration = "nbinom"),
        trees_offset = off, plot_progress = FALSE)
    m$history[[1]][3]
  }
  expect_false(isTRUE(all.equal(loss1(1e-5), loss1(1.0))))   # the offset is actually used
  expect_false(isTRUE(all.equal(loss1(0.001, "nbinom"), loss1(1.0, "nbinom"))))
  # NULL leaves each family on the constant it used before the argument existed
  expect_equal(loss1(NULL), loss1(1e-5))
  expect_equal(loss1(NULL, "nbinom"), loss1(0.001, "nbinom"))
})

test_that("a model that never went through fit() uses the same death draw as a fitted one (#70)", {
  skip_if_not(torch::torch_is_installed())
  library(data.table)
  sites <- 1:3
  Nsp   <- max(fia_obs_dt$species)
  init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
  # an unfitted model has learned no climate scaling, so give it z-scored climate
  env   <- copy(fia_env_dt[siteID %in% sites])[, `:=`(temp = as.numeric(scale(temp)), prec = as.numeric(scale(prec)))]
  mk <- function() {
    FINN.seed(1)
    finn(N_species = Nsp, recruits_dbh = 12.9,
         growth_process       = createProcess(~ temp + prec, FINN::growth),
         mortality_process    = createProcess(~ temp + prec, FINN::mortality),   # default mortality_draw = "binomial"
         regeneration_process = createProcess(~ temp + prec, FINN::regeneration))
  }
  sim <- function(m) {
    FINN.seed(2)
    as.data.table(simulateForest(m, env = env, init_cohort = init, patches = 4, patch_size = 0.06)$wide$site)
  }
  m <- mk()
  expect_true(isTRUE(m$rb_gradient))   # the default of fit(rb_gradient =), set at build time

  # the unfitted model simulates exactly as with the flag set by hand (the exact binomial) ...
  m_true <- mk(); m_true$rb_gradient <- TRUE
  expect_equal(sim(m), sim(m_true))
  # ... and not as the Gamma draw that a missing flag used to select
  m_gamma <- mk(); m_gamma$rb_gradient <- FALSE
  expect_false(isTRUE(all.equal(sim(m)$trees, sim(m_gamma)$trees)))
})
