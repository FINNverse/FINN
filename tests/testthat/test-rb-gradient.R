# rb_gradient keeps the forward pass identical and only changes the gradient path.
# Without this check the straight-through substitution could silently alter the
# simulated stand, which is the one thing it must not do (FINNverse/FINN#53).

test_that("rb_gradient leaves the forward pass and the default behaviour untouched", {
  skip_if_not(torch::torch_is_installed())
  library(data.table)
  sites <- 1:3
  Nsp   <- max(fia_obs_dt$species)
  init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
  mk <- function() {
    FINN.seed(1)
    finn(N_species = Nsp, recruits_dbh = 12.9,
         growth_process       = createProcess(~ temp + prec, FINN::growth),
         mortality_process    = createProcess(~ temp + prec, FINN::mortality),
         regeneration_process = createProcess(~ temp + prec, FINN::regeneration))
  }
  run <- function(rb) {
    m <- mk()
    FINN.seed(2)
    fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
        init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
        rb_gradient = rb, plot_progress = FALSE)
    as.data.table(do.call(rbind, m$history))
  }
  expect_equal(run(TRUE)[1, ], run(FALSE)[1, ], tolerance = 1e-6)   # same forward, epoch 1
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
