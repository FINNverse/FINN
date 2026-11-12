## The Beverton-Holt cap is FINN's to create and fit()'s to initialise: the user
## declares nothing, `cap` sets the resolution, and `prior = "max_obs_recruits"`
## can only be resolved once the data are there (FINNverse/FINN#58).
test_that("finn() creates the cap itself and cap = sets its resolution", {
  skip_if_no_torch()
  m <- finn(N_species = 3)
  expect_false(is.null(m$reg_logK))
  expect_equal(m$reg_logK$numel(), 1L)          # cap = "all" is the default
  expect_equal(m$reg_cap, "all")

  ms <- finn(N_species = 3, regeneration_process =
               createProcess(~1, func = FINN::regeneration_capped, cap = "species"))
  expect_equal(ms$reg_logK$numel(), 3L)

  mu <- finn(N_species = 3, regeneration_process =
               createProcess(~1, func = FINN::regeneration))
  expect_null(mu$reg_logK)                      # the uncapped process has no cap
})

test_that("the prior is validated at construction", {
  skip_if_no_torch()
  expect_error(createProcess(func = FINN::regeneration_capped, prior = -1), "max_obs_recruits")
  expect_error(createProcess(func = FINN::regeneration_capped, prior_lambda = -1), "prior_lambda")
  expect_equal(createProcess(func = FINN::regeneration_capped, prior = 500)$prior, 500)
})

test_that("fit() resolves max_obs_recruits and starts the cap there", {
  skip_if_no_torch()
  library(data.table)
  sites <- 1:3
  Nsp   <- max(fia_obs_dt$species)
  init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
  d     <- fia_obs_dt[siteID %in% sites]
  run <- function(...) {
    FINN.seed(1)
    m <- finn(N_species = Nsp, recruits_dbh = 12.9, ...)
    FINN.seed(2)
    fit(m, data = d, env = fia_env_dt[siteID %in% sites], init_cohort = init,
        epochs = 0L, patches = 4, patch_size = 0.06, plot_progress = FALSE)
    m
  }
  m <- run()
  expect_equal(m$reg_prior_value, max(d$reg, na.rm = TRUE))
  # epochs = 0 is set-up only, so the cap still sits exactly on the prior
  expect_equal(as.numeric(m$reg_logK)[1], log(m$reg_prior_value), tolerance = 1e-6)

  mf <- run(regeneration_process = createProcess(~1, func = FINN::regeneration_capped, prior = 500))
  expect_equal(mf$reg_prior_value, 500)
})
