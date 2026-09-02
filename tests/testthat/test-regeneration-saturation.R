library(testthat)
library(FINN)
library(data.table)

# The FINN::regeneration_saturation process (Beverton-Holt density cap). Its
# carrying capacity K = exp(reg_logK) is a custom process parameter declared via
# createProcess(custom_parameters=): a length-1 init -> shared K, a
# length-N_species init -> per-species K. Torch-gated; tiny FIA subset so it runs
# fast. The per-species case is a regression guard for the tensor-broadcast bug.

build_sat_model <- function(logK_len, Nsp) {
  FINN.seed(1)
  finn(
    N_species            = Nsp, recruits_dbh = 12.9,
    competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    growth_process       = createProcess(~ temp + prec, FINN::growth,    optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = createProcess(~ temp + prec, FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
    regeneration_process = createProcess(
      ~ temp + prec, FINN::regeneration_saturation,
      custom_parameters = list(reg_logK = rep(log(800), logK_len)),
      optimizeSpecies = TRUE, optimizeEnv = TRUE)
  )
}

sat_data <- function(n_sites = 4L) {
  ext <- function(f) system.file("extdata", f, package = "FINN")
  obs <- data.table::fread(ext("fia_obs_dt.csv"))
  env <- data.table::fread(ext("fia_env_dt.csv"))
  it  <- data.table::fread(ext("fia_init_trees.csv"))
  sites <- utils::head(sort(unique(obs$siteID)), n_sites)
  Nsp <- max(obs$species)
  list(obs = obs[siteID %in% sites], env = env[siteID %in% sites],
       init = makeInitCohorts(it[siteID %in% sites], Nspecies = Nsp), Nsp = Nsp)
}

test_that("per-species K (reg_logK length N_species) fits without crashing", {
  skip_if_no_torch()
  d <- sat_data()
  m <- build_sat_model(d$Nsp, d$Nsp)                       # per-species K

  expect_false(is.null(m$reg_logK))                        # parameter registered
  expect_true(m$reg_logK$requires_grad)                    # and trainable
  expect_length(as.numeric(m$reg_logK), d$Nsp)             # one K per species

  before <- as.numeric(m$reg_logK)
  suppressMessages(
    fit(m, env = d$env, data = d$obs, init_cohort = d$init, device = "cpu",
        epochs = 3L, patches = 2, patch_size = 0.06, lr = 0.02, plot_progress = FALSE)
  )
  after <- as.numeric(m$reg_logK)
  expect_length(after, d$Nsp)
  expect_true(all(is.finite(after)))
  expect_false(isTRUE(all.equal(before, after)))           # K was optimised
})

test_that("shared K (reg_logK length 1) fits and stays scalar", {
  skip_if_no_torch()
  d <- sat_data()
  m <- build_sat_model(1L, d$Nsp)                          # shared K
  expect_length(as.numeric(m$reg_logK), 1L)
  suppressMessages(
    fit(m, env = d$env, data = d$obs, init_cohort = d$init, device = "cpu",
        epochs = 3L, patches = 2, patch_size = 0.06, lr = 0.02, plot_progress = FALSE)
  )
  expect_length(as.numeric(m$reg_logK), 1L)
  expect_true(is.finite(as.numeric(m$reg_logK)))
})

test_that("regeneration_saturation without reg_logK errors clearly", {
  skip_if_no_torch()
  d <- sat_data()
  m <- finn(
    N_species = d$Nsp, recruits_dbh = 12.9,
    competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    regeneration_process = createProcess(~ temp + prec, FINN::regeneration_saturation,
                                         optimizeSpecies = TRUE, optimizeEnv = TRUE))
  expect_error(
    suppressMessages(
      fit(m, env = d$env, data = d$obs, init_cohort = d$init, device = "cpu",
          epochs = 1L, patches = 2, patch_size = 0.06, plot_progress = FALSE)),
    "reg_logK")
})
