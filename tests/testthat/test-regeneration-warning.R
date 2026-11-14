library(testthat)
library(FINN)
library(data.table)

# fit() warns when the regeneration process is calibrated while no observed
# recruitment reaches the loss. Two ways that happens, and both must warn:
# the term is switched off, or the family is set but the data carry no finite
# `reg` in the scored cells. No other response can train regeneration -- recruits
# reach dbh through a near-step sigmoid and trees through a rounding that is 0
# below 0.5 recruits per patch and species (FINNverse/FINN#43).

reg_warn_setup <- function() {
  obs <- data.table::as.data.table(FINN::fia_obs_dt)
  env <- data.table::as.data.table(FINN::fia_env_dt)
  it  <- data.table::as.data.table(FINN::fia_init_trees)
  keep <- utils::head(sort(unique(obs$siteID)), 3L)
  list(obs = obs[siteID %in% keep], env = env[siteID %in% keep],
       it = it[siteID %in% keep], Nsp = max(obs$species))
}

reg_warn_fit <- function(d, data, family) {
  FINN.seed(1)
  m <- finn(N_species = d$Nsp, recruits_dbh = 12.9,
            competition_process  = createProcess(~0, FINN::competition),
            growth_process       = createProcess(~ temp, FINN::growth),
            regeneration_process = createProcess(~ temp, FINN::regeneration, optimizeSpecies = TRUE),
            mortality_process    = createProcess(~ temp, FINN::mortality))
  suppressMessages(fit(m, env = d$env, data = data,
                       init_cohort = makeInitCohorts(d$it, Nspecies = d$Nsp),
                       epochs = 1L, patches = 4, patch_size = 0.06, lr = 0,
                       loss_family = family, plot_progress = FALSE))
}

LOSS_REG <- c(dbh = "mse", ba = "mse", trees = "mse", growth = "mse",
              mortality = "mse", regeneration = "mse")

test_that("the term switched off warns", {
  skip_if_no_torch()
  d <- reg_warn_setup()
  expect_warning(reg_warn_fit(d, d$obs, replace(LOSS_REG, "regeneration", "none")),
                 "loss_family\\['regeneration'\\]")
})

test_that("a family with no finite observation warns too", {
  skip_if_no_torch()
  d <- reg_warn_setup()
  no_reg <- data.table::copy(d$obs)[, reg := NA_real_]
  expect_warning(reg_warn_fit(d, no_reg, LOSS_REG), "no finite `reg` value")
})

test_that("an observed recruitment does not warn", {
  skip_if_no_torch()
  d <- reg_warn_setup()
  expect_no_warning(reg_warn_fit(d, d$obs, LOSS_REG), message = "regeneration process is being calibrated")
})
