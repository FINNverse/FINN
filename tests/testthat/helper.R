# Shared test helpers (sourced automatically by testthat before the test files).

# Torch guard: the libtorch backend is downloaded at runtime and is absent on
# CRAN / most CI, so torch-dependent tests skip there. Modelled on cito's
# skip_if_no_torch().
skip_if_no_torch <- function() {
  testthat::skip_if_not_installed("torch")
  testthat::skip_if_not(torch::torch_is_installed(),
                        "libtorch backend not available")
}

# A tiny, hand-checkable tree list for the data-preparation functions. One site,
# one patch, two inventories (2000 and 2010, a 10-year interval), four trees:
#   t1 oak  : alive both visits, dbh 20 -> 24         (grows; survives)
#   t2 oak  : alive 2000 (dbh 30), dead 2010          (at risk; dies)
#   t3 pine : alive both visits, dbh 40 -> 45         (grows; survives)
#   t4 pine : recruit ("new") in 2010, dbh 12.9       (regeneration; not at risk)
# These give exactly predictable growth / mortality / regeneration values.
toy_tree_dt <- function() {
  data.table::data.table(
    siteName     = "s1",
    patchName    = "p1",
    treeName     = c("t1", "t1", "t2", "t2", "t3", "t3", "t4"),
    year         = c(2000, 2010, 2000, 2010, 2000, 2010, 2010),
    species_name = c("oak", "oak", "oak", "oak", "pine", "pine", "pine"),
    dbh          = c(20,   24,   30,   NA,   40,   45,   12.9),
    status       = c("alive","alive","alive","dead","alive","alive","new"),
    living       = c(TRUE,  TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
  )
}

# The bundled FIA tables on a THREE-STEP axis, one step per inventory: env
# years 0, 1, 2 and the two remeasurements at years 1 and 2, without
# `period_length`. The datasets themselves are annual (env 1..22, observations
# at years 10 and 20 with period_length 10); this compressed copy is only for the
# tests of the engine's step bookkeeping (anchors, windows, observed recruits,
# input checks), which name individual timesteps and are about the mechanics,
# not about the data. The climate is static, so any three rows of env are the
# same climate.
fia_census_steps <- function(sites = 1:4) {
  obs <- data.table::copy(FINN::fia_obs_dt[siteID %in% sites])
  obs[, year := data.table::frank(year, ties.method = "dense"), by = siteID]
  obs[, period_length := NULL]
  env <- data.table::copy(FINN::fia_env_dt[siteID %in% sites & year %in% 1:3])
  env[, year := year - 1L]
  list(obs = obs, env = env, it = data.table::copy(FINN::fia_init_trees[siteID %in% sites]),
       Nsp = max(FINN::fia_obs_dt$species))
}

# Fit a tiny FINN model on a 4-site subset of the bundled FIA datasets (the real
# input format). Torch-gated: callers MUST skip_if_no_torch() first. Returns the
# fitted model together with the env / init_cohort it was fit on, so the
# integration test and the xAI smoke tests can share one small, ~1s setup.
fit_toy_finn <- function(epochs = 2L, n_sites = 4L) {
  obs <- FINN::fia_obs_dt
  env <- FINN::fia_env_dt
  it  <- FINN::fia_init_trees

  sites <- utils::head(sort(unique(obs$siteID)), n_sites)
  obs <- obs[siteID %in% sites]
  env <- env[siteID %in% sites]
  it  <- it[siteID %in% sites]
  Nsp <- max(obs$species)

  FINN.seed(1)
  init <- makeInitCohorts(it, Nspecies = Nsp)
  m <- finn(
    N_species            = Nsp, recruits_dbh = 12.9,
    competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    growth_process       = createProcess(~ temp + prec, FINN::growth,       optimizeSpecies = TRUE, optimizeEnv = TRUE),
    regeneration_process = createProcess(~ temp + prec, FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = createProcess(~ temp + prec, FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
  )
  suppressMessages(
    fit(m, env = env, data = obs, init_cohort = init, device = "cpu",
        epochs = epochs, patches = 4, patch_size = 0.06, lr = 0.01, plot_progress = FALSE)
  )
  list(model = m, env = env, init = init, obs = obs, Nsp = Nsp)
}
