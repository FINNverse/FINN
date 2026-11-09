## Fixture shared by test-identity.R and dev/make_identity_snapshots.R: four
## small fits that between them touch the constant-period path, the per-site
## path, the `nonfinite` switch and the per-tree aggregation.
Nsp_i <- 2; Nsites_i <- 3; Tmax_i <- 6; SEED_i <- 11
LOSS_i <- c(dbh = "mse", ba = "mse", trees = "mse", growth = "mse",
            mortality = "mse", regeneration = "mse")

mk_model_i <- function(sample_reg = FALSE) finn(N_species = Nsp_i, recruits_dbh = 5,
  competition_process  = createProcess(~0, func = FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~1 + env1, initEnv = matrix(c(0.5, 0.5, 0.4, -0.4), Nsp_i, 2),
                                       initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp_i, 2),
                                       func = FINN::growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(-2.5, -2.5, 0.3, -0.3), Nsp_i, 2),
                                       initSpecies = matrix(c(0.2, -0.2, 0.1, 0.1, 0, 0), Nsp_i, 3),
                                       func = FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(1.5, 1.0, 0.5, -0.5), Nsp_i, 2),
                                       initSpecies = c(0.1, 0.1), func = FINN::regeneration,
                                       sample_regeneration = sample_reg, optimizeSpecies = TRUE, optimizeEnv = TRUE))

identity_inputs <- function() {
  it <- data.table::data.table(
    siteID = rep(1:Nsites_i, each = 8), patchID = rep(rep(1:2, each = 4), Nsites_i),
    species = rep(1:2, 12), dbh = rep(c(10, 20, 30, 40), 6) + rep(0:2, each = 8) * 5,
    treeName = paste0("t", 1:24), trees = 3)
  env <- data.table::data.table(expand.grid(siteID = 1:Nsites_i, year = 1:Tmax_i))
  env[, env1 := (siteID - 2) * 0.5]
  obs <- data.table::data.table(expand.grid(siteID = 1:Nsites_i, year = c(3L, 6L), species = 1:Nsp_i))
  obs[, `:=`(dbh = 20, ba = 10, trees = 300, growth = 0.02, mort = 0.02, reg = 5,
             species_name = paste0("sp", species))]
  data.table::setorder(obs, siteID, year, species)
  tr <- data.table::data.table(expand.grid(siteID = 1:Nsites_i, year = c(3L, 6L),
                                           patchID = 1:2, slot = 1:4))
  tr[, `:=`(species = rep_len(1:Nsp_i, .N), dbh = rep_len(c(8, 16, 24, 32), .N),
            trees = 3, growth = 0.02, died = rep_len(c(0, 0, 1, 0), .N))]
  tr[, slot := NULL]
  list(init = makeInitCohorts(it, Nspecies = Nsp_i), env = env, obs = obs, tree = tr)
}

identity_configs <- function() list(
  constant   = list(period = 2L),
  per_site   = list(period = NULL),
  skip_step  = list(period = 2L, nonfinite = "skip_step"),
  tree_mort  = list(period = NULL, loss_aggregation = c(mortality = "tree"), tree_data = TRUE),
  ## the only config that draws recruits: guards the default regeneration sampler as well
  sample_reg = list(period = 2L, sample_reg = TRUE),
  ## unequal scored units per observation year. The other per-site configs carry the
  ## same number of cells in every year, so the per-unit pooling of the per-site loss
  ## and the former sum of per-year means agree there and the snapshot cannot see a
  ## change to the weighting (FINNverse/FINN#41). Here year 3 is one species-census
  ## short, so the weights are unequal and the snapshot moves with them.
  per_site_unbal = list(period = NULL, unbalanced = TRUE))

identity_run <- function(cfg) {
  spec <- identity_configs()[[cfg]]
  L <- identity_inputs()
  obs <- data.table::copy(L$obs)
  obs[, period_length := if (is.null(spec$period)) data.table::fifelse(siteID == 1L, 3L, 2L) else spec$period]
  # one species-census short in the first scored year: 5 cells against 6
  if (isTRUE(spec$unbalanced)) obs <- obs[!(siteID == 1L & year == 3L & species == 2L)]
  args <- list(env = L$env, data = obs, init_cohort = L$init, patches = 2L, patch_size = 0.1,
               env_autoscale = FALSE, epochs = 3L, lr = 0.01, weights = rep(1, 6),
               batchsize = Nsites_i, shuffle = FALSE, plot_progress = FALSE, device = "cpu",
               loss_family = LOSS_i)
  if (!is.null(spec$nonfinite)) args$nonfinite <- spec$nonfinite
  if (!is.null(spec$loss_aggregation)) args$loss_aggregation <- spec$loss_aggregation
  if (isTRUE(spec$tree_data)) args$tree_data <- L$tree
  FINN.seed(SEED_i)
  m <- mk_model_i(isTRUE(spec$sample_reg))
  suppressMessages(suppressWarnings(do.call(fit, c(list(m), args))))
  list(history = as.numeric(unlist(m$history)),
       pars    = unlist(lapply(m$parameters[order(names(m$parameters))],
                               function(p) as.numeric(p$detach()$cpu()))),
       weights = as.numeric(m$loss_weights), baseline = as.numeric(m$loss_baseline))
}
