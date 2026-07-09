# dev/precompute_vignettes.R  (STEP 3 of the data pipeline; NOT shipped)
# Precompute heavy (torch) results for the quarto vignettes, so the vignettes
# render on CRAN WITHOUT a libtorch backend. Inputs: the CSVs from
# dev/make_extdata.R (STEP 1) and the .pt pre-fits from dev/train_fia_model.R
# (STEP 2). See data-raw/README.md for the full chain. Run locally after
# retraining / changing a vignette's analysis, then reinstall the package so
# system.file() picks up the refreshed caches:
#
#   Rscript dev/precompute_vignettes.R
#   R CMD INSTALL --no-multiarch --no-docs --no-byte-compile .
#
# Each vignette loads its cache via
#   readRDS(system.file("extdata", "vig_<name>.rds", package = "FINN"))
# and the torch code itself is shown there with `#| eval: false`.

suppressMessages({library(FINN); library(torch); library(data.table)})
ext   <- function(f) system.file("extdata", f, package = "FINN")
outdir <- "inst/extdata"

# ---------------------------------------------------------------------------
# fit-to-fia
# ---------------------------------------------------------------------------
local({
  obs_dt     <- fread(ext("fia_obs_dt.csv"))
  env_dt     <- fread(ext("fia_env_dt.csv"))
  init_trees <- fread(ext("fia_init_trees.csv"))
  species_dt <- fread(ext("fia_species_dt.csv"))
  init_cohorts <- makeInitCohorts(init_trees, Nspecies = max(obs_dt$species))
  m        <- torch_load(ext("fia_process_finn.pt"))
  m_hybrid <- torch_load(ext("fia_hybrid_finn.pt"))

  # convergence
  loss <- data.table(epoch = seq_along(m$history), loss = sapply(m$history, sum))

  # assess (process)
  sim  <- predict(m, env = env_dt, init_cohort = init_cohorts,
                  patches = 4, patch_size = 0.06, device = "cpu")
  pred <- sim$long$site[, .(pred = mean(value)), by = .(siteID, year, species, variable)]
  pred[, variable := as.character(variable)]
  obs_m <- melt(obs_dt, id.vars = c("siteID", "year", "species", "species_name"),
                variable.name = "variable", value.name = "obs")
  obs_m[, variable := as.character(variable)]
  cmp <- merge(obs_m, pred, by = c("siteID", "year", "species", "variable"))
  fit_metrics <- cmp[!is.na(obs) & !is.na(pred),
                     .(n = .N, spearman = round(cor(obs, pred, method = "spearman"), 2)),
                     by = variable][order(-spearman)]

  # niches (process growth response to climate)
  scaling <- as.data.table(m$env_scaling); env_vars <- scaling$variable
  niche <- rbindlist(lapply(c("temp", "prec"), function(target) {
    s    <- scaling[variable == target]
    grid <- seq(min(env_dt[[target]]), max(env_dt[[target]]), length.out = 60)
    X <- matrix(0, length(grid), length(env_vars) + 1); X[, 1] <- 1
    X[, which(env_vars == target) + 1] <- (grid - s$center) / s$scale
    eff <- as.matrix(m$nn_growth(torch_tensor(X, dtype = torch_float32())))
    rbindlist(lapply(seq_len(ncol(eff)), function(j)
      data.table(var = target, species = j, x = grid, growth_mult = exp(eff[, j]))))
  }))
  niche <- merge(niche, species_dt, by = "species")

  # process vs hybrid performance
  metrics_of <- function(model) {
    model$eval()
    sim <- predict(model, env = env_dt, init_cohort = init_cohorts,
                   patches = 4, patch_size = 0.06, device = "cpu")
    p <- sim$long$site[, .(pred = mean(value)), by = .(siteID, year, species, variable)]
    p[, variable := as.character(variable)]
    cc <- merge(obs_m, p, by = c("siteID", "year", "species", "variable"))
    cc[!is.na(obs) & !is.na(pred),
       .(spearman = round(cor(obs, pred, method = "spearman"), 2),
         rmse     = round(sqrt(mean((obs - pred)^2)), 2)), by = variable]
  }
  comparison <- merge(
    metrics_of(m)[,        .(variable, spearman_proc = spearman, rmse_proc = rmse)],
    metrics_of(m_hybrid)[, .(variable, spearman_hyb  = spearman, rmse_hyb  = rmse)],
    by = "variable")[order(-spearman_proc)]

  # ALE growth response (process vs hybrid)
  FINN.seed(1)
  invisible(capture.output({
    ale_proc <- ALE(m,        env_dt, init_cohorts)
    ale_hyb  <- ALE(m_hybrid, env_dt, init_cohorts)
  }))
  growth <- rbind(
    data.table(ale_proc$growth)[, model := "Process (mechanistic)"],
    data.table(ale_hyb$growth)[,  model := "Hybrid (growth = NN)"])

  saveRDS(list(loss = loss, fit_metrics = fit_metrics, cmp = cmp, niche = niche,
               comparison = comparison, growth = growth),
          file.path(outdir, "vig_fit_to_fia.rds"), version = 2)
  cat("wrote vig_fit_to_fia.rds\n")
})

# ---------------------------------------------------------------------------
# succession-demo  (mirrors the parameter chunks of succession-demo.qmd)
# ---------------------------------------------------------------------------
local({
  FINN.seed(1234)
  Ntimesteps <- 300; Nsites <- 1; Nsp <- 4; patch_size <- 0.1
  shadeSP <- c(0.60, 0.42, 0.25, 0.08)
  parReg  <- shadeSP
  parRegEnv <- list(matrix(c(c(3.2, 2.6, 2.0, 1.2), rep(0, Nsp)), Nsp, 2))
  parGrowth <- matrix(c(shadeSP, c(0.040, 0.038, 0.034, 0.030)), Nsp, 2)
  parGrowthEnv <- list(matrix(c(c(1.00, 0.65, 0.30, -0.20), rep(0, Nsp)), Nsp, 2))
  parMort <- matrix(c(c(-1.0, -0.8, -0.5, -0.3), rep(0, Nsp), c(-0.3, -0.3, -0.3, -0.3)), Nsp, 3)
  parMortEnv <- list(matrix(c(c(-0.69, -1.40, -2.30, -3.30), rep(0, Nsp)), Nsp, 2))
  env_dt <- data.table(expand.grid(siteID = 1:Nsites, year = 1:Ntimesteps)); env_dt$env1 <- 0

  m <- finn(
    N_species            = Nsp,
    competition_process  = createProcess(~0,      FINN::competition),
    growth_process       = createProcess(~1+env1, FINN::growth,       initSpecies = parGrowth, initEnv = parGrowthEnv),
    regeneration_process = createProcess(~1+env1, FINN::regeneration, initSpecies = parReg,    initEnv = parRegEnv),
    mortality_process    = createProcess(~1+env1, FINN::mortality,    initSpecies = parMort,   initEnv = parMortEnv)
  )
  sim <- predict(m, init_cohort = NULL, env = env_dt, disturbance = NULL,
                 patches = 100, device = "cpu")
  p_dat <- sim$long$site[, .(value = mean(value)), by = .(year, species, variable)]
  p_dat[variable %in% c("ba", "trees"), value := value / patch_size]
  sp_labs <- c("1 pioneer", "2 early-mid", "3 mid-late", "4 climax")

  saveRDS(list(p_dat = p_dat, sp_labs = sp_labs),
          file.path(outdir, "vig_succession.rds"), version = 2)
  cat("wrote vig_succession.rds\n")
})

# ---------------------------------------------------------------------------
# Introduction-to-FINN  (mirrors the setup chunks of Introduction-to-FINN.qmd;
# parameters are random draws, so any valid seeded run makes a fine demo)
# ---------------------------------------------------------------------------
local({
  FINN.seed(1234)
  Ntimesteps <- 500; Nsites <- 1; patch_size <- 0.1; Nsp <- 5
  shadeSP <- c(0.1, 0.2, 0.5, 0.5, 0.7)
  parReg  <- shadeSP
  parRegEnv <- list(matrix(c(c(1, 2, 3, 3, 5), runif(Nsp, -2, 2)), Nsp, 2))
  parGrowth <- matrix(c(shadeSP, c(0.04, 0.05, 0.05, 0.06, 0.1)), Nsp, 2)
  parGrowthEnv <- list(matrix(c(c(0.2, 0.3, 0.5, 1, 1) * 0.5, runif(Nsp, -2, -0.5)), Nsp, 2))
  parMort <- matrix(c(as.numeric(scale(shadeSP)), as.numeric(scale(parGrowth[, 2])), rep(0, Nsp)), Nsp, 3)
  parMortEnv <- list(matrix(c(runif(Nsp, -3, -2), runif(Nsp, -3, -2)), Nsp, 2))

  env_dt <- data.table(expand.grid(list(siteID = 1:Nsites, year = 1:Ntimesteps)))
  dist_dt <- env_dt
  env_dt$env1 <- rep(0, Ntimesteps)
  disturbance_frequency <- 0.05
  disturbance_intensity <- rbinom(Ntimesteps * Nsites, 1, 0.2) * runif(Ntimesteps * Nsites, 0.5, 1)
  dist_dt$intensity <- rbinom(Ntimesteps * Nsites, 1, disturbance_frequency) * disturbance_intensity

  mk <- function() finn(N_species = Nsp,
    competition_process  = createProcess(~0, func = FINN::competition),
    growth_process       = createProcess(~1+env1, initEnv = parGrowthEnv, initSpecies = parGrowth, func = FINN::growth),
    regeneration_process = createProcess(~1+env1, initEnv = parRegEnv,    initSpecies = parReg,    func = FINN::regeneration),
    mortality_process    = createProcess(~1+env1, initEnv = parMortEnv,   initSpecies = parMort,   func = FINN::mortality))
  s1   <- simulateForest(mk(), init_cohort = NULL, env = env_dt, disturbance = dist_dt, device = "cpu", patches = 1)
  s100 <- simulateForest(mk(), init_cohort = NULL, env = env_dt, disturbance = dist_dt, device = "cpu", patches = 100)

  saveRDS(list(patches_1 = s1$long$site, patches_100 = s100$long$site),
          file.path(outdir, "vig_intro.rds"), version = 2)
  cat("wrote vig_intro.rds\n")
})
