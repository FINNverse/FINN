# ============================================================================
# Path B -- fit ONE condition of the BCI soft-membership experiment
# ============================================================================
# Runs one of three species-resolution FINN models on BCI, all via
# finn_membership so the ONLY difference is how per-species demographic
# parameters are grouped:
#
#   free    : frozen identity membership  -> free per-species params  (baseline)
#   ruger   : frozen Rueger PFT one-hot   -> params tied by Rueger's 5 PFTs
#   learned : learned soft membership K=5 -> FINN discovers the grouping
#
# Held out: a fixed set of SITES, so we can measure prediction on species where
# training data is scarce (the rare-species gain). One condition per run -> maps
# to a cluster array task. Writes a self-contained result RDS.
#
# Config via env vars (cluster overrides): COND, EPOCHS, DEVICE, PATCHES, TAG.
# ============================================================================

# prefer the cluster's FINN 0.2.0 lib if present (no-op on other machines)
local({ l <- file.path(path.expand("~"), "Rlib-finn020"); if (dir.exists(l)) .libPaths(c(l, .libPaths())) })
suppressMessages({ library(FINN); library(torch); library(data.table) })
source("dev/pft-bci/finn_membership.R")

COND    <- Sys.getenv("COND",    "learned")               # free | ruger | learned
EPOCHS  <- as.integer(Sys.getenv("EPOCHS",  "8000"))
DEVICE  <- Sys.getenv("DEVICE",  "cpu")
PATCHES <- as.integer(Sys.getenv("PATCHES", "25"))
LR      <- as.numeric(Sys.getenv("LR", "0.01"))
BATCH   <- as.integer(Sys.getenv("BATCH", "250"))
SAT     <- as.integer(Sys.getenv("SAT", "1"))            # 1 = regeneration_saturation ON (annual-timestep fix)
TAG     <- Sys.getenv("TAG", "")
SEED    <- as.integer(Sys.getenv("SEED", "42"))
D       <- "dev/pft-bci"; RES <- file.path(D, "results"); dir.create(RES, showWarnings = FALSE)
set.seed(SEED); FINN.seed(SEED)

# pft5 = the PUBLISHED approach (data aggregated to 5 Rueger PFTs), refit with the
# CURRENT FINN as a same-version anchor. All other conditions use species-resolution.
if (COND == "pft5") {
  FE   <- file.path(D, "data/pft5")               # bundled published 5-PFT data
  obs  <- fread(file.path(FE, "obs_dt.csv"))
  env  <- fread(file.path(FE, "env_dt.csv"))
  coh  <- fread(file.path(FE, "initial_cohorts1985.csv"))
  pftm <- data.table(species = 1:5, PFT_2axes = 1:5)
} else {
  obs  <- fread(file.path(D, "data/obs_species.csv"))
  env  <- fread(file.path(D, "data/env.csv"))
  coh  <- fread(file.path(D, "data/initial_cohorts1985.csv"))
  pftm <- fread(file.path(D, "data/species_pft.csv"))
}
Nsp  <- max(obs$species)
predictors <- c("Prec","SR_kW_m2","RH_prc","T_max","T_min","swp")

# ---- held-out sites --------------------------------------------------------
all_sites  <- sort(unique(obs$siteID))
test_sites <- all_sites[seq(2, length(all_sites), by = 4)]   # 5 of 20, deterministic
train_sites <- setdiff(all_sites, test_sites)
sub <- function(dt, sites) dt[siteID %in% sites]
reindex <- function(dt) { dt <- copy(dt); dt[, siteID := as.integer(as.factor(siteID))]; dt }
obs_tr <- reindex(sub(obs, train_sites)); env_tr <- reindex(sub(env, train_sites)); coh_tr <- reindex(sub(coh, train_sites))
obs_te <- reindex(sub(obs, test_sites));  env_te <- reindex(sub(env, test_sites));  coh_te <- reindex(sub(coh, test_sites))

# ---- Rueger membership vector (length Nsp); "other"/unclassified -> own group
ruger <- rep(max(pftm$PFT_2axes) + 1L, Nsp)
ruger[pftm$species] <- pftm$PFT_2axes

# ---- build the model for this condition ------------------------------------
reg_sat <- if (SAT) list(K_init = 200, shared = TRUE, bounds = c(1, 3000)) else NULL
mk <- function(membership = NULL, K = 5L) finn_membership(
  mm_K = K, mm_membership = membership, N_species = Nsp, recruits_dbh = 1.0,
  regeneration_saturation = reg_sat,
  competition_process  = createProcess(~0, FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createProcess(~., FINN::growth,       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~., FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~., FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE))
m <- switch(COND,
  free    = mk(membership = seq_len(Nsp)),          # frozen identity -> free per-species
  ruger   = mk(membership = ruger),                 # frozen Rueger PFT one-hot
  learned = mk(K = 5L),                             # learned soft membership
  pft5    = mk(membership = seq_len(Nsp)),          # 5 PFT-units, free (published anchor)
  stop("unknown COND: ", COND))
cat("condition:", COND, " | N_species:", Nsp, " | K:", m$mm_K,
    " | train sites:", length(train_sites), " test:", length(test_sites),
    " | epochs:", EPOCHS, " | device:", DEVICE, "\n")

cohort_tr <- FINN::CohortMat$new(obs_df = coh_tr, sp = Nsp)
# FINN's dataloader batches over SITES and drops the last incomplete batch, so
# batchsize must not exceed the number of training sites (else zero batches).
BATCH <- min(BATCH, length(train_sites))
t0 <- proc.time()[3]
fit(m, data = obs_tr, env = env_tr, init_cohort = cohort_tr, device = DEVICE,
    epochs = EPOCHS, batchsize = BATCH, patches = PATCHES, lr = LR,
    optimizer = torch::optim_adam,                      # cluster torch lacks optim_ignite_adam
    env_autoscale = TRUE, plot_progress = FALSE,
    weights = c(0.1, 10, 1.0, 1, 1, 1),                 # annual-timestep working config
    loss = c(dbh="mse", ba="mse", trees="nbinom", growth="mse", mortality="mse", regeneration="nbinom"))
cat("fit done in", round(proc.time()[3] - t0, 1), "s\n")

# ---- held-out prediction on test sites -------------------------------------
cohort_te <- FINN::CohortMat$new(obs_df = coh_te, sp = Nsp)
m$eval()
sim <- predict(m, env = env_te, init_cohort = cohort_te, patches = PATCHES, device = DEVICE)
pred <- sim$long$site
setDT(pred)

# training observation count per species (rare-species definition)
nobs <- obs_tr[!is.na(growth) & trees > 0, .(n_train = .N), by = species]

out <- list(
  cond = COND, K = m$mm_K, epochs = EPOCHS, seed = SEED,
  test_sites = test_sites, train_sites = train_sites,
  A = as.matrix(m$mm_A), ruger = ruger, nobs = nobs,
  pred = pred, obs_te = obs_te,
  history = tryCatch(m$history, error = function(e) NULL))
saveRDS(out, file.path(RES, paste0("pathB_", COND, TAG, ".rds")))
cat("saved", file.path(RES, paste0("pathB_", COND, TAG, ".rds")), "\n")
