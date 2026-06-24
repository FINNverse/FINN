# dev/train_fia_model.R  (STEP 2 of the data pipeline; NOT shipped — dev/ is .Rbuildignore'd)
# Reads the bundled CSVs made by dev/make_extdata.R (STEP 1) and trains the two
# pre-fits shipped with the fit-to-fia vignette, saving them to inst/extdata/:
#   - fia_process_finn.pt : Process-FINN (mechanistic growth)
#   - fia_hybrid_finn.pt  : Hybrid-FINN  (growth replaced by a neural network)
# Downstream: dev/precompute_vignettes.R (STEP 3) caches the vignette results.
# See data-raw/README.md for the full chain. Run from the package root (needs
# a torch backend):  Rscript dev/train_fia_model.R [epochs]
suppressMessages({library(FINN); library(torch); library(data.table)})

epochs <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(epochs)) epochs <- 3000L

ext <- function(f) file.path("inst/extdata", f)
obs_dt  <- fread(ext("fia_obs_dt.csv"))
env_dt  <- fread(ext("fia_env_dt.csv"))      # RAW (untransformed) environment
init    <- fread(ext("fia_init_trees.csv"))
init_cohorts <- makeInitCohorts(init, Nspecies = max(obs_dt$species))

train_and_save <- function(label, growth_process, file) {
  FINN.seed(42)
  m <- finn(
    N_species            = uniqueN(obs_dt$species),
    recruits_dbh         = 12.9,
    competition_process  = createProcess(~0, FINN::competition,  optimizeSpecies = TRUE),
    growth_process       = growth_process,
    regeneration_process = createProcess(~., FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = createProcess(~., FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
  )
  cat(sprintf("[%s] training %d epochs on %d sites, %d species ...\n",
              label, epochs, uniqueN(obs_dt$siteID), uniqueN(obs_dt$species)))
  t0 <- Sys.time()
  m$fit(
    env = env_dt, data = obs_dt, init_cohort = init_cohorts, device = "cpu",
    epochs = epochs, batchsize = 40L, patches = 4, patch_size = 0.06,
    weights = c(0.1, 10, 1.0, 10.0, 1, 1), lr = 0.01,
    env_autoscale = TRUE       # FINN learns + stores the standardization (default)
  )
  cat(sprintf("[%s] done in %.1f min\n", label,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  for (f in c("param_history", "loss_raw")) if (f %in% names(m)) m[[f]] <- NULL
  torch::torch_save(m, ext(file))
  m2 <- torch::torch_load(ext(file))   # reload check
  cat(sprintf("[%s] saved %s (%.0f KB); reload class: %s\n", label, file,
              file.info(ext(file))$size / 1024, paste(class(m2), collapse = ", ")))
}

# 1) Process-FINN: mechanistic growth
train_and_save(
  "process",
  createProcess(~., FINN::growth, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  "fia_process_finn.pt"
)

# 2) Hybrid-FINN: growth replaced by a small (plain MLP) neural network.
#    transformer = FALSE keeps the shipped model small (~2.4 MB vs ~4.3 MB).
train_and_save(
  "hybrid",
  createHybrid(~., hidden = c(20L, 20L), transformer = FALSE),
  "fia_hybrid_finn.pt"
)
