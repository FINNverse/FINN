# dev/train_fia_model.R  (NOT shipped; dev/ is .Rbuildignore'd)
# Trains the Process-FINN pre-fit shipped with the fit-to-fia vignette and saves
# it to inst/extdata/fia_process_finn.pt. Run once from the package root:
#   Rscript dev/train_fia_model.R [epochs]
suppressMessages({library(FINN); library(torch); library(data.table)})
FINN.seed(42)

epochs <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(epochs)) epochs <- 3000L

ext <- function(f) file.path("inst/extdata", f)
obs_dt  <- fread(ext("fia_obs_dt.csv"))
env_dt  <- fread(ext("fia_env_dt.csv"))
init    <- fread(ext("fia_init_trees.csv"))

init_cohorts <- makeInitCohorts(init, Nspecies = max(obs_dt$species))

m <- finn(
  N_species            = uniqueN(obs_dt$species),
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0, FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createProcess(~., FINN::growth,       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~., FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~., FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)

cat(sprintf("Training %d epochs on %d sites, %d species ...\n",
            epochs, uniqueN(obs_dt$siteID), uniqueN(obs_dt$species)))
t0 <- Sys.time()
m$fit(
  env         = env_dt,
  data        = obs_dt,
  init_cohort = init_cohorts,
  device      = "cpu",
  epochs      = epochs,
  batchsize   = 40L,
  patches     = 4,
  patch_size  = 0.06,
  weights     = c(0.1, 10, 1.0, 10.0, 1, 1),
  lr          = 0.01
)
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# drop bulky per-epoch training traces the vignette does not need (keep $history)
for (f in c("param_history", "loss_raw")) if (f %in% names(m)) m[[f]] <- NULL

torch::torch_save(m, ext("fia_process_finn.pt"))
cat("saved", ext("fia_process_finn.pt"),
    sprintf("(%.0f KB)\n", file.info(ext("fia_process_finn.pt"))$size / 1024))

# quick reload check (does the saved object come back usable?)
m2 <- torch::torch_load(ext("fia_process_finn.pt"))
cat("reload OK; classes:", paste(class(m2), collapse = ", "), "\n")
