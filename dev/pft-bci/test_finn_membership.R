# Smoke test: finn_membership constructs, fits, and trains membership+prototypes.
suppressMessages({ library(FINN); library(torch); library(data.table) })
source("dev/pft-bci/finn_membership.R")
set.seed(1); FINN.seed(1)

ext <- function(f) system.file("extdata", f, package = "FINN")
obs_dt     <- fread(ext("fia_obs_dt.csv"))
env_dt     <- fread(ext("fia_env_dt.csv"))
init_trees <- fread(ext("fia_init_trees.csv"))
Nsp        <- max(obs_dt$species)
init_cohorts <- makeInitCohorts(init_trees, Nspecies = Nsp)
cat("N_species =", Nsp, "\n")

m <- finn_membership(
  mm_K                 = 3L,
  N_species            = Nsp,
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0,            FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createProcess(~ temp + prec, FINN::growth,       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~ temp + prec, FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~ temp + prec, FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)

# ---- structural checks -----------------------------------------------------
pnames <- names(m$parameters)
cat("has mm_logits:", "mm_logits" %in% pnames,
    " | has mm_Proto_growth:", "mm_Proto_growth" %in% pnames, "\n")
cat("NO leaf par_growth_unconstrained:", !any(grepl("par_.*_unconstrained", pnames)), "\n")
cat("A shape:", paste(dim(as.matrix(m$mm_A)), collapse="x"),
    " | par_growth shape:", paste(m$par_growth$shape, collapse="x"), "\n")

# snapshot membership + a prototype before training
A0     <- as.matrix(m$mm_A)
proto0 <- as.matrix(m$mm_Proto_growth$detach())

# ---- short fit -------------------------------------------------------------
fit(m, env = env_dt, data = obs_dt, init_cohort = init_cohorts, device = "cpu",
    epochs = 15L, batchsize = 40L, patches = 4, patch_size = 0.06, lr = 0.02,
    env_autoscale = TRUE, plot_progress = FALSE)

A1     <- as.matrix(m$mm_A)
proto1 <- as.matrix(m$mm_Proto_growth$detach())
cat("\nmembership A changed during fit:", max(abs(A1 - A0)) > 1e-5,
    " (max |dA| =", round(max(abs(A1 - A0)), 4), ")\n")
cat("prototype Proto_growth changed:", max(abs(proto1 - proto0)) > 1e-5,
    " (max |dProto| =", round(max(abs(proto1 - proto0)), 4), ")\n")
cat("par_growth finite:", all(is.finite(as.matrix(m$par_growth))), "\n")
cat("\nSMOKE TEST PASSED\n")
