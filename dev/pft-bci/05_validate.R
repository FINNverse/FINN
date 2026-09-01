# ============================================================================
# Path B -- paper-style validation of the fitted models
# ============================================================================
# For every fitted condition (from its latest checkpoint model) computes:
#
#  (A) EQUILIBRIUM SIMULATION (as FINNetAl code/05_equilibrium.R): replay the env
#      ~20x with a small stochastic disturbance regime, run the forest to
#      equilibrium, and compare the simulated long-run pattern (mean ba/trees/dbh
#      per species) to the OBSERVED pattern -> does the model reproduce the
#      observed equilibrium composition/structure.
#
#  (B) HELD-OUT PREDICTION metrics on the test sites (from the saved predictions):
#      R2, RMSE and Spearman's rho per response variable (the paper reports rho +
#      RMSE; R2 added per request).
#
# Loads models from results/ckpt_<cond><tag>/epoch_<max>model.pt, so it works on
# partial checkpoints too. Writes results/validation_summary.rds + figures.
#
# Env vars: REPLAYS (default 20), EQ_PATCHES (default 10), EQ_LAST (years of the
# tail used as 'equilibrium', default 150).
# ============================================================================

local({ l <- file.path(path.expand("~"), "Rlib-finn020"); if (dir.exists(l)) .libPaths(c(l, .libPaths())) })
suppressWarnings(suppressMessages({ library(data.table); library(FINN); library(ggplot2) }))
source("dev/pft-bci/finn_membership.R")
D <- "dev/pft-bci"; RES <- file.path(D, "results"); set.seed(123)
REPLAYS   <- as.integer(Sys.getenv("REPLAYS", "20"))
EQ_PATCHES<- as.integer(Sys.getenv("EQ_PATCHES", "10"))
EQ_LAST   <- as.integer(Sys.getenv("EQ_LAST", "150"))
vars      <- c("dbh", "ba", "trees")

metrics <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs); pred <- pred[ok]; obs <- obs[ok]
  if (length(obs) < 3 || stats::sd(obs) == 0) return(data.table(R2=NA_real_, RMSE=NA_real_, rho=NA_real_, n=length(obs)))
  data.table(R2   = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
             RMSE = sqrt(mean((pred - obs)^2)),
             rho  = suppressWarnings(stats::cor(pred, obs, method = "spearman")),
             n    = length(obs))
}
latest_ckpt <- function(dir) {
  fs <- list.files(dir, pattern = "epoch_[0-9]+model.pt", full.names = TRUE)
  if (!length(fs)) return(NA_character_)
  fs[which.max(as.integer(gsub("\\D", "", basename(fs))))]
}
# data source for a condition name (annual vs _p7; species vs pft5)
data_paths <- function(nm) {
  p7 <- grepl("_p7$", nm); is5 <- grepl("^pft5", nm)
  if (is5) list(obs = file.path(D, if (p7) "data/pft5_p7/obs_dt.csv" else "data/pft5/obs_dt.csv"),
                env = file.path(D, if (p7) "data/pft5_p7/env_dt.csv" else "data/pft5/env_dt.csv"),
                coh = file.path(D, if (p7) "data/pft5_p7/initial_cohorts1985.csv" else "data/pft5/initial_cohorts1985.csv"))
  else     list(obs = file.path(D, if (p7) "data/obs_species_p7.csv" else "data/obs_species.csv"),
                env = file.path(D, if (p7) "data/env_p7.csv" else "data/env.csv"),
                coh = file.path(D, "data/initial_cohorts1985.csv"))
}

equilibrium_sim <- function(m, env, coh, Npatch) {
  # long env: 20 replays of the env block, contiguous increasing years
  ny <- max(env$year); env_long <- rbindlist(lapply(0:(REPLAYS-1), function(k) {
    e <- copy(env); e[, year := year + k*ny]; e }))
  # small stochastic disturbance regime (as paper): ~1%/yr, 0.4-1.6% patches
  d <- unique(env_long[, .(siteID, year)])
  d[, intensity := stats::rbinom(.N, 1, 1) * stats::runif(.N, 0.0043, 0.016)]
  ch <- FINN::CohortMat(obs_df = coh, sp = m$N_species)
  m$eval()
  sim <- m$simulate(env = env_long, disturbance = d, init_cohort = ch,
                    patches = Npatch, device = "cpu")
  w <- as.data.table(sim$wide$site)
  if ("reg" %in% names(w)) w[, reg := reg / 0.1]
  # 'equilibrium' = mean over the tail years, per species
  w[year > max(year) - EQ_LAST, lapply(.SD, mean, na.rm = TRUE),
    by = species, .SDcols = intersect(vars, names(w))]
}

res_files <- list.files(RES, pattern = "^pathB_.*\\.rds$", full.names = TRUE)
res_files <- res_files[!grepl("smoke|bench|verify|fixtest", res_files)]
cat("conditions found:", paste(sub("^pathB_", "", sub("\\.rds$", "", basename(res_files))), collapse = ", "), "\n\n")

heldout <- list(); equil <- list(); eqplot <- list()
for (rf in res_files) {
  nm <- sub("^pathB_", "", sub("\\.rds$", "", basename(rf)))
  R  <- readRDS(rf)

  ## (B) held-out prediction metrics
  pred <- as.data.table(R$pred)
  obs  <- melt(R$obs_te, id.vars = c("siteID","year","species"), measure.vars = vars,
               variable.name = "variable", value.name = "obs")[, variable := as.character(variable)]
  mB <- merge(pred[variable %in% vars], obs, by = c("siteID","year","species","variable"))
  hb <- mB[, metrics(value, obs), by = variable][, `:=`(cond = nm, kind = "heldout")]
  heldout[[nm]] <- hb

  ## (A) equilibrium simulation vs observed pattern
  ck <- latest_ckpt(file.path(RES, paste0("ckpt_", nm)))
  if (is.na(ck)) { cat(sprintf("[%s] no checkpoint -> skip equilibrium\n", nm)); next }
  dp <- data_paths(nm)
  m  <- tryCatch(torch::torch_load(ck), error = function(e) NULL)
  if (is.null(m)) { cat(sprintf("[%s] checkpoint load failed -> skip equilibrium\n", nm)); next }
  env <- fread(dp$env); coh <- fread(dp$coh); obsf <- fread(dp$obs)
  simeq <- tryCatch(equilibrium_sim(m, env, coh, EQ_PATCHES), error = function(e) { cat(sprintf("[%s] sim error: %s\n", nm, conditionMessage(e))); NULL })
  if (is.null(simeq)) next
  obseq <- obsf[, lapply(.SD, mean, na.rm = TRUE), by = species, .SDcols = vars]   # observed pattern
  cmp <- merge(melt(simeq, id.vars="species", variable.name="variable", value.name="sim"),
               melt(obseq, id.vars="species", variable.name="variable", value.name="obs"),
               by = c("species","variable"))
  ea <- cmp[, metrics(sim, obs), by = variable][, `:=`(cond = nm, kind = "equilibrium")]
  equil[[nm]] <- ea
  eqplot[[nm]] <- cmp[, cond := nm]
  cat(sprintf("[%s] epoch=%s | equilibrium ba: R2=%.2f rho=%.2f | held-out ba: R2=%.2f rmse=%.3g\n",
              nm, gsub("\\D","",basename(ck)),
              ea[variable=="ba"]$R2, ea[variable=="ba"]$rho,
              hb[variable=="ba"]$R2, hb[variable=="ba"]$RMSE))
}

summary_dt <- rbindlist(c(heldout, equil), fill = TRUE)
saveRDS(list(summary = summary_dt, equilibrium_points = rbindlist(eqplot, fill = TRUE)),
        file.path(RES, "validation_summary.rds"))

cat("\n==================== VALIDATION SUMMARY ====================\n")
cat("\n-- Equilibrium match (simulated long-run vs observed pattern) --\n")
print(dcast(summary_dt[kind=="equilibrium"], cond ~ variable, value.var = "R2"), digits = 2)
cat("\n-- Held-out prediction R2 --\n")
print(dcast(summary_dt[kind=="heldout"], cond ~ variable, value.var = "R2"), digits = 2)
cat("\n-- Held-out prediction RMSE --\n")
print(dcast(summary_dt[kind=="heldout"], cond ~ variable, value.var = "RMSE"), digits = 3)

## figure: simulated equilibrium vs observed (ba), per condition
pts <- rbindlist(eqplot, fill = TRUE)
if (nrow(pts)) {
  p <- ggplot(pts[variable=="ba" & is.finite(sim) & is.finite(obs)], aes(obs, sim)) +
    geom_abline(slope=1, intercept=0, colour="grey50", linewidth=0.3) +
    geom_point(alpha=0.4, size=0.8) + facet_wrap(~cond, scales="free") +
    scale_x_log10() + scale_y_log10() +
    labs(title="Equilibrium validation: simulated vs observed basal area per species",
         x="observed mean BA", y="simulated equilibrium BA") + theme_minimal(base_size=10)
  ggsave(file.path(D, "validation_equilibrium_ba.png"), p, width=9, height=6, dpi=130)
  cat("\nfigure: validation_equilibrium_ba.png\n")
}
