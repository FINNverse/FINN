# ============================================================================
# Diagnostic: is the ANNUAL simulation off BETWEEN the observed census years?
# ============================================================================
# The annual (period35) fit only sees observations at census years 5,10,...,30.
# This runs the fitted annual model over the full 30 annual steps and overlays
# the observed census points, to see whether the in-between years are sensible
# or spike (which would explain the equilibrium drift).
# Output: per-year stand trajectory (mean over sites) vs observed census points.
# ============================================================================
local({ l <- file.path(path.expand("~"), "Rlib-finn020"); if (dir.exists(l)) .libPaths(c(l, .libPaths())) })
suppressWarnings(suppressMessages({ library(data.table); library(FINN); library(ggplot2) }))
source("dev/pft-bci/finn_membership.R")
D <- "dev/pft-bci"; RES <- file.path(D, "results"); set.seed(1)

conds <- c(pft5 = "data/pft5", learned = "data")     # annual pft5 (5 sp) + learned (257 sp)
traj <- list(); obsl <- list()
for (nm in names(conds)) {
  ck <- file.path(RES, paste0("ckpt_", nm), "epoch_4000model.pt")
  if (!file.exists(ck)) { cat(nm, "no epoch_4000 ckpt\n"); next }
  m   <- torch::torch_load(ck)
  env <- fread(file.path(D, conds[[nm]], if (nm=="pft5") "env_dt.csv" else "env.csv"))
  coh <- fread(file.path(D, conds[[nm]], if (nm=="pft5") "initial_cohorts1985.csv" else "initial_cohorts1985.csv"))
  obs <- fread(file.path(D, if (nm=="pft5") "data/pft5/obs_dt.csv" else "data/obs_species.csv"))
  ch  <- FINN::CohortMat(obs_df = coh, sp = m$N_species); m$eval()
  sim <- predict(m, env = env, init_cohort = ch, patches = 10L, device = "cpu")
  w   <- as.data.table(sim$long$site); w[variable=="reg", value := value/0.1]
  # stand total per year (sum species) then mean over sites
  st  <- w[variable %in% c("ba","trees","reg"), .(value = sum(value, na.rm=TRUE)), by=.(siteID, year, variable)]
  st  <- st[, .(sim = mean(value, na.rm=TRUE)), by=.(year, variable)][, cond := nm]
  ob  <- melt(obs, id.vars=c("siteID","year","species"), measure.vars=c("ba","trees","reg"),
              variable.name="variable", value.name="value")
  ob  <- ob[, .(value=sum(value,na.rm=TRUE)), by=.(siteID,year,variable)][, .(obs=mean(value,na.rm=TRUE)), by=.(year,variable)][, cond := nm]
  traj[[nm]] <- st; obsl[[nm]] <- ob

  # spikiness: max between-census year vs the census-year values (trees, ba)
  cens <- sort(unique(obs$year))
  for (v in c("ba","trees")) {
    sv <- st[variable==v]; between <- sv[!year %in% cens]; atc <- sv[year %in% cens]
    cat(sprintf("[%s %s] max between-census sim=%.3g | max at-census sim=%.3g | ratio=%.2f | obs max=%.3g\n",
                nm, v, max(between$sim), max(atc$sim), max(between$sim)/max(atc$sim),
                max(ob[variable==v]$obs)))
  }
}
TR <- rbindlist(traj); OB <- rbindlist(obsl)
p <- ggplot(TR, aes(year, sim)) +
  geom_line(colour="steelblue") + geom_point(data=OB, aes(year, obs), colour="firebrick", size=1.6) +
  facet_grid(variable ~ cond, scales="free_y") +
  labs(title="Annual simulation (line) vs observed census points (red)",
       subtitle="stand total, mean over sites -- are the in-between years sensible?",
       y="value", x="year") + theme_minimal(base_size=11)
ggsave(file.path(D, "annual_trajectory_check.png"), p, width=8, height=6, dpi=130)
cat("\nfigure: annual_trajectory_check.png\n")
