# =============================================================================
# End-to-end: FIT FINN (through fit(), not a regression) to the NW-FVA Fichte
# yield table, then simulate and check it reproduces Dg(age) and N(age).
#
# Setup: each yield class = one "site"; site index (Ekl) is the environment.
# growth_env (site-responsive size-decline) + mortality (free, to absorb the
# thinning-driven density decline the table folds into removal); regeneration
# suppressed (even-aged). We report whether the simulated stand tracks the table.
# =============================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages({library(data.table); library(ggplot2)})
FINN.seed(1)

gy <- as.data.table(nwfva_yield_tables())[Art == "Fichte"][order(Ekl, Alter)]
classes <- sort(unique(gy$Ekl))
patch_size <- 1.0                      # work in per-ha units directly

obs <- list(); init <- list(); envm <- list()
for (i in seq_along(classes)) {
  d <- gy[Ekl == classes[i]]; y0 <- d$Alter[1]; d[, yr := Alter - y0]
  obs[[i]]  <- data.table(siteID = i, year = d$yr[-1], species = 1L,
                          dbh = d$Dg[-1], ba = d$G[-1], trees = d$N[-1],
                          growth = NA_real_, mort = NA_real_, reg = NA_real_)
  init[[i]] <- data.table(siteID = i, patchID = 1L, species = 1L,
                          dbh = d$Dg[1], trees = d$N[1])
  envm[[i]] <- data.table(siteID = i, ekl = as.numeric(classes[i]))
}
obs_dt <- rbindlist(obs)
Tmax   <- max(obs_dt$year)
env_dt <- merge(CJ(siteID = seq_along(classes), year = 1:Tmax),
                rbindlist(envm), by = "siteID")
init_cohort <- makeInitCohorts(rbindlist(init), Nspecies = 1L)

# --- model: growth responds to site (growth_env); regeneration switched off ----
m <- finn(
  N_species            = 1L,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~ ekl, FINN::growth_env, optimizeSpecies = TRUE,
                                       optimizeEnv = TRUE, custom_parameters = list(k_env = 0)),
  mortality_process    = createProcess(~ ekl, FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~1, FINN::regeneration, initEnv = list(matrix(-12, 1, 1)),
                                       optimizeSpecies = FALSE, optimizeEnv = FALSE)  # ~0 recruits
)

fit(m, data = obs_dt, env = env_dt, init_cohort = init_cohort,
    patches = 1, patch_size = patch_size, device = "cpu",
    epochs = 400L, lr = 0.02, batchsize = length(classes),
    weights = c(dbh = 5, ba = 1, trees = 1, growth = 0, mortality = 0, regeneration = 0),
    env_autoscale = TRUE, plot_progress = FALSE)

# --- simulate the fitted model and compare to the table ------------------------
sim <- predict(m, init_cohort = init_cohort, env = env_dt, patches = 1,
               patch_size = patch_size, device = "cpu")
s <- as.data.table(sim$long$site)
sim_dbh <- s[variable == "dbh", .(siteID, year, sim = value)]
sim_n   <- s[variable == "trees", .(siteID, year, sim = value)]

cmp <- copy(obs_dt)[, .(siteID, year, Dg_obs = dbh, N_obs = trees)]
cmp <- merge(cmp, sim_dbh, by = c("siteID","year")); setnames(cmp, "sim", "Dg_sim")
cmp <- merge(cmp, sim_n,   by = c("siteID","year")); setnames(cmp, "sim", "N_sim")
cmp[, Ekl := classes[siteID]][, age := year + gy[, .(a = min(Alter)), Ekl][match(cmp$Ekl, Ekl)]$a]

r2 <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
cat(sprintf("\n=== FIT QUALITY ===\n  Dg  R2 = %.3f\n  N   R2 = %.3f\n",
            r2(cmp$Dg_obs, cmp$Dg_sim), r2(cmp$N_obs, cmp$N_sim)))

pl <- melt(cmp, id.vars = c("Ekl","age"), measure.vars = c("Dg_obs","Dg_sim"))
pl[, kind := ifelse(grepl("obs", variable), "yield table", "FINN (fitted)")]
p <- ggplot(pl, aes(age, value, colour = factor(Ekl))) +
  geom_point(data = pl[kind == "yield table"], size = 1, alpha = .8) +
  geom_line(data = pl[kind == "FINN (fitted)"], linewidth = .7) +
  scale_colour_viridis_d(option = "D", end = .9, name = "Ertragsklasse") +
  labs(title = "FINN fitted to the NW-FVA Fichte yield table",
       subtitle = sprintf("points = table Dg(age); lines = fitted FINN simulation (Dg R² = %.3f)",
                          r2(cmp$Dg_obs, cmp$Dg_sim)),
       x = "Alter (yr)", y = expression(D[g]~"(cm)")) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("dev/pilots/fit_yieldtable.png", p, width = 8, height = 5, dpi = 130)
cat("wrote dev/pilots/fit_yieldtable.png\n")
