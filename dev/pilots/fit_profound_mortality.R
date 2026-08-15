# =============================================================================
# Fit FINN's NATURAL mortality to PROFOUND (Bily Kriz): an unmanaged Norway-spruce
# monitoring stand, 1997-2015, stems 602 -> 313 by self-thinning. This is the
# natural-death signal the yield tables cannot give (they fold death into harvest).
# We fit growth + mortality (regeneration off) and check FINN reproduces N(year)
# and Dg(year).
# =============================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages({library(data.table); library(ggplot2)})

b  <- readRDS(path.expand("~/working-directory/archive/FINN-PROFOUND/raw-data/bily_kriz.rds"))
tr <- as.data.table(b$trees)[!is.na(dbh1_cm)]
agg <- tr[, .(N = .N, Dg = sqrt(mean(dbh1_cm^2)), G = sum(pi * (dbh1_cm/200)^2)), by = year][order(year)]
agg[, yr := year - min(year)]                                   # 0..18

# init cohorts = the 1997 diameter distribution (2 cm bins), duplicated to 2 sites
# (FINN squeezes a single site; two identical replicates sidestep it).
tr97 <- tr[year == min(year)][, dbh := round(dbh1_cm / 2) * 2]
ci <- tr97[, .(trees = .N), by = dbh][order(dbh)][, `:=`(species = 1L)]
init_trees <- rbindlist(lapply(1:2, function(s) copy(ci)[, `:=`(siteID = s, patchID = 1L)]))
init_trees[, cohortID := seq_len(.N), by = .(siteID, patchID)]
init_cohort <- makeInitCohorts(init_trees, Nspecies = 1L)

obs_dt <- rbindlist(lapply(1:2, function(s)
  data.table(siteID = s, year = agg$yr[-1], species = 1L,
             dbh = agg$Dg[-1], ba = agg$G[-1], trees = agg$N[-1],
             growth = NA_real_, mort = NA_real_, reg = NA_real_)))
Tmax <- max(obs_dt$year)
env_dt <- CJ(siteID = 1:2, year = 1:Tmax)[, const := 0]

FINN.seed(1)
m <- finn(N_species = 1L,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~1, FINN::growth,      optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~1, FINN::mortality,   optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~1, FINN::regeneration, initEnv = list(matrix(-12, 1, 1)),
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE))
fit(m, data = obs_dt, env = env_dt, init_cohort = init_cohort,
    patches = 1, patch_size = 1.0, device = "cpu", epochs = 500L, lr = 0.02, batchsize = 2L,
    weights = c(dbh = 3, ba = 1, trees = 5, growth = 0, mortality = 0, regeneration = 0),
    env_autoscale = TRUE, plot_progress = FALSE)

sim <- predict(m, init_cohort = init_cohort, env = env_dt, patches = 1, patch_size = 1.0, device = "cpu")
s <- as.data.table(sim$long$site)[siteID == 1]
simN  <- s[variable == "trees", .(yr = year, N_sim = value)]
simDg <- s[variable == "dbh",   .(yr = year, Dg_sim = value)]
cmp <- merge(agg[yr > 0, .(yr, N_obs = N, Dg_obs = Dg)], simN, by = "yr")
cmp <- merge(cmp, simDg, by = "yr")[, age := yr + min(agg$year)]
r2 <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
cat(sprintf("\n=== PROFOUND Bily Kriz natural-mortality fit ===\n  N(stems) R2 = %.3f\n  Dg      R2 = %.3f\n",
            r2(cmp$N_obs, cmp$N_sim), r2(cmp$Dg_obs, cmp$Dg_sim)))
print(cmp[, .(year = age, N_obs, N_sim = round(N_sim), Dg_obs = round(Dg_obs,1), Dg_sim = round(Dg_sim,1))])

pl <- rbind(
  data.table(year = agg$year, val = agg$N, kind = "observed", var = "N (stems)"),
  data.table(year = cmp$age,  val = cmp$N_sim, kind = "FINN", var = "N (stems)"),
  data.table(year = agg$year, val = agg$Dg, kind = "observed", var = "Dg (cm)"),
  data.table(year = cmp$age,  val = cmp$Dg_sim, kind = "FINN", var = "Dg (cm)"))
p <- ggplot(pl, aes(year, val, colour = kind)) +
  geom_point(data = pl[kind == "observed"], size = 1.3) +
  geom_line(data = pl[kind == "FINN"], linewidth = .8) +
  facet_wrap(~var, scales = "free_y") +
  scale_colour_manual(values = c(observed = "#1c6b4a", FINN = "#b25a1c"), name = NULL) +
  labs(title = "FINN's natural mortality fitted to PROFOUND (Bily Kriz spruce)",
       subtitle = "unmanaged self-thinning 1997-2015; points = census, lines = fitted FINN",
       x = "year", y = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("dev/pilots/fit_profound_mortality.png", p, width = 8, height = 4, dpi = 130)
cat("wrote dev/pilots/fit_profound_mortality.png\n")
