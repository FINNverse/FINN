# Illustrative managed-vs-unmanaged demonstration of the management tool.
# NOTE: parameters are hand-set (NOT fitted) purely to produce a plausible stand
# so the mechanism is visible; the real calibration is the NW-FVA yield-table fit.
suppressMessages(devtools::load_all(".", quiet = TRUE)); suppressMessages(library(data.table)); suppressMessages(library(ggplot2))
FINN.seed(1)
H <- 120
env <- data.table(expand.grid(siteID = 1:2, year = 1:H)); env[, climate := 0]
# Hand-set to a realistic stand (dbh ~30->47 cm, unmanaged BA building to ~40 m2/ha).
m <- finn(N_species = 1L,
  competition_process  = createProcess(~0, FINN::competition,  initSpecies = matrix(c(0.25, 0.10), 1, 2)),
  growth_process       = createProcess(~1, FINN::growth,       initSpecies = matrix(c(0.05, 0.05), 1, 2), initEnv = list(matrix(-1.0, 1, 1))),
  regeneration_process = createProcess(~1, FINN::regeneration, initSpecies = matrix(0.05, 1, 1),          initEnv = list(matrix(1.6, 1, 1))),
  mortality_process    = createProcess(~1, FINN::mortality,    initSpecies = matrix(c(0,0,0), 1, 3),       initEnv = list(matrix(-3.2, 1, 1))))

ha <- function(sim, v) { d <- setDT(sim$long$site); d[variable == v, .(val = sum(value)/0.1/2), by = year] }  # per-ha, avg 2 sites

# unmanaged
base <- predict(m, init_cohort = NULL, env = env, patches = 30, patch_size = 0.1, device = "cpu")
b_ba <- ha(base, "ba"); setnames(b_ba, "val", "ba")

# managed: graduated thinning (below->above) every 10 yr from age 20, target 45 cm
sched <- data.frame(year = seq(20, 110, by = 10))
sched$prescription <- I(replicate(nrow(sched),
  management_params(thinning_intensity = 0.22,
                    thinning_size_bias = graduated_bias(from = -0.6, to = 0.8, ref = 52),
                    target_diameter = 52, target_harvest_fraction = 0.10,
                    target_diameter_spread = 4), simplify = FALSE))
run <- simulate_managed(m, env = env, patches = 30, schedule = sched, patch_size = 0.1, device = "cpu")
m_ba <- setDT(run$trajectory)[variable == "ba", .(ba = sum(value)/0.1/length(unique(siteID))), by = year]

cat("=== harvest log (site 1) ===\n"); print(head(run$harvest, 12))
cat(sprintf("\nunmanaged BA @120 = %.1f ; managed BA @120 = %.1f m2/ha\n",
            b_ba[year==H]$ba, m_ba[year==max(m_ba$year)]$ba))

df <- rbind(data.frame(year = b_ba$year, ba = b_ba$ba, regime = "unmanaged"),
            data.frame(year = m_ba$year, ba = m_ba$ba, regime = "managed (graduated thinning + Zielstaerkennutzung)"))
p <- ggplot(df, aes(year, ba, colour = regime)) +
  geom_vline(xintercept = seq(20,100,10), colour = "grey85", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("unmanaged" = "#1b9e77", "managed (graduated thinning + Zielstaerkennutzung)" = "#d95f02"), name = NULL) +
  labs(title = "FINN as a management simulator: managed vs unmanaged stand basal area",
       subtitle = "Illustrative (hand-set parameters). Grey lines = thinning entries. Driver: simulate_managed().",
       x = "stand age (yr)", y = expression("stand basal area (m"^2*" ha"^-1*")")) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("dev/pilots/demo_managed.png", p, width = 8, height = 4.6, dpi = 130)
cat("wrote dev/pilots/demo_managed.png\n")
