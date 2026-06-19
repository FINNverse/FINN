# Scratch: tune 4-species succession parameters for vignettes/succession-demo.qmd
# NOT shipped (dev/ is .Rbuildignore'd). Goal: pioneers peak early, climax late.
#
# KEY MECHANIC (verified in R/Processes.R):
#   growth/regeneration gate on `light` via sigmoid centered at the species'
#   shade parameter. A species only grows/regenerates where light > shade_param.
#   => LOW shade param  = shade-tolerant = CLIMAX (works under closed canopy)
#      HIGH shade param = light-demanding = PIONEER (needs open ground)
#   (The plan's table has this inverted; we follow the source.)
library(FINN)
library(data.table)
library(ggplot2)

FINN.seed(1234)
Ntimesteps <- 300
Nsites     <- 1
Nsp        <- 4          # pioneer -> early-mid -> mid-late -> climax
patch_size <- 0.1

# ---- species trade-offs (sp1 pioneer ... sp4 climax) -----------------------
# shade = light fraction a species needs to grow/regenerate (pioneer high)
shadeSP <- c(0.60, 0.42, 0.25, 0.08)

# regeneration light threshold tracks shade tolerance; fecundity via regEnv int.
parReg    <- shadeSP
parRegEnv <- list(matrix(c(
  c(3.2, 2.6, 2.0, 1.2),   # intercept: seed rain, pioneer prolific, climax sparse
  rep(0, Nsp)              # env slope = 0 (constant-env demo)
), Nsp, 2))

# growth: col1 = shade threshold, col2 = size-dependent growth decay.
#   Keep col2 gentle & near-equal so it doesn't accidentally handicap pioneers;
#   the speed trade-off is carried by the growth-rate intercept below.
parGrowth <- matrix(c(
  shadeSP,
  c(0.040, 0.038, 0.034, 0.030)
), Nsp, 2)
parGrowthEnv <- list(matrix(c(
  c(1.00, 0.65, 0.30, -0.20),  # intercept: overall growth rate (pioneer ~4x climax);
  rep(0, Nsp)                   # scaled down vs before to slow canopy closure
), Nsp, 2))                     # -> stretches the early-seral phase

# mortality varies through the run via two ecological dependences:
#   col1 (light): negative -> mortality RISES as a cohort is overtopped (low
#                 light). Steep for pioneers (shade-intolerant), ~flat for climax.
#   col3 (growth): negative -> suppressed, slow-growing trees die more (vigour).
#   col2 (dbh): left 0.
# mort = sigmoid(intercept + col1*light + col3*growth)
parMort <- matrix(c(
  c(-1.0, -0.8, -0.5, -0.3),   # col1 light: pioneer strongly shade-intolerant
  rep(0, Nsp),                 # col2 dbh: off
  c(-0.3, -0.3, -0.3, -0.3)    # col3 growth: vigour lowers mortality
), Nsp, 3)
parMortEnv <- list(matrix(c(
  c(-0.69, -1.40, -2.30, -3.30),  # intercept: per-species baseline level
  rep(0, Nsp)
), Nsp, 2))

# competition: col1 = height allometry (climax tallest -> overtops late),
#              col2 = competition strength
parComp <- matrix(c(
  c(0.50, 0.55, 0.62, 0.70),
  c(0.30, 0.25, 0.20, 0.15)
), Nsp, 2)

# ---- constant environment, no disturbance (clean primary succession) -------
env_dt <- data.table(expand.grid(siteID = 1:Nsites, year = 1:Ntimesteps))
env_dt$env1 <- 0

# ---- model + simulate ------------------------------------------------------
m <- finn(
  N_species            = Nsp,
  competition_process  = createProcess(~0,      FINN::competition),
  growth_process       = createProcess(~1+env1, FINN::growth,       initSpecies = parGrowth, initEnv = parGrowthEnv),
  regeneration_process = createProcess(~1+env1, FINN::regeneration, initSpecies = parReg,    initEnv = parRegEnv),
  mortality_process    = createProcess(~1+env1, FINN::mortality,    initSpecies = parMort,   initEnv = parMortEnv)
)

sim <- m$simulate(init_cohort = NULL, env = env_dt, disturbance = NULL,
                  patches = 100, device = "cpu")

p_dat <- sim$long$site[, .(value = mean(value)), by = .(year, species, variable)]
p_dat[variable %in% c("ba","trees"), value := value / patch_size]

ba <- p_dat[variable == "ba"]
cat("\n=== BA peak year per species ===\n")
print(ba[, .(peak_year = year[which.max(value)], peak_ba = round(max(value),2),
             final_ba = round(value[year == max(year)],2)), by = species][order(species)])

cat("\n=== mortality range per species (should now VARY) ===\n")
print(p_dat[variable == "mort", .(mort_min = round(min(value),3), mort_max = round(max(value),3),
             mort_mean = round(mean(value),3)), by = species][order(species)])

sp_labs <- c("1 pioneer","2 early-mid","3 mid-late","4 climax")

# --- headline: single basal-area panel --------------------------------------
p_ba <- ggplot(ba, aes(year, value, colour = factor(species))) +
  geom_line(linewidth = 0.9) +
  scale_colour_viridis_d(name = "Species", labels = sp_labs) +
  labs(x = "Year", y = expression(Basal~area~(m^2/ha)),
       title = "Primary succession: basal area by species") +
  theme_minimal()
ggsave("dev/scratch_succession.png", p_ba, width = 7, height = 4.5, dpi = 110)

# --- classic multi-panel: all status variables ------------------------------
panel <- copy(p_dat)
# show regeneration as the mean per ha (r_mean_ha); drop per-patch counts (reg)
panel[, variable2 := factor(
  variable,
  levels = c("dbh","ba","trees","AL","growth","mort","r_mean_ha"),
  labels = c("avg. DBH [cm]","Basal Area [m²/ha]","Trees [N/ha]",
             "Available Light [%]","Growth [cm/yr]","Mortality [%]",
             "Regeneration [N/ha]"))]
p_multi <- ggplot(panel[!is.na(variable2)],
                  aes(year, value, colour = factor(species))) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~variable2, scales = "free_y", ncol = 2, strip.position = "left") +
  coord_cartesian(ylim = c(0, NA)) +
  scale_colour_viridis_d(name = "Species", labels = sp_labs) +
  labs(x = "Year", y = NULL, title = "Succession across all state variables") +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90),
        axis.title.y = element_blank())
ggsave("dev/scratch_succession_multi.png", p_multi, width = 8, height = 9, dpi = 110)

cat("\nsaved dev/scratch_succession.png + dev/scratch_succession_multi.png\n")
cat("\n=== variables present in sim$long$site ===\n")
print(sort(unique(p_dat$variable)))
