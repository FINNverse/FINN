# =============================================================================
# One worked example per WET2024 type, for the catalogue.
#
# Same illustrative stand for all types (a single, productive species, hand-set
# — NOT fitted); only the MANAGEMENT PARAMETERS differ between panels, so each
# figure isolates what that WET's treatment does. Species-composition effects
# (Mischungsregulierung, different Leitbaumarten) are not shown by a single-species
# stand; they are documented in the catalogue table.
#
# Y axis = standing volume (Vorrat, Vfm/ha) = G * H(Dg) * form factor, the German-
# forestry stocking measure. Realistic reference: German forests average ~335
# Vfm/ha (BWI4 2022); managed conifer stands hold ~350-550 Vfm/ha.
# =============================================================================
suppressMessages(devtools::load_all("../..", quiet = TRUE))   # run from dev/wet2024/
suppressMessages({library(data.table); library(ggplot2)})

wet <- fread("wet2024_parameters.csv")
H <- 150L
FINN.seed(1)
env <- data.table(expand.grid(siteID = 1:2, year = 1:H))[, climate := 0]

# Productive illustrative stand that can sustain tending at high stocking.
mk <- function() finn(N_species = 1L,
  competition_process  = createProcess(~0, FINN::competition,  initSpecies = matrix(c(.25, .09), 1, 2)),
  growth_process       = createProcess(~1, FINN::growth,       initSpecies = matrix(c(.05, .030), 1, 2), initEnv = list(matrix(-0.95, 1, 1))),
  regeneration_process = createProcess(~1, FINN::regeneration, initSpecies = matrix(.05, 1, 1),          initEnv = list(matrix(2.0, 1, 1))),
  mortality_process    = createProcess(~1, FINN::mortality,    initSpecies = matrix(c(0,0,0), 1, 3),      initEnv = list(matrix(-2.4, 1, 1))))
m <- mk()

# stand-level series: per-ha basal area, mean diameter, height, standing volume.
FORM <- 0.47                                   # Derbholz-Formzahl (conifer ~0.45-0.5)
Hcurve <- function(dbh) 1.3 + 38 * (1 - exp(-0.033 * dbh))   # height (m) from Dg (cm)
stand <- function(x) {
  d <- setDT(x)
  ba <- d[variable == "ba",  .(BA = sum(value)/.1/uniqueN(siteID)), by = year]
  dg <- d[variable == "dbh", .(Dg = mean(value)), by = year]
  m  <- merge(ba, dg, by = "year")
  m[, Vfm := BA * Hcurve(Dg) * FORM][]
}
base <- stand(predict(m, init_cohort = NULL, env = env, patches = 30, patch_size = .1, device = "cpu")$long$site)

# Applied intensities are scaled down from the catalogue values: a "massige
# Durchforstung" removes ~15-20% of basal area per entry, and holds the stand at
# high stocking rather than knocking it down. THIN_SCALE / HARV_SCALE encode that
# the catalogue's first-pass intensities were too heavy for a 5-7 yr Turnus.
THIN_SCALE <- 0.85; HARV_SCALE <- 0.15
dir.create("examples", showWarnings = FALSE)
summ <- list()
for (i in seq_len(nrow(wet))) {
  w <- wet[i]
  managed_ok <- is.finite(w$turnus_a)                 # Moorwald (m): no regular management
  if (managed_ok) {
    yrs <- seq(35, H - 15, by = as.integer(w$turnus_a))   # tending after canopy build-up
    par <- management_params(
      target_diameter         = if (is.finite(w$zielstaerke_cm)) w$zielstaerke_cm else Inf,
      target_harvest_fraction = w$nutzungsstaerke * HARV_SCALE,
      target_diameter_spread  = if (is.finite(w$streubreite_cm)) w$streubreite_cm else 5,
      thinning_intensity      = w$durchf_staerke * THIN_SCALE,
      thinning_size_bias      = w$durchf_art,
      risk_level              = if (is.finite(w$risikostufe)) w$risikostufe else 0)
    sched <- data.frame(year = yrs)
    sched$prescription <- I(replicate(length(yrs), par, simplify = FALSE))
    run <- simulate_managed(m, env = env, patches = 30, schedule = sched, patch_size = .1, device = "cpu")
    ms <- stand(run$trajectory)
    resid <- round(mean(tail(ms$Vfm, 6)))
    dgaus <- round(mean(run$harvest$Dg_aus, na.rm = TRUE), 0)
    entries <- length(yrs)
  } else {
    ms <- base; resid <- round(mean(tail(base$Vfm, 6))); dgaus <- NA; entries <- 0
  }
  summ[[i]] <- data.table(code = w$code, name = w$name, residual_vfm = resid,
                          mean_Dg_aus = dgaus, entries = entries)

  df <- rbind(data.frame(year = base$year, Vfm = base$Vfm, regime = "unmanaged"),
              data.frame(year = ms$year,   Vfm = ms$Vfm,   regime = "under this WET"))
  p <- ggplot(df, aes(year, Vfm, colour = regime)) +
    { if (managed_ok) geom_vline(xintercept = yrs, colour = "grey90", linewidth = .25) } +
    geom_line(linewidth = .8) +
    scale_colour_manual(values = c("unmanaged" = "#9bb0a4", "under this WET" = "#b25a1c"), name = NULL) +
    coord_cartesian(ylim = c(0, max(base$Vfm) * 1.02)) +
    labs(x = "stand age (yr)", y = expression("Vorrat (Vfm ha"^-1*")"),
         subtitle = sprintf("%s (%s)", w$name, w$code)) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", legend.margin = margin(t = -6),
          plot.subtitle = element_text(face = "bold", size = 9))
  ggsave(sprintf("examples/%s.png", w$code), p, width = 4.4, height = 2.5, dpi = 130)
}
summ <- rbindlist(summ)
fwrite(summ, "examples/summary.csv")
cat(sprintf("unmanaged final Vorrat: %.0f Vfm/ha\n", mean(tail(base$Vfm, 6))))
cat("=== WET example summary (residual Vorrat Vfm/ha, mean Dg_aus) ===\n"); print(summ)
