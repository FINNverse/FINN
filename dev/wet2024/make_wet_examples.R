# =============================================================================
# One worked example per WET2024 type, for the catalogue.
#
# Same illustrative stand for all types (a single, fast-growing species, hand-set
# — NOT fitted); only the MANAGEMENT PARAMETERS differ between panels, so each
# figure isolates what that WET's treatment does. Species-composition effects
# (Mischungsregulierung, different Leitbaumarten) are not shown by a single-species
# stand; they are documented in the catalogue table.
# =============================================================================
suppressMessages(devtools::load_all("../..", quiet = TRUE))   # run from dev/wet2024/
suppressMessages({library(data.table); library(ggplot2)})

wet <- fread("wet2024_parameters.csv")
H <- 150L
FINN.seed(1)
env <- data.table(expand.grid(siteID = 1:2, year = 1:H))[, climate := 0]

# common illustrative stand: reaches ~70 cm, ~44 m2/ha unmanaged
mk <- function() finn(N_species = 1L,
  competition_process  = createProcess(~0, FINN::competition,  initSpecies = matrix(c(.25, .09), 1, 2)),
  growth_process       = createProcess(~1, FINN::growth,       initSpecies = matrix(c(.05, .035), 1, 2), initEnv = list(matrix(-1.05, 1, 1))),
  regeneration_process = createProcess(~1, FINN::regeneration, initSpecies = matrix(.05, 1, 1),          initEnv = list(matrix(2.0, 1, 1))),  # active ingrowth (continuous-cover)
  mortality_process    = createProcess(~1, FINN::mortality,    initSpecies = matrix(c(0,0,0), 1, 3),      initEnv = list(matrix(-2.6, 1, 1))))
m <- mk()

ha <- function(x) { d <- setDT(x); d[variable == "ba", .(ba = sum(value)/.1/uniqueN(siteID)), by = year] }
base <- ha(predict(m, init_cohort = NULL, env = env, patches = 30, patch_size = .1, device = "cpu")$long$site)

dir.create("examples", showWarnings = FALSE)
summ <- list()
for (i in seq_len(nrow(wet))) {
  w <- wet[i]
  managed_ok <- is.finite(w$turnus_a)                 # Moorwald (m): no regular management
  if (managed_ok) {
    # Tending begins after canopy build-up (age 35), as in practice, so the stand
    # accumulates first and thinning then maintains it (rather than being knocked
    # down at establishment). Parameters are the catalogue values.
    yrs <- seq(35, H - 15, by = as.integer(w$turnus_a))
    par <- management_params(
      target_diameter         = if (is.finite(w$zielstaerke_cm)) w$zielstaerke_cm else Inf,
      target_harvest_fraction = w$nutzungsstaerke,
      target_diameter_spread  = if (is.finite(w$streubreite_cm)) w$streubreite_cm else 5,
      thinning_intensity      = w$durchf_staerke,
      thinning_size_bias      = w$durchf_art,
      risk_level              = if (is.finite(w$risikostufe)) w$risikostufe else 0)
    sched <- data.frame(year = yrs)
    sched$prescription <- I(replicate(length(yrs), par, simplify = FALSE))
    run <- simulate_managed(m, env = env, patches = 30, schedule = sched, patch_size = .1, device = "cpu")
    mba <- ha(run$trajectory)
    resid <- round(mean(tail(mba$ba, 6)), 1)
    dgaus <- round(mean(run$harvest$Dg_aus, na.rm = TRUE), 0)
    entries <- length(yrs)
  } else {
    mba <- base; resid <- round(mean(tail(base$ba, 6)), 1); dgaus <- NA; entries <- 0
  }
  summ[[i]] <- data.table(code = w$code, name = w$name, residual_ba = resid,
                          mean_Dg_aus = dgaus, entries = entries)

  df <- rbind(data.frame(year = base$year, ba = base$ba, regime = "unmanaged"),
              data.frame(year = mba$year,  ba = mba$ba,  regime = "under this WET"))
  p <- ggplot(df, aes(year, ba, colour = regime)) +
    { if (managed_ok) geom_vline(xintercept = yrs, colour = "grey90", linewidth = .25) } +
    geom_line(linewidth = .8) +
    scale_colour_manual(values = c("unmanaged" = "#9bb0a4", "under this WET" = "#b25a1c"), name = NULL) +
    coord_cartesian(ylim = c(0, max(base$ba) * 1.02)) +
    labs(x = "stand age (yr)", y = expression("BA (m"^2*" ha"^-1*")"),
         subtitle = sprintf("%s (%s)", w$name, w$code)) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", legend.margin = margin(t = -6),
          plot.subtitle = element_text(face = "bold", size = 9))
  ggsave(sprintf("examples/%s.png", w$code), p, width = 4.4, height = 2.5, dpi = 130)
}
summ <- rbindlist(summ)
fwrite(summ, "examples/summary.csv")
cat("=== WET example summary (residual BA, mean Dg_aus of removals) ===\n"); print(summ)
cat("\nwrote examples/<code>.png for", nrow(wet), "WETs\n")
