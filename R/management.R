# =============================================================================
# FINN | R/management.R
#
# Prescriptive forest management for FINN simulations, as ONE continuous
# parametric operator.
#
# -----------------------------------------------------------------------------
# Why this module exists
# -----------------------------------------------------------------------------
# FINN's `disturbance` argument is a per-patch Bernoulli event that kills every
# tree in the patch -- an all-or-nothing, stand-replacing clearfell. It cannot
# express the size- and species-selective partial removals that make up most
# silviculture (thinning from below/above, Zielstärkennutzung, continuous-cover
# selection, species conversion). This module supplies that, but NOT as a set of
# discrete named regimes: every strategy -- including all the WET2024 forest
# development types -- is a POINT in one common, continuous parameter set (the
# "management genome"). See dev/wet2024/ for the catalogue and
# dev/gy_growth_parameterisation.md for how it is calibrated to the NW-FVA yield
# tables and the FVA-BW Durchforstungshilfe.
#
# The operator turns a parameter set into a per-cohort, per-entry REMOVAL FRACTION
#   r(dbh, s) = clamp(
#       target_harvest_fraction * logistic((dbh - target_diameter_s)/spread)  # Zielstärkennutzung
#     + thinning_intensity * size_weight(dbh; thinning_size_bias) * species_factor(s)  # Durchforstung
#     , 0, 1)
# and applies it to the cohort snapshot. `planting_rate` acts on regeneration, not
# removal. `risk_level` is a single dial that lowers the target diameter.
#
# -----------------------------------------------------------------------------
# Interpretation caveat -- state this whenever managed runs are reported
# -----------------------------------------------------------------------------
# A FINN model fitted to unmanaged (or observationally managed) demography has not
# "seen" the altered light regime, planting or regeneration that follow a real
# thinning. A prescriptive removal answers "given this fitted demography, what
# happens if we take out these trees on this schedule" -- a legitimate scenario,
# not a claim that the model has learned managed-stand behaviour.
# =============================================================================

# Cohort snapshot contract: a data.frame with at least `siteID, patchID, species,
# cohortID, dbh, trees`, matching simulateForest(..., return_cohorts = TRUE). dbh
# in cm, trees in stems ha-1. Per-species parameters are indexed by `species`.

# Risk dial: at risk_level = 1 the target diameter is lowered by this fraction
# (WET2024 lowers Zielstärke under higher risk, e.g. Fichte 50 -> 45 cm ~ 10%).
.RISK_ZIELSTAERKE_DROP <- 0.15


#' Basal area of a cohort snapshot (m2 ha-1)
#' @keywords internal
#' @noRd
cohort_basal_area <- function(cohorts, by = NULL) {
  ba <- pi * (cohorts$dbh / 200)^2 * cohorts$trees      # dbh cm -> m radius
  if (is.null(by)) return(sum(ba, na.rm = TRUE))
  stats::aggregate(list(ba = ba), by = cohorts[by], FUN = sum, na.rm = TRUE)
}

# Quadratic mean diameter (Dg, cm) of a snapshot.
.cohort_qmd <- function(cohorts) {
  w <- cohorts$trees
  if (sum(w, na.rm = TRUE) <= 0) return(NA_real_)
  sqrt(sum(w * cohorts$dbh^2, na.rm = TRUE) / sum(w, na.rm = TRUE))
}

# Resolve a scalar-or-per-species parameter to one value per cohort.
.per_species <- function(x, species) {
  if (length(x) == 1L) return(rep(x, length(species)))
  x[species]
}


#' Management parameters (the common "genome")
#'
#' Builds one point in FINN's continuous management parameter space. Every
#' strategy -- thinning, target-diameter harvest, species conversion, and each of
#' the WET2024 forest development types -- is expressed as a value of these
#' parameters, never a discrete case. Names carry the WET vocabulary.
#'
#' @param target_diameter `Zielstärke` (cm): the dbh at which a tree is harvested.
#'   Scalar or a per-species vector. `Inf` (default) means no target-diameter
#'   harvest.
#' @param target_harvest_fraction `Nutzungsstärke` (0-1): share of at/above-target
#'   stems removed per entry.
#' @param target_diameter_spread `Zielstärken-Streubreite` (cm, > 0): how sharp
#'   (small) vs. gradual (large) the target-diameter threshold is.
#' @param thinning_intensity `Durchforstungsstärke` (0-1): mean fraction removed
#'   per entry among the tended (sub-target) trees.
#' @param thinning_size_bias `Durchforstungsart`: which sizes the thinning takes.
#'   `< 0` from below (Niederdurchforstung), `> 0` from above (Hochdurchforstung),
#'   `0` neutral. Either a scalar OR a `function(stand)` returning a scalar, where
#'   `stand` is a list with `Dg` (quadratic mean dbh) and `target` (per-cohort
#'   target diameter) -- use a function to let the bias swing from below to above
#'   as the stand matures (the NW-FVA graduated regime; see [graduated_bias()]).
#' @param species_removal_pref `Mischungsregulierung / Baumartenwahl`: scalar or
#'   per-species. `0` neutral, `> 0` removed preferentially (convert away, e.g.
#'   spruce), `< 0` retained/favoured. Enters as a factor `1 + pref`.
#' @param planting_rate `Verjüngung / Vorbau` (0-1): advance-planting rate of
#'   target species. Acts on regeneration, not removal (used by the driver).
#' @param entry_interval_years `Eingriffsturnus`: years between interventions.
#' @param risk_level `Risikostufe` (0-1): a single dial that lowers
#'   `target_diameter` (by up to 15% at `risk_level = 1`).
#' @return An object of class `management_params`.
#' @seealso [apply_management()], [graduated_bias()], [thin_from_below()],
#'   [target_diameter_harvest()]
#' @export
management_params <- function(target_diameter = Inf,
                              target_harvest_fraction = 0,
                              target_diameter_spread = 5,
                              thinning_intensity = 0,
                              thinning_size_bias = 0,
                              species_removal_pref = 0,
                              planting_rate = 0,
                              entry_interval_years = 5,
                              risk_level = 0) {
  stopifnot(all(target_diameter > 0), target_diameter_spread > 0,
            target_harvest_fraction >= 0, target_harvest_fraction <= 1,
            thinning_intensity >= 0, thinning_intensity <= 1,
            risk_level >= 0, risk_level <= 1, entry_interval_years > 0,
            planting_rate >= 0, planting_rate <= 1)
  if (!is.function(thinning_size_bias) && length(thinning_size_bias) != 1L)
    stop("thinning_size_bias must be a scalar or a function(stand).", call. = FALSE)
  structure(list(target_diameter = target_diameter,
                 target_harvest_fraction = target_harvest_fraction,
                 target_diameter_spread = target_diameter_spread,
                 thinning_intensity = thinning_intensity,
                 thinning_size_bias = thinning_size_bias,
                 species_removal_pref = species_removal_pref,
                 planting_rate = planting_rate,
                 entry_interval_years = entry_interval_years,
                 risk_level = risk_level),
            class = "management_params")
}

#' @export
print.management_params <- function(x, ...) {
  bias <- if (is.function(x$thinning_size_bias)) "<function of stand>" else x$thinning_size_bias
  cat("<management_params>  (one point in FINN's management genome)\n")
  cat(sprintf("  Zielstärke (target_diameter)      : %s cm\n", paste(x$target_diameter, collapse = ", ")))
  cat(sprintf("  Nutzungsstärke (target_harvest)   : %s\n", x$target_harvest_fraction))
  cat(sprintf("  Zielstärken-Streubreite (spread)  : %s cm\n", x$target_diameter_spread))
  cat(sprintf("  Durchforstungsstärke (intensity)  : %s\n", x$thinning_intensity))
  cat(sprintf("  Durchforstungsart (size_bias)     : %s\n", bias))
  cat(sprintf("  Mischungsregulierung (species)    : %s\n", paste(x$species_removal_pref, collapse = ", ")))
  cat(sprintf("  Verjüngung (planting_rate)        : %s\n", x$planting_rate))
  cat(sprintf("  Eingriffsturnus (interval, yr)    : %s\n", x$entry_interval_years))
  cat(sprintf("  Risikostufe (risk_level)          : %s\n", x$risk_level))
  invisible(x)
}


#' Per-cohort removal fraction for a management parameter set
#'
#' Evaluates the operator's removal-fraction field on a cohort snapshot: the
#' Zielstärkennutzung (logistic ramp above the target diameter) plus the
#' Durchforstung (size- and species-weighted tending), clamped to \[0, 1\].
#'
#' @param cohorts a cohort snapshot (see the module contract).
#' @param params a [management_params] object.
#' @return numeric vector, one removal fraction per cohort row.
#' @export
management_removal <- function(cohorts, params) {
  stopifnot(inherits(params, "management_params"))
  dbh <- cohorts$dbh; trees <- cohorts$trees; sp <- cohorts$species
  n <- length(dbh)
  target <- .per_species(params$target_diameter, sp)
  pref   <- .per_species(params$species_removal_pref, sp)
  target <- target * (1 - params$risk_level * .RISK_ZIELSTAERKE_DROP)   # Risikostufe

  Dg <- .cohort_qmd(cohorts)

  # Durchforstung: size-weighted, normalised so the trees-weighted mean weight is
  # 1 (=> average tended removal ~ thinning_intensity), times the species factor.
  bias <- if (is.function(params$thinning_size_bias))
    params$thinning_size_bias(list(Dg = Dg, target = target, dbh = dbh)) else params$thinning_size_bias
  sw <- if (is.finite(Dg) && Dg > 0) (dbh / Dg)^bias else rep(1, n)
  wbar <- sum(trees * sw, na.rm = TRUE) / sum(trees, na.rm = TRUE)
  if (is.finite(wbar) && wbar > 0) sw <- sw / wbar
  r_thin <- params$thinning_intensity * sw * (1 + pref)

  # Zielstärkennutzung: harvest ramps up as dbh crosses the target diameter.
  r_harv <- params$target_harvest_fraction *
    stats::plogis((dbh - target) / params$target_diameter_spread)

  pmin(pmax(r_thin + r_harv, 0), 1)
}


#' Apply a management prescription to a cohort snapshot
#'
#' Runs a prescription and returns the thinned snapshot plus the removal it
#' implied (stems, basal area, and the mean diameter of the removed trees, so a
#' run can be reported and validated against a yield table's `Dg_aus`). A
#' prescription is either a [management_params] object (the parametric operator)
#' or, for full control, a `function(cohorts)` returning a modified snapshot.
#'
#' @param cohorts a cohort snapshot (see the module contract).
#' @param prescription a [management_params] object or a `function(cohorts)`.
#' @return a list with `cohorts` (after), `removed_trees`, `removed_ba`,
#'   `Dg_removed` (mean dbh of removed trees, cm) and `Dg_aus_rel`
#'   (`Dg_removed / Dg_before`, < 1 from below, > 1 from above).
#' @export
apply_management <- function(cohorts, prescription) {
  before_trees <- sum(cohorts$trees, na.rm = TRUE)
  before_ba <- cohort_basal_area(cohorts)
  Dg_before <- .cohort_qmd(cohorts)

  if (inherits(prescription, "management_params")) {
    r <- management_removal(cohorts, prescription)
    removed_n_i <- cohorts$trees * r
    after <- cohorts
    after$trees <- cohorts$trees - removed_n_i
  } else if (is.function(prescription)) {
    after <- prescription(cohorts)
    if (!identical(dim(after), dim(cohorts)))
      stop("a prescription function must not add or drop cohort rows.", call. = FALSE)
    removed_n_i <- cohorts$trees - after$trees
  } else {
    stop("prescription must be a management_params object or a function(cohorts).",
         call. = FALSE)
  }

  rem_tot <- sum(removed_n_i, na.rm = TRUE)
  Dg_removed <- if (rem_tot > 0)
    sqrt(sum(removed_n_i * cohorts$dbh^2, na.rm = TRUE) / rem_tot) else NA_real_
  list(cohorts = after,
       removed_trees = rem_tot,
       removed_ba = before_ba - cohort_basal_area(after),
       Dg_removed = Dg_removed,
       Dg_aus_rel = Dg_removed / Dg_before)
}


#' A size-dependent thinning bias (graduated Durchforstungsart)
#'
#' Returns a `function(stand)` for `management_params(thinning_size_bias = ...)`
#' that swings the thinning from below to above as the stand matures -- the NW-FVA
#' "gestaffelte Durchforstung" / FVA-BW Durchforstungshilfe pattern, where young
#' stands are tended from below and, as the mean diameter approaches the target,
#' harvesting shifts to the strong (from above). The bias interpolates smoothly
#' from `from` (at `Dg = 0`) to `to` (at `Dg >= target`).
#'
#' @param from bias at the start of stand development (default -0.6, from below).
#' @param to bias near the target diameter (default 0.8, from above).
#' @param ref reference diameter for the swing; defaults to the target diameter.
#' @return a `function(stand)` returning a scalar bias.
#' @export
graduated_bias <- function(from = -0.6, to = 0.8, ref = NULL) {
  function(stand) {
    r <- if (is.null(ref)) stand$target else ref
    r <- suppressWarnings(min(r[is.finite(r)], na.rm = TRUE))
    if (!is.finite(r) || r <= 0) return(from)
    x <- max(0, min(1, stand$Dg / r))       # 0 (young) .. 1 (at target)
    from + (to - from) * (3 * x^2 - 2 * x^3) # smoothstep
  }
}


# -----------------------------------------------------------------------------
# Presets: each returns a management_params (a point in the genome), not a case.
# -----------------------------------------------------------------------------

#' Thinning from below (Niederdurchforstung) -- a parameter preset
#' @param intensity `Durchforstungsstärke` (0-1) removed per entry.
#' @param ... further [management_params] overrides (e.g. `entry_interval_years`).
#' @return a [management_params] object.
#' @export
thin_from_below <- function(intensity = 0.2, ...)
  management_params(thinning_intensity = intensity, thinning_size_bias = -1, ...)

#' Thinning from above (Hochdurchforstung) -- a parameter preset
#' @inheritParams thin_from_below
#' @return a [management_params] object.
#' @export
thin_from_above <- function(intensity = 0.2, ...)
  management_params(thinning_intensity = intensity, thinning_size_bias = 1, ...)

#' Target-diameter harvest (Zielstärkennutzung) -- a parameter preset
#' @param target_diameter `Zielstärke` (cm), scalar or per-species.
#' @param fraction `Nutzungsstärke` (0-1) of at/above-target stems taken per entry.
#' @param spread `Zielstärken-Streubreite` (cm).
#' @param ... further [management_params] overrides.
#' @return a [management_params] object.
#' @export
target_diameter_harvest <- function(target_diameter, fraction = 1, spread = 5, ...)
  management_params(target_diameter = target_diameter,
                    target_harvest_fraction = fraction,
                    target_diameter_spread = spread, ...)

#' Species-selective removal (conversion / sanitation) -- a parameter preset
#'
#' Removes a fraction of the named species and protects the rest, via
#' `Mischungsregulierung`. Requires the total number of species so the per-species
#' preference vector can be built.
#' @param species integer index/indices to remove.
#' @param n_species total number of species in the model.
#' @param fraction fraction of the targeted species removed per entry (0-1).
#' @param ... further [management_params] overrides.
#' @return a [management_params] object.
#' @export
species_removal <- function(species, n_species, fraction = 1, ...) {
  pref <- rep(-1, n_species)     # protect others (factor 1 + (-1) = 0)
  pref[species] <- 0             # neutral on the targeted species (factor 1)
  management_params(thinning_intensity = fraction, species_removal_pref = pref,
                    thinning_size_bias = 0, ...)
}

#' Prescribed clearfell -- a parameter preset
#'
#' Removes (nearly) all stems in one entry. A deterministic, scheduled counterpart
#' to FINN's stochastic `disturbance` clearfell.
#' @param ... further [management_params] overrides.
#' @return a [management_params] object.
#' @export
clearfell <- function(...)
  management_params(thinning_intensity = 1, thinning_size_bias = 0, ...)


# -----------------------------------------------------------------------------
# Stepwise driver -- contract fixed, implementation still pending
# -----------------------------------------------------------------------------
#' Simulate a fitted FINN model under a management schedule
#'
#' Threads a management schedule through a fitted model by simulating in segments:
#' run to the next scheduled entry with `simulateForest(return_cohorts = TRUE)`,
#' apply the prescription to the final cohort snapshot with [apply_management()],
#' re-initialise from the thinned state, and continue. The post-thinning light
#' release and growth response are then produced by the fitted demography itself.
#'
#' @param model a fitted `finn_class` model.
#' @param env the environment table (as for `simulateForest`).
#' @param patches number of patches per site.
#' @param schedule a data.frame with `year` and a list-column `prescription`
#'   ([management_params] objects or functions).
#' @param ... passed to `simulateForest()`.
#' @return a list with the stitched `$trajectory` and a `$harvest` log.
#' @export
simulate_managed <- function(model, env, patches, schedule, ...) {
  stop("simulate_managed(): segmented driver not yet implemented on this branch. ",
       "The parametric operator (management_params/management_removal/apply_management) ",
       "is ready; this driver is the next build (segment the horizon at schedule$year, ",
       "re-init from the thinned cohort snapshot between segments).", call. = FALSE)
}
