# =============================================================================
# FINN | R/management.R
#
# Prescriptive forest management for FINN simulations.
#
# -----------------------------------------------------------------------------
# Why this module exists
# -----------------------------------------------------------------------------
# FINN already accepts a `disturbance` argument in fit()/predict()/simulateForest().
# Internally that is a per-patch, per-year Bernoulli event that, when it fires,
# multiplies *every* tree in the patch by zero -- an all-or-nothing, STAND-REPLACING
# kill. It can represent clearfelling on a fixed rotation, or catastrophic natural
# disturbance (windthrow, bark beetle), but it *cannot* represent the size-selective
# partial removals that make up most silviculture:
#
#   * thinning from below / from above (Nieder- / Hochdurchforstung)
#   * target-diameter harvest           (Zielstärkennutzung)
#   * continuous-cover / selection cuts (Plenterung, Dauerwald)
#   * species conversion                (Umbau, e.g. spruce -> mixed)
#
# These act on the cohort state -- which trees, of which size and species, are
# removed -- so they need a size-selective operator, not a scalar probability.
# This module supplies that operator layer as a *prescriptive* (rule-based, not
# fitted) intervention applied to the cohort snapshot between simulation steps.
#
# -----------------------------------------------------------------------------
# Interpretation caveat -- state this whenever managed runs are reported
# -----------------------------------------------------------------------------
# A FINN model fitted to unmanaged (or observationally managed) demography has not
# "seen" the altered light regime, planting, or regeneration response that follow a
# real thinning. A prescriptive removal answers "given this fitted demography, what
# happens to the stand if we take out these trees on this schedule" -- a legitimate
# scenario, but not a claim that the model has learned managed-stand behaviour.
# See docs/harvest_process_design.md in the FINN-bwi project for the alternative
# (harvest as a *fitted* process / env covariate) and when to prefer it.
#
# Design status: the prescription constructors and `apply_management()` below are
# complete and unit-testable (pure functions on a cohort data.frame). The stepwise
# driver `simulate_managed()` that threads them through a fitted model is the piece
# under construction on this branch -- its contract is fixed here.
# =============================================================================


# -----------------------------------------------------------------------------
# Cohort snapshot contract
# -----------------------------------------------------------------------------
# A "cohort snapshot" is a data.frame with (at least) these columns, matching the
# cohort table returned by `simulateForest(..., return_cohorts = TRUE)`:
#
#   siteID   integer   site identifier
#   patchID  integer   patch (replicate stand) within the site
#   species  integer   species index (1..n_species)
#   cohortID integer   cohort identifier within the patch
#   dbh      numeric   cohort mean diameter at breast height (cm)
#   trees    numeric   stems per hectare represented by the cohort
#
# A management prescription is a function
#     f(cohorts) -> cohorts
# that returns the snapshot with `trees` reduced by whatever it removed (a cohort
# thinned to zero stems is kept as a zero row, not dropped, so downstream indexing
# is stable). It must not change dbh/species/identity columns.


#' Basal area of a cohort snapshot
#'
#' Cross-sectional area at breast height, m2 ha-1, summed over the snapshot (or
#' over the grouping given by `by`). `dbh` is in cm, `trees` in stems ha-1.
#'
#' @param cohorts a cohort snapshot (see the module contract).
#' @param by optional character vector of grouping columns; `NULL` totals everything.
#' @return a numeric total, or a data.frame with one `ba` value per group.
#' @keywords internal
#' @noRd
cohort_basal_area <- function(cohorts, by = NULL) {
  ba_stem <- pi * (cohorts$dbh / 200)^2          # m2 per stem (dbh cm -> m radius)
  ba <- ba_stem * cohorts$trees                  # m2 ha-1 per cohort
  if (is.null(by)) return(sum(ba, na.rm = TRUE))
  stats::aggregate(list(ba = ba), by = cohorts[by], FUN = sum, na.rm = TRUE)
}


#' Thinning from below (Niederdurchforstung)
#'
#' Removes the smallest-diameter stems first until the residual basal area is at
#' most `residual_ba` (m2 ha-1), applied within each patch. The classic even-aged
#' tending cut: it releases the dominant crop trees by taking the suppressed ones.
#'
#' @param residual_ba target basal area to leave standing, m2 ha-1.
#' @param min_dbh only stems at or above this dbh (cm) are eligible; smaller
#'   regeneration is left. Default `0` (all eligible).
#' @return a management prescription (a function of a cohort snapshot).
#' @export
thin_from_below <- function(residual_ba, min_dbh = 0) {
  stopifnot(residual_ba >= 0, min_dbh >= 0)
  .thin_directional(residual_ba, min_dbh, from = "below")
}


#' Thinning from above (Hochdurchforstung)
#'
#' Removes the largest-diameter stems first until the residual basal area is at
#' most `residual_ba`. Favours the mid-storey; used to harvest value early or to
#' break up a closing canopy.
#'
#' @inheritParams thin_from_below
#' @return a management prescription.
#' @export
thin_from_above <- function(residual_ba, min_dbh = 0) {
  stopifnot(residual_ba >= 0, min_dbh >= 0)
  .thin_directional(residual_ba, min_dbh, from = "above")
}


# Shared engine for directional thinning to a residual basal area, per patch.
.thin_directional <- function(residual_ba, min_dbh, from = c("below", "above")) {
  from <- match.arg(from)
  function(cohorts) {
    parts <- split(cohorts, list(cohorts$siteID, cohorts$patchID), drop = TRUE)
    out <- lapply(parts, function(p) {
      elig <- p$dbh >= min_dbh & p$trees > 0
      if (!any(elig)) return(p)
      # cumulative BA walking from the removal end inwards
      ord <- order(p$dbh[elig], decreasing = (from == "above"))
      idx <- which(elig)[ord]
      ba_stem <- pi * (p$dbh[idx] / 200)^2
      ba_here <- ba_stem * p$trees[idx]
      total <- sum(ba_stem * p$trees)            # eligible BA only
      to_remove <- max(0, total - residual_ba)
      keep_frac <- .partial_removal(ba_here, to_remove)
      p$trees[idx] <- p$trees[idx] * keep_frac
      p
    })
    do.call(rbind, unname(out))[order(cohorts$siteID, cohorts$patchID), , drop = FALSE]
  }
}


# Given per-cohort removable BA (in removal order) and a BA amount to remove,
# return the per-cohort *surviving* fraction. The cohort that straddles the target
# is thinned proportionally; everything past it is untouched.
.partial_removal <- function(ba_here, to_remove) {
  keep <- rep(1, length(ba_here))
  remaining <- to_remove
  for (i in seq_along(ba_here)) {
    if (remaining <= 0) break
    if (ba_here[i] <= remaining) {
      keep[i] <- 0
      remaining <- remaining - ba_here[i]
    } else {
      keep[i] <- 1 - remaining / ba_here[i]
      remaining <- 0
    }
  }
  keep
}


#' Target-diameter harvest (Zielstärkennutzung)
#'
#' Removes a fraction `p` of the stems in every cohort whose diameter has reached
#' or exceeded a species-specific target. The backbone of continuous-cover
#' forestry: trees are harvested individually once they are big enough, and the
#' stand is never cleared.
#'
#' @param dbh_target target diameter (cm). Either a single number applied to all
#'   species, or a named/indexed numeric vector giving a target per species.
#' @param p fraction of the at-or-above-target stems removed per entry (0-1).
#'   Default `1` (take all mature stems this step).
#' @return a management prescription.
#' @export
target_diameter_harvest <- function(dbh_target, p = 1) {
  stopifnot(p >= 0, p <= 1, all(dbh_target > 0))
  function(cohorts) {
    tgt <- if (length(dbh_target) == 1L) rep(dbh_target, nrow(cohorts))
           else dbh_target[cohorts$species]
    mature <- cohorts$dbh >= tgt
    cohorts$trees[mature] <- cohorts$trees[mature] * (1 - p)
    cohorts
  }
}


#' Species-selective removal (conversion / sanitation)
#'
#' Removes a fraction `p` of the stems of the named species, optionally only above
#' a diameter. Represents species conversion (e.g. taking spruce out of a
#' vulnerable stand to favour broadleaves) or sanitation felling.
#'
#' @param species integer species index or vector of indices to target.
#' @param p fraction removed per entry (0-1). Default `1`.
#' @param min_dbh only stems at or above this dbh (cm) are taken. Default `0`.
#' @return a management prescription.
#' @export
species_removal <- function(species, p = 1, min_dbh = 0) {
  stopifnot(p >= 0, p <= 1, min_dbh >= 0)
  function(cohorts) {
    hit <- cohorts$species %in% species & cohorts$dbh >= min_dbh
    cohorts$trees[hit] <- cohorts$trees[hit] * (1 - p)
    cohorts
  }
}


#' Prescribed clearfell
#'
#' Removes every stem (optionally only of some species). A deterministic,
#' scheduled counterpart to FINN's stochastic `disturbance` clearfell -- use this
#' when the *timing* is a decision (a rotation-end cut on chosen years) rather than
#' a random hazard.
#'
#' @param species optional integer indices; `NULL` (default) fells all species.
#' @return a management prescription.
#' @export
clearfell <- function(species = NULL) {
  function(cohorts) {
    hit <- if (is.null(species)) rep(TRUE, nrow(cohorts))
           else cohorts$species %in% species
    cohorts$trees[hit] <- 0
    cohorts
  }
}


#' Apply a management prescription to a cohort snapshot
#'
#' Runs the prescription and returns both the thinned snapshot and the removal it
#' implied (stems and basal area, so a run can report its harvest). This is the
#' single point where a prescription touches state; the stepwise driver calls it
#' once per scheduled entry.
#'
#' @param cohorts a cohort snapshot (see the module contract).
#' @param prescription a function returned by one of the constructors above.
#' @return a list with `cohorts` (after), `removed_trees`, `removed_ba`.
#' @export
apply_management <- function(cohorts, prescription) {
  stopifnot(is.function(prescription))
  before_trees <- sum(cohorts$trees, na.rm = TRUE)
  before_ba <- cohort_basal_area(cohorts)
  after <- prescription(cohorts)
  if (!identical(dim(after), dim(cohorts)))
    stop("a prescription must not add or drop cohort rows", call. = FALSE)
  list(cohorts = after,
       removed_trees = before_trees - sum(after$trees, na.rm = TRUE),
       removed_ba = before_ba - cohort_basal_area(after))
}


# -----------------------------------------------------------------------------
# Stepwise driver -- contract fixed here, implementation in progress
# -----------------------------------------------------------------------------
#' Simulate a fitted FINN model under a management schedule
#'
#' Threads prescriptive management through a fitted model by simulating in
#' segments: run the model to the next scheduled entry with
#' `simulateForest(return_cohorts = TRUE)`, apply the prescription to the final
#' cohort snapshot with [apply_management()], re-initialise from the thinned state,
#' and continue. Because the removal enters the state the model then propagates,
#' the post-thinning light release and growth response are produced by the fitted
#' demography itself -- no core change and no refit.
#'
#' @param model a fitted `finn_class` model.
#' @param env the environment table to simulate over (as for `simulateForest`).
#' @param patches number of patches (replicate stands) per site.
#' @param schedule a data.frame with columns `year` and `prescription` (a
#'   list-column of prescription functions), giving what to do and when.
#' @param ... passed to `simulateForest()` (e.g. `device`).
#' @return a list with the stitched `$trajectory` (per-year stand state) and a
#'   `$harvest` log (removed stems/BA per entry).
#' @export
simulate_managed <- function(model, env, patches, schedule, ...) {
  stop("simulate_managed(): segmented driver not yet implemented on this branch. ",
       "The prescription operators and apply_management() are ready; this driver ",
       "is the next build (segment the horizon at schedule$year, re-init from the ",
       "thinned cohort snapshot between segments).", call. = FALSE)
}
