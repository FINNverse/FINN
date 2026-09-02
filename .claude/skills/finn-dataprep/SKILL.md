---
name: finn-dataprep
description: How to turn raw forest-inventory data into the tables FINN needs (obs_dt, env_dt, init cohorts) using makeObsData -> resolveSiteIDs -> makeInitCohorts. Use whenever preparing, loading, aligning, or debugging FINN input data, or when column/unit questions about obs_dt / env_dt / init_trees / dist_dt come up.
---

# Preparing data for FINN

FINN calibrates against clearly defined tables. Build them from two raw
ingredients — a **tree list** (one row per tree × inventory) and a **site
environment table** — with the pipeline
`makeObsData()` -> `resolveSiteIDs()` -> `makeInitCohorts()` (in `R/createInputs.R`).
The `C-Data_preparation` vignette is the worked reference.

## The tables FINN consumes

- **`obs_dt`** (one row per site × year × species): `siteID, year, ba, dbh, trees,
  growth, mort, n_at_risk, n_died, reg, species, species_name`. Year starts at 1;
  the initial state (year 0) is **not** in `obs_dt`. NAs allowed in any response.
  Units: `ba` m², `dbh` cm, `trees` count, `growth` relative growth rate, `reg`
  trees ha⁻¹. Mortality is a **pair of counts** — `n_at_risk` (alive at interval
  start) and `n_died` — which is what the default `mortality = "binomial"` loss
  consumes; `mort = n_died/n_at_risk` is derived for convenience.
- **`env_dt`**: `siteID, year, <env vars…>` in **natural units**; must cover every
  site × year in `obs_dt`. **Do not pre-scale** — FINN z-standardises internally
  (`env_autoscale`, default) and reuses the constants at prediction time.
- **`init_trees`** (→ `makeInitCohorts`): `siteID, patchID, year, species, dbh,
  trees, living` (the year-0 slice).
- **`dist_dt`** (optional): `siteID, year, intensity` (fraction of patches
  disturbed that year).

## Raw tree list -> obs_dt

A raw tree record needs `siteName, patchName, treeName, year` (calendar year),
`species_name`, `dbh` (cm), and a `status` of `new`/`alive`/`dead`.
`makeObsData()` reconstructs each tree's previous state itself by sorting on
`treeName`/`year` — there is deliberately **no tree-level `mort` flag**; mortality
is derived as a closed cohort (alive at start, counted dead at end), which avoids
mis-booking deaths when a tree's species id changes between visits.

```r
obs_list <- makeObsData(
  tree_dt           = tree_dt[complete == TRUE],  # only fully re-measured plots
  plotsize          = 0.06,        # ha per patch (scales recruitment to per-ha)
  aggregate_by_site = FALSE,       # keep patch-level resolution
  fix_period_length = 10,          # drop sites whose remeasure gap != 10 yr
  dbh_growth_thresh = c(-10, 50),  # drop implausible dbh change
  NspeciesQuantile  = 0.98,        # keep species covering 98% of stems; rest -> "other"
  minNyears         = 2)
obs_dt <- obs_list$obs_dt
```

Use `Nspecies` (hard cap) or `NspeciesQuantile` (smallest set covering a stem
fraction) to control species count; the rest become `"other"`, still tracked as a
competitor.

## Align IDs and build cohorts

Raw tables use names; FINN works with integer indices aligned across tables. Pass
the **processed** tree list `obs_list$tree_dt` (same species lumping as obs_dt),
not the raw one:

```r
inputs <- resolveSiteIDs(tree_dt = obs_list$tree_dt, obs_dt = obs_dt,
                         env_dt = env_dt, createInitCohorts = TRUE)
init_cohorts <- inputs$initCohorts     # a CohortMat, dims c(sites, patches, cohorts)
```

To build cohorts directly from a year-0 slice: `makeInitCohorts(init_trees,
Nspecies = uniqueN(obs_dt$species))` (expects `siteID, patchID, species, dbh,
treeName`, optional `trees`).

## Checklist before fitting

- Period length consistent (`fix_period_length`); FINN assumes a fixed remeasure interval.
- `env_dt` covers **every** `siteID × year` in `obs_dt` (missing rows break the fit).
- Units: `dbh` cm, `ba` m², `reg` trees ha⁻¹; env in natural units (never pre-scaled).
- Species IDs contiguous `1..N`, aligned across tables — always go through
  `resolveSiteIDs()`, never hand-code them.
- Acquisition steps that need external downloads (rFIA, WorldClim rasters) are
  shown `eval = FALSE` in the vignette — do not run them in package code/tests.
