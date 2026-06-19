# FINN — vignette build (Claude Code project)

## Mission
Author three package vignettes for the FINN R package that let a new reader go
from install to a real result quickly. The three vignettes are independent but
ordered for a reader:

1. `vignettes/succession-demo.qmd` — **plausible succession from fictional
   species.** Define 3–5 species with contrasting demographic trade-offs and
   show FINN produces ecologically sensible succession (pioneers peak early,
   shade-tolerants dominate late). No data, no fitting — pure simulation.
2. `vignettes/data-preparation.qmd` — **how to prepare data for FINN.**
   Standalone reference: the exact tables FINN needs, their columns/units, and
   the `makeObsData()` → `resolveSiteIDs()` → `makeInitCohorts()` pipeline.
3. `vignettes/fit-to-fia.qmd` — **fit FINN to data (Oregon FIA subset).** Load a
   small shipped subset, build the FINN object, fit, inspect fit + parameters.

## Hard constraints (do not violate)
- These are **package vignettes** (in `vignettes/`, shipped, run during
  `R CMD check`). They must:
  - knit from a **clean install** with **no external downloads** and **no reads
    outside the package** (bundle every input as a small data file);
  - run **fast** (target < ~30 s each). Guard any heavy fit with
    `eval = FALSE` and load a pre-fitted model object instead.
- Engine is **Quarto**: YAML uses `%\VignetteEngine{quarto::html}`. Match the
  style of the existing `vignettes/Introduction-to-FINN.qmd`.
- Do **not** change exported function signatures or the package API. If a
  vignette needs a helper, put it inline in the chunk.
- Keep shipped data tiny. CRAN source tarball limit is 5 MB; aim for each
  bundled dataset well under 1 MB (subsample aggressively).
- Do not edit files under `playground/` or `FINNetAl-main/` — reference only.

## FINN API quick reference (verified against R/ source)
Model construction:
```r
m <- finn(
  N_species            = Nsp,
  recruits_dbh         = 1.0,                 # starting dbh of recruits
  competition_process  = createProcess(~0,        FINN::competition),
  growth_process       = createProcess(~1+env1,   FINN::growth,        initSpecies = parGrowth, initEnv = parGrowthEnv),
  regeneration_process = createProcess(~1+env1,   FINN::regeneration,  initSpecies = parReg,    initEnv = parRegEnv),
  mortality_process    = createProcess(~1+env1,   FINN::mortality,     initSpecies = parMort,   initEnv = parMortEnv)
)
```
- `createProcess(formula, func, initSpecies=, initEnv=, optimizeSpecies=, optimizeEnv=, hidden=, NN=, ...)`.
  - `formula` is the environment model, e.g. `~0` (none), `~1+env1`, `~.` (all env cols).
  - For **fitting**, set `optimizeSpecies=TRUE`/`optimizeEnv=TRUE`. For **simulation
    demos**, pass fixed `initSpecies`/`initEnv` and leave optimize flags off.
  - Replace a process with a neural net via `hidden=` / `NN=` (hybrid; not needed
    for these three vignettes but available).

Simulation (no fitting):
```r
sim <- m$simulate(init_cohort = NULL, env = env_dt, disturbance = dist_dt,
                  device = "cpu", patches = 100)
```
- `sim$long$site` is a long data.table with columns `year, species, variable, value`.
  Variables: `dbh, ba, trees, AL, growth, mort, reg, r_mean_ha`.
- `ba`/`trees` are per patch_size; divide by `patch_size` (default 0.1) for per-ha.
- Alternatively `pred2DF(sim)` returns a tidy data.frame.

Fitting:
```r
m$fit(env = env_dt, data = obs_dt, init_cohort = init_cohorts,
      device = "cpu", epochs = 500L, batchsize = 500L,
      patches = 4, patch_size = 0.06,
      weights = c(0.1, 10, 1, 10, 1, 1), lr = 0.01)
# or the exported wrapper: fit(m, data=obs_dt, env=env_dt, init_cohort=..., ...)
```
- `loss` default: `c(dbh="mse", ba="mse", trees="poisson", growth="mse",
  mortality="mse", regeneration="nbinom")`.

Data builders (see `R/createInputs.R`):
- `makeObsData(tree_dt, plotsize, aggregate_by_site=, fix_period_length=,
  dbh_growth_thresh=, Npatches=, Nspecies= | NspeciesQuantile=, minNyears=)`
  → list with `$obs_dt` and `$tree_dt`. Lumps rare species into `"other"`.
- `resolveSiteIDs(tree_dt, obs_dt, env_dt, createInitCohorts=)` → aligns
  site/patch/species IDs across the three tables.
- `makeInitCohorts(init_trees, Nspecies=)` → `CohortMat` for `init_cohort=`.

## Data formats (column contracts)
- **obs_dt** (one row per site×year×species):
  `siteID, year, ba, dbh, trees, growth, mort, reg, species, species_name`.
  Year starts at 1; the init state (year 0) is NOT in obs_dt. NAs allowed in any
  response. Units: ba m²/patch, dbh cm, trees count, growth rel. growth rate,
  mort rate, reg trees/ha.
- **env_dt**: `siteID, year, <env vars…>` — must cover every siteID×year in
  obs_dt. Z-scale env vars; keep the scaling table to back-transform for plots.
- **init_trees** (→ makeInitCohorts): `siteID, patchID, year, species, dbh,
  trees, living` (year == 0 rows).
- **dist_dt** (optional): `siteID, year, intensity` (fraction of patches
  disturbed that year).

## Reference files (read, don't rewrite)
- `vignettes/Introduction-to-FINN.qmd` — source for vignette 1's mechanics.
- `vignettes/fit-fia-data/*.csv` — full Oregon subset (too big to ship as-is;
  subsample for vignette 3's bundled data).
- `playground/calibration/FIA/FIA_dataprep_vignette.R` — the real raw→FINN
  pipeline; basis for vignette 2's narrative.

## Verification loop (run before declaring done)
1. `quarto render vignettes/<file>.qmd` succeeds with no errors.
2. `R CMD build .` then `R CMD check --no-manual` is clean (vignettes build).
3. Eyeball the figures: succession curves are monotone-sensible; fit converges.
4. Confirm no absolute paths and no network/file reads outside the package.
See `dev/vignette-plan.md` for the full spec and acceptance criteria per vignette.
