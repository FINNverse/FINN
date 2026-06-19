# Vignette build plan

Detailed spec for the three FINN vignettes. Read alongside `../CLAUDE.md`
(API + data contracts). Each section gives purpose, narrative arc, the chunks to
write, data needs, and acceptance criteria.

---

## Design decisions (made up front)

1. **All three are package vignettes**, not website-only articles — they ship in
   `vignettes/` and must run during `R CMD check`. Consequence: everything is
   small, fast, and self-contained.
2. **Heavy fitting is precomputed.** The FIA fit (vignette 3) is NOT trained at
   knit time. We train it once with a dev script, save the fitted `finn` object
   (`torch::torch_save()`), ship that object, and the vignette `torch_load()`s it
   to show parameters/fit. The training code is shown with `eval = FALSE`.
3. **Bundled data lives in `inst/extdata/`** and is read with
   `system.file("extdata", "<file>", package = "FINN")`. Never read from
   `vignettes/...` relative paths in a shipped vignette.
4. **Oregon FIA subset is subsampled to ~30–50 sites** so each shipped CSV is
   < ~300 KB. The full CSVs in `vignettes/fit-fia-data/` are the source.
5. **Data-prep vignette runs the FINN-side pipeline for real** on a tiny shipped
   raw `tree_dt` sample; the rFIA/WorldClim *acquisition* steps are shown with
   `eval = FALSE` (they need external downloads).

---

## Vignette 1 — `succession-demo.qmd`
**Title:** "Plausible succession from a handful of species"

### Purpose
Show that FINN, given only 3–5 fictional species with contrasting demographic
trade-offs, produces ecologically sensible succession — pioneers rise and fall,
shade-tolerant species take over late. This is the "wow, it just works" entry
point. No data, no fitting.

### The trade-off design (the scientific core)
Define 4 species along a classic life-history axis. Parameterize so the
trade-offs are explicit and visible in the output:

> **Shade parameter direction (verified against `R/Processes.R`).** In FINN the
> shade parameter (`parGrowth[,1]` / `parReg`) is the *light threshold* a species
> needs: growth/regeneration gate on a sigmoid centered at `light = shade_param`,
> so a species only succeeds where `light > shade_param`. Therefore **HIGH shade
> param = light-demanding = pioneer**, and **LOW shade param = shade-tolerant =
> climax**. (Earlier drafts of this table had it inverted.)

| Species | Strategy        | shade (light threshold) | growth rate | mortality | regeneration |
|--------:|-----------------|------------------------:|------------:|----------:|-------------:|
| 1 | pioneer            | high (0.60, needs light) | fast | high | high/early |
| 2 | early-mid          | mid–high (0.42)          | fast-ish | mid | high |
| 3 | mid-late           | mid (0.25)               | moderate | low | moderate |
| 4 | climax shade-tol.  | low (0.08)               | slow | very low | low but persistent |

Key levers (from Introduction-to-FINN.qmd parameter layout):
- `shadeSP` = light fraction a species needs (pioneers low, climax high).
- `parGrowth[,2]` size-dependent growth (pioneers higher).
- `parMort` baseline + size/growth dependence (pioneers higher baseline).
- `parReg` regeneration tied to shade; pioneers regenerate prolifically in light.
- Use a **constant environment** (`env1 = 0`) so succession is driven by the
  trade-offs, not the environment.

### Chunks to write
1. `setup` — `library(FINN); library(data.table); library(ggplot2); FINN.seed(...)`,
   set `Nsp=4`, `Ntimesteps=300`, `Nsites=1`, `patch_size=0.1`.
2. `species-parameters` — build `parReg, parGrowth, parGrowthEnv, parMort,
   parMortEnv, parRegEnv, parComp` per the table; print a tidy `pars_dt` so the
   reader sees the trade-offs. Add 1–2 sentences interpreting each parameter.
3. `environment` — constant `env_dt` (env1 = 0) and an optional light
   `dist_dt` (low disturbance, e.g. freq 0.02) to show gap dynamics.
4. `simulate` — `finn(...)` with fixed params (optimize flags OFF) then
   `m$simulate(init_cohort = NULL, env = env_dt, disturbance = dist_dt,
   patches = 100, device = "cpu")`. Use 100 patches for smooth means.
5. `plot-succession` — basal area (and trees) per species over time, faceted;
   reuse the Intro vignette's ggplot block. **This is the headline figure.**
6. `interpretation` (prose) — point out pioneer peak-then-decline and climax
   takeover; tie each pattern back to the trade-off parameters.
7. Optional `bare-ground vs disturbance` comparison (two panels).

### Acceptance criteria
- Knits in < ~30 s on CPU.
- Figure clearly shows successional turnover (pioneer early peak, climax late
  dominance). If it doesn't, tune `parMort`/`parReg`, don't fake it.
- Prose explains *why* in terms of the parameters.

---

## Vignette 2 — `data-preparation.qmd`
**Title:** "Preparing your data for FINN"

### Purpose
A standalone reference a user consults to get their own inventory data into
FINN. Documents the required tables, columns, units, and the builder pipeline.
Should be runnable on a tiny shipped raw sample so readers see real output.

### Narrative arc (mirrors playground/.../FIA_dataprep_vignette.R)
1. **What FINN needs** (prose + a diagram-ish list): three/four tables —
   `obs_dt`, `env_dt`, `init_trees`, optional `dist_dt`. Give the column
   contract table for each (copy from CLAUDE.md "Data formats").
2. **Start from a tree list.** Show the expected raw `tree_dt` shape:
   `siteName, patchName, treeName, year, species_name, dbh (cm), status
   (new/alive/dead), mort_cause, ...`. Ship a tiny example
   (`inst/extdata/example_tree_dt.csv`, a few hundred rows).
3. **Define status flags** — `reg` (status=="new"), `mort` (dead after alive),
   `living`. Show the `data.table` recipe.
4. **Build obs_dt** with `makeObsData(tree_dt, plotsize, fix_period_length,
   dbh_growth_thresh, NspeciesQuantile|Nspecies, Npatches, minNyears)`. Explain
   each argument and the **"other" species lumping** behavior. Run it; show
   `$obs_dt` head + `summary()`.
5. **Environmental data** — `env_dt` must cover every site×year; z-scale and keep
   a scaling table for back-transformed plots. Show the recipe (acquisition from
   WorldClim shown `eval=FALSE`; ship pre-extracted env values).
6. **Resolve IDs** — `resolveSiteIDs(tree_dt, obs_dt, env_dt)` to get integer
   `siteID/patchID/species` aligned across tables. Show before/after.
7. **Initial cohorts** — pull `year == 0` rows into `init_trees`, set
   `trees := 0` for non-living, then `makeInitCohorts(init_trees, Nspecies=)`.
   Show the resulting `CohortMat`.
8. **Checklist** (prose) — common pitfalls: inconsistent period length, init
   state left in obs_dt, env not covering all site×year, units (dbh cm, ba m²).

### Data needs
- `inst/extdata/example_tree_dt.csv` — small raw tree list (subsample of
  `vignettes/fit-fia-data/full_tree_dt.csv` or `init_trees.csv` reduced to a few
  sites, keeping the raw columns `status/mort_cause/...`).
- A few pre-extracted env rows so the env section runs without WorldClim.

### Acceptance criteria
- Every "how to" step that can run, runs on the shipped sample (makeObsData,
  resolveSiteIDs, makeInitCohorts all execute and print output).
- Acquisition-only steps (rFIA, terra/WorldClim) are `eval=FALSE` with a clear
  note that they need external data.
- Column-contract tables are correct and match `R/createInputs.R`.

---

## Vignette 3 — `fit-to-fia.qmd`
**Title:** "Fitting FINN to forest inventory data (Oregon, US FIA)"

### Purpose
Show the end-to-end calibration: load a prepared Oregon subset, build the FINN
object with optimizable processes, fit, and inspect fit quality + parameters.

### Narrative arc
1. **The data** — load the shipped subset (`obs_dt`, `env_dt`, `init_trees`,
   `species_dt`, `env_scales_dt`) via `system.file(...)`. Plot a couple of
   observed series (ba/trees by year) and the env distributions
   (back-transformed using `env_scales_dt`). Reuse blocks from the existing
   `fit-FIA.qmd`.
2. **Initial cohorts** — `makeInitCohorts(init_trees, Nspecies=...)`.
3. **Build Process-FINN** — `finn(N_species=..., recruits_dbh=12.9,
   competition=createProcess(~0, FINN::competition, optimizeSpecies=TRUE),
   growth/regeneration/mortality = createProcess(~., FINN::<f>,
   optimizeSpecies=TRUE, optimizeEnv=TRUE))`.
4. **Fit** — show the `m$fit(env, data, init_cohort, epochs, batchsize, patches,
   patch_size, weights, lr)` call with `eval=FALSE`; in the live chunk
   `torch_load()` a shipped pre-fitted model (`inst/extdata/fia_process_finn.pt`).
5. **Assess fit** — predicted vs observed (Spearman correlation per variable),
   convergence/loss curve if recorded.
6. **Inspect parameters / niches** — pull species×env effects; optionally `ALE()`
   to show an environmental response curve (the niche). This is the bridge to the
   paper's niche-inference aim.
7. **Simulate forward** from the fitted model (optional) and from bare ground.

### Data needs (ship in inst/extdata/, subsampled to ~30–50 sites)
- `fia_obs_dt.csv`, `fia_env_dt.csv`, `fia_env_scales_dt.csv`,
  `fia_init_trees.csv`, `fia_species_dt.csv`.
- `fia_process_finn.pt` — pre-fitted model (train via `dev/train_fia_model.R`,
  keep that script in `dev/`, not shipped). If the .pt is too large or
  torch-version-fragile, fall back to shipping a saved predictions data.frame and
  skip live model loading.

### Acceptance criteria
- Vignette knits without training (loads precomputed objects), < ~30 s.
- Shows at least one quantitative fit metric and one parameter/niche figure.
- All inputs via `system.file("extdata", ...)`; no `vignettes/` relative reads.

---

## Cross-cutting tasks
- [ ] Create `inst/extdata/` and the subsampled datasets; add a `dev/make_extdata.R`
      script that regenerates them from `vignettes/fit-fia-data/` (keep in `dev/`).
- [ ] Register all three in `_pkgdown.yml` (already stubbed) and confirm
      `VignetteBuilder: quarto` + Quarto CLI available.
- [ ] Add `^dev$` to `.Rbuildignore` (done) so dev scripts don't ship.
- [ ] Final: `R CMD check` clean with vignettes building.

## Open choices to confirm with Yannek
- Vignette 1 species count: 4 (default) vs 3 or 5.
- Whether to ship the pre-fitted `.pt` (preferred) or precomputed predictions
  only (more robust to torch version drift).
- Env variables to feature in the FIA fit (temp, prec are the safe pair).
