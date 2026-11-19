# Contributing to FINN

## Workflow

1. **Stage specific paths** (`git add R/foo.R man/foo.Rd tests/...`), never
   `git add -A`, `git add --all` or `git add .`.
2. **Regenerate the docs before committing R changes:**
   `Rscript -e 'roxygen2::roxygenise()'`, and commit the resulting `man/*.Rd` and
   `NAMESPACE` with the code. Never edit `man/*.Rd` or `NAMESPACE` by hand.
3. **Bundle commits.** Few commits that each make one meaningful change, no
   fix-up commits; the body explains the why.
4. **CI must be green** (`R CMD check` on GitHub Actions) before a change lands
   on `main`. Pushing and CRAN submissions are the maintainer's call.
5. **Vignettes are precompiled:** edit `vignettes/*.Rmd.orig`, re-knit with
   `Rscript vignettes/build.R` (needs a torch backend) and commit the generated
   `.Rmd` and figures with it.
6. **Only package material in the repository.** Scratch scripts, experiments and
   personal configuration stay outside it (or untracked).

```bash
Rscript -e 'devtools::test()'                       # tests (torch-gated ones skip without a backend)
R CMD build . && R CMD check --as-cran FINN_*.tar.gz # full CRAN check
```

`dev/` holds the scripts that rebuild shipped material: `make_extdata.R` (the
bundled FIA datasets from `data-raw/`), `make_identity_snapshots.R` (the
reference values of `tests/testthat/test-identity.R`) and
`check_install_fresh.R`.

<!-- The table below is refreshed from git history; last refreshed 2026-11-08. -->

## Contributors

- **Yannek Käber** — `y.kaeber@posteo.de` (CRAN maintainer). git aliases: `ykaeber`, `Yannek Kaeber`.
- **Maximilian Pichler** — `maximilian.pichler@biologie.uni-regensburg.de`. git aliases: `MaximilianPi`, `dwjak123lkdmaKOP`, `Maximilian Pichler`.

## Contributions by theme

Best-effort attribution from git history.

| theme / aspect | primarily | notes |
|---|---|---|
| Hybrid / neural-network processes (Level-1 & Level-2, transformer/DNN) | Maximilian | core of the MEE 2026 method paper |
| Core model & mechanistic processes (`finn-class`, `Processes`) | both | growth/mortality/regeneration/competition + forward pass |
| Calibration features (`env_autoscale`, `custom_parameters`, saturating regeneration, `recruit_obs_weight`, `growth_rate`, per-site aggregation windows for mixed remeasurement intervals) | Yannek | the 0.2.0 additions; per-site windows PR #17 |
| Observation/loss API (`loss_family`, `loss_aggregation`: per-tree growth & mortality, size classes, dbh quantiles; `predictTrees`; persistent cohort ids) and `regeneration_adult` | Yannek | branch `loss-aggregation`, merged into `main` with `reanchor`/`teacher_forcing`/`loss_resolution` (design study of the training regimes: Max) |
| Trajectory regimes and loss resolution: the design study behind `reanchor` / `teacher_forcing` / `loss_resolution` and the cell-vs-per-tree loss axis | Maximilian | branch `fitting-interface`; the study itself lives outside this repo |
| Per-tree accounting: retention of simulated-dead tracked cohorts, the forced per-site path under a per-tree term, `$tree_coverage` → `$tree_deaths`, `createProcess(light_mask = )` | Maximilian | branch `fitting-interface` |
| Stochastic draws: `createProcess(mortality_draw = , regeneration_sampler = )` and `createProcess(growth_scale = )` | Maximilian | branch `fitting-interface`; opt-in, defaults unchanged |
| `fit()` set-up performance (O(n) site/year lookup, tree tensors on the trees' own patches, `epochs = 0`) and `mortality(growth_scale = 0.01)` | Yannek | branch `fast-setup` |
| CRAN packaging & release (DESCRIPTION, NEWS, cran-comments, review fixes) | Yannek | maintainer |
| Locale-stable ids: radix ordering in `resolveSiteIDs()`, and the regression test behind it | Maximilian | the fourth change that does not reproduce 0.1.0 off a `C` locale |
| Per-site loss pooling: one mean per batch instead of a sum of per-year means | Yannek | FINNverse/FINN#29, partial — the unit-count half |
| 0.2.0 release readiness: `R CMD check --as-cran` clean-up, the argument-validation batch, `.Rbuildignore`, and the two test repairs | Maximilian | PR #32 review round |
| Data interface (`makeObsData`/`resolveSiteIDs`/`makeInitCohorts`, extdata) and the bundled `fia_*` datasets | both (datasets: Yannek) | FIA sample as lazy-loaded data |
| Vignettes & docs | both | precompiled vignettes, pkgdown, README |