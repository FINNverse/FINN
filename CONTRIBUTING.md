# Contributing to FINN

**Workflow rules** (branch/PR, `git add` discipline, roxygenise-before-commit,
`dev/` layout, commit style) live in [CLAUDE.md](CLAUDE.md) — the single source of
truth. This file records **who contributed what**, and is maintained
automatically (see below), not by hand.

<!-- AUTO-MAINTAINED: the `project-progress` skill refreshes the table below from
     git history. Do not hand-edit; if something looks wrong, fix it in a session
     and it will be kept current. Last refreshed: 2026-11-08. -->

## Contributors

- **Yannek Käber** — `y.kaeber@posteo.de` (CRAN maintainer). git aliases: `ykaeber`, `Yannek Kaeber`.
- **Maximilian Pichler** — `maximilian.pichler@biologie.uni-regensburg.de`. git aliases: `MaximilianPi`, `dwjak123lkdmaKOP`, `Maximilian Pichler`.

## Contributions by theme

Best-effort attribution from git history; refined over time by the skill.

| theme / aspect | primarily | notes |
|---|---|---|
| Hybrid / neural-network processes (Level-1 & Level-2, transformer/DNN) | Maximilian | core of the MEE 2026 method paper |
| Core model & mechanistic processes (`finn-class`, `Processes`) | both | growth/mortality/regeneration/competition + forward pass |
| Calibration features (`env_autoscale`, `custom_parameters`, saturating regeneration, `recruit_obs_weight`, `growth_period_scale`, per-site aggregation windows for mixed remeasurement intervals) | Yannek | the 0.2.0 additions; per-site windows PR #17 |
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
| Collaboration / Claude Code setup | Yannek | shared CLAUDE.md, skills, hooks, claude-config |
