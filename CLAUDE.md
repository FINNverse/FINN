# FINN — shared Claude Code instructions

FINN is a hybrid dynamic forest (gap) model in R, built on `torch`. Demographic
processes (growth, mortality, regeneration, competition) can be fully
mechanistic, replaced by neural networks, or any mix, and are calibrated
end-to-end by gradient descent via `fit()`. This file is auto-loaded for **every**
Claude session in this repo, so it is the shared brain for all collaborators
(currently Yannek Käber and Max Pichler). Keep it committed and current — this is
how each person's Claude stays in sync. Personal, private notes go in
`~/.claude/` or a gitignored `CLAUDE.local.md`, never here.

## Golden rules (do not violate)

1. **Never `git add -A`, `git add --all`, or `git add .`** in this repo. It has
   swept unpublished manuscripts and personal scratch into release commits.
   Always stage specific paths (`git add R/foo.R man/foo.Rd tests/...`). A
   PreToolUse hook enforces this.
2. **Regenerate docs before committing R changes.** After editing any roxygen,
   run `Rscript -e 'roxygen2::roxygenise()'` and stage the resulting `man/*.Rd`
   and `NAMESPACE` in the *same* commit. Never hand-edit `man/*.Rd` or `NAMESPACE`.
3. **Do not push to `origin` or submit to CRAN without the maintainer's explicit
   go-ahead.** Outward-facing actions are the maintainer's call (Yannek is the
   CRAN maintainer). Local commits are fine; pushing is not, unless asked.
4. **Branch per feature; merge to `main` via PR.** Never commit directly to a
   shared `main` once it is pushed. CI (`R CMD check`) must be green before merge.
5. **Keep personal scratch in your own `dev/` subfolder** (`dev/yannek/`,
   `dev/max/`) — never at the repo root and never in `R/`, `tests/`, `man/`,
   `vignettes/`, or a package commit. See CONTRIBUTING.md.
6. **Commit as yourself** (your own git `user.name`/`user.email`). Team
   convention: no AI co-author trailer; imperative subject line; the body
   explains the *why*.

## Commands

```bash
# run the test suite (torch-gated tests skip without a backend)
Rscript -e 'devtools::test()'
# regenerate man/ + NAMESPACE from roxygen (ALWAYS before committing R changes)
Rscript -e 'roxygen2::roxygenise()'
# full CRAN check (what CI runs; do this before a release PR)
R CMD build . && R CMD check --as-cran FINN_*.tar.gz
# load the package from source for a quick interactive check
Rscript -e 'pkgload::load_all("."); ...'
```

Vignettes are **precompiled** (a torch backend is needed to knit the fits): edit
the `vignettes/*.Rmd.orig`, then re-knit with `vignettes/build.R` and commit the
generated `.Rmd`. Do not rely on vignettes knitting during `R CMD check`.

## Repository layout

- `R/` — package source (13 files). Core: `finn-class.R` (the `finn()`/`fit()`
  model + forward pass), `Processes.R` (the mechanistic process functions:
  `growth`, `mortality`, `regeneration`, `regeneration_saturation`,
  `competition`, and their `_hybrid` variants), `Processes_utils.R`
  (`createProcess`, `createHybrid`, parameter init/bounds, env scaling),
  `xAI.R` (ALE, feature importance), `createInputs.R` (data builders).
- `tests/testthat/` — torch-gated with `skip_if_no_torch()`; shared fixtures in
  `helper.R` (`fit_toy_finn()`).
- `man/`, `NAMESPACE` — **generated** by roxygen; never edit by hand.
- `vignettes/` — `*.Rmd.orig` sources + committed knitted `*.Rmd` (precompiled).
- `inst/extdata/` — small shipped datasets (FIA subset) used by tests/vignettes.
- `dev/` — scratch and tooling; **personal work goes in `dev/<name>/`**. Shared
  dev tooling (`make_extdata.R`, plans) stays at `dev/` top level.

## Key API facts (verify against source before relying on them)

- `finn(N_species, growth_process=, mortality_process=, regeneration_process=,
  competition_process=, recruits_dbh=, recruit_obs_weight=, growth_period_scale=)`.
- `createProcess(formula, func, optimizeSpecies=, optimizeEnv=, initSpecies=,
  initEnv=, custom_parameters=, ...)`. A custom `func` is bound to `self` and can
  declare extra trainable parameters via `custom_parameters` (registered,
  optimized by `fit()`, saved/loaded).
- `createHybrid(...)` replaces a whole process with a DNN.
- Density-dependent recruitment is the `FINN::regeneration_saturation` process
  (Beverton-Holt, `K = exp(reg_logK)` declared via `custom_parameters`), a drop-in
  alternative to `FINN::regeneration` — **not** a `finn()` option.

## Release state

- **CRAN: 0.1.0** (published; maintainer Yannek Käber). **`main` is 0.2.0,
  pushed to GitHub but not yet submitted to CRAN.** 0.2.0 changes are additive/backward-compatible
  (per-predictor `env_autoscale`, `custom_parameters`, saturating regeneration,
  `recruit_obs_weight`, `growth_period_scale`). See `NEWS.md`.
- Anything you want the *other* collaborator's Claude to know must live in this
  file (or a committed doc) — personal Claude memory does not travel.
