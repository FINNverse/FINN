# FINN — shared Claude Code instructions

FINN is a hybrid dynamic forest (gap) model in R, built on `torch`. Demographic
processes (growth, mortality, regeneration, competition) can be fully
mechanistic, replaced by neural networks, or any mix, and are calibrated
end-to-end by gradient descent via `fit()`. This file is auto-loaded for **every**
Claude session in this repo, so it is the shared brain for all collaborators
(currently Yannek Käber and Max Pichler). Keep it committed and current — this is
how each person's Claude stays in sync. Personal, private notes go in
`~/.claude/` or a gitignored `CLAUDE.local.md`, never here.

## Shared collaboration setup (how the pieces fit — read this first)

Yannek and Max collaborate on this and related projects with Claude Code. The
setup spans this repo plus a shared config repo, all committed so a fresh session
discovers it automatically:

- **This repo** carries the auto-loaded brain: this `CLAUDE.md` (rules + map),
  `.claude/settings.json` + `.claude/hooks/` (shared permissions + a hook that
  blocks `git add -A`), repo-specific `.claude/commands/` (`/finn-check`,
  `/finn-release-check`) and `.claude/skills/` (`finn-process`, `finn-dataprep`,
  `finn-evaluation`, `finn-vignettes`), and `CONTRIBUTING.md` (skill-maintained
  record of who contributed what).
- **Shared config repo `FINNverse/claude-config`**, cloned wherever you keep your
  project repos — **the path differs per machine**, so nothing here hardcodes it
  and the skills resolve it themselves (`$CLAUDE_CONFIG_REPO`, their own symlink
  target, or a sibling of this repo). It holds cross-project **skills**
  — `project-progress` (auto-logs progress + refreshes
  `CONTRIBUTING.md`) and `doc-freshness` (periodically flags stale docs and asks
  you) — plus new-project templates and this project's **progress history** at
  `claude-config/projects/FINN/progress.md`.
- **First-time setup on a new machine** (each of us, once): clone `claude-config`
  next to your projects, then link the two shared skills either into
  `~/.claude/skills/` (available in every repo on the machine) or into this repo's
  `.claude/skills/` (scoped to this project; hide them with `.git/info/exclude`,
  not the tracked `.gitignore`). See the claude-config README. Without either link
  the cross-project skills won't load, but this repo's own commands/skills still do.
- **Do not hand-maintain** the progress history or `CONTRIBUTING.md` — the skills
  do. Volatile facts (versions, release status) are NOT stored here; read them
  from `DESCRIPTION`/`NEWS.md`.

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

## Where to find things (don't restate these here — they drift)

- **Current version / what changed** → `DESCRIPTION` (`Version:`) and `NEWS.md`.
  CRAN maintainer is Yannek Käber. Don't hard-code version/release status in this
  file; read it from those.
- **Process API & how to add one** → the `finn-process` skill (`createProcess`,
  custom processes, `custom_parameters`, the extend-never-mutate rule).
- **Data prep** → the `finn-dataprep` skill. **Validation/comparison** → the
  `finn-evaluation` skill. **Vignettes** → the `finn-vignettes` skill.
- **Per-project progress history** → `claude-config/projects/FINN/progress.md`
  (maintained by the `project-progress` skill). **Who contributed what** →
  `CONTRIBUTING.md` (also skill-maintained). Do not hand-maintain these.
- Anything you want the *other* collaborator's Claude to know must live in this
  file or another committed doc — personal Claude memory does not travel.
