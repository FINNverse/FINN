# Contributing to FINN

**Workflow rules** (branch/PR, `git add` discipline, roxygenise-before-commit,
`dev/` layout, commit style) live in [CLAUDE.md](CLAUDE.md) — the single source of
truth. This file records **who contributed what**, and is maintained
automatically (see below), not by hand.

<!-- AUTO-MAINTAINED: the `project-progress` skill refreshes the table below from
     git history. Do not hand-edit; if something looks wrong, fix it in a session
     and it will be kept current. Last refreshed: 2026-09-02. -->

## Contributors

- **Yannek Käber** — `y.kaeber@posteo.de` (CRAN maintainer). git aliases: `ykaeber`, `Yannek Kaeber`.
- **Maximilian Pichler** — `maximilian.pichler@biologie.uni-regensburg.de`. git aliases: `MaximilianPi`, `dwjak123lkdmaKOP`, `Maximilian Pichler`.

## Contributions by theme

Best-effort attribution from git history; refined over time by the skill.

| theme / aspect | primarily | notes |
|---|---|---|
| Hybrid / neural-network processes (Level-1 & Level-2, transformer/DNN) | Maximilian | core of the MEE 2026 method paper |
| Core model & mechanistic processes (`finn-class`, `Processes`) | both | growth/mortality/regeneration/competition + forward pass |
| Calibration features (`env_autoscale`, `custom_parameters`, saturating regeneration, `recruit_obs_weight`, `growth_period_scale`) | Yannek | the 0.2.0 additions |
| CRAN packaging & release (DESCRIPTION, NEWS, cran-comments, review fixes) | Yannek | maintainer |
| Data interface (`makeObsData`/`resolveSiteIDs`/`makeInitCohorts`, extdata) | both | |
| Vignettes & docs | both | precompiled vignettes, pkgdown, README |
| Collaboration / Claude Code setup | Yannek | shared CLAUDE.md, skills, hooks, claude-config |
