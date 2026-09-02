# Contributing to FINN

For Yannek, Max, and their Claude sessions. Machine-facing rules live in
[CLAUDE.md](CLAUDE.md); this is the human-facing workflow.

## Branch & PR model

- `main` is the integration branch and is protected. **Never commit directly to
  it** once it is on `origin`; open a pull request.
- One branch per feature/fix, named `topic-short-desc` (e.g. `pft-soft-membership`,
  `fix-reg-saturation`). Keep branches short-lived and PRs small — AI generates
  large diffs fast, and long-lived branches make merges painful.
- We do **not** use branch protection, so this is convention, not enforcement:
  merge a PR only after **CI (`R CMD check`) is green** and, ideally, the other
  person has glanced at it. Run `/finn-check` (or `R CMD check`) locally before
  merging — with no protection, that local discipline is the only guard against
  either of us pushing a broken `main`.
- Coordinate ownership via GitHub Issues so you are not editing the same files at
  once. Rough split: assign each feature an issue and an owner before starting.

## Before every commit

1. If you touched roxygen/exports: `Rscript -e 'roxygen2::roxygenise()'` and stage
   the regenerated `man/*.Rd` + `NAMESPACE` **in the same commit**.
2. `Rscript -e 'devtools::test()'` passes locally (or you know why a torch test skips).
3. Stage **specific paths** — never `git add -A` (a hook blocks it anyway).
4. Commit as yourself, imperative subject, body explains *why*, no AI co-author trailer.

## Generated / conflict-prone files

`man/*.Rd`, `NAMESPACE`, and precompiled vignette `.Rmd` are generated. **Do not
hand-merge them** — on a conflict, take either side, then re-run
`roxygenise()` / the vignette build and commit the regenerated result. `NEWS.md`
uses `merge=union` (see `.gitattributes`) so both entries survive a merge.

## The `dev/` directory

`dev/` is shared but namespaced. **Put your own scratch in `dev/<yourname>/`**
(`dev/yannek/`, `dev/max/`). Shared tooling (data builders, plans) stays at the
`dev/` top level. Nothing in `dev/` ships in the package — it is listed in
`.Rbuildignore`. Never let `dev/` content into an `R/`, `man/`, or release commit.

## Releases (maintainer only)

CRAN submission and pushing release tags are the maintainer's (Yannek's) call.
Checklist for a release PR: bump `DESCRIPTION` Version, update `NEWS.md` and
`cran-comments.md`, `roxygenise()`, `R CMD check --as-cran` clean, then the
maintainer submits.

## Claude Code setup

The committed `CLAUDE.md` and `.claude/settings.json` configure both our Claude
sessions identically. Personal overrides go in gitignored `CLAUDE.local.md` /
`.claude/settings.local.json`. Cross-project templates and shared slash commands
live in the separate `claude-config` repo (symlink its `shared/` into `~/.claude/`).
