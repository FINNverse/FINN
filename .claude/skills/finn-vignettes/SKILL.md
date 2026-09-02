---
name: finn-vignettes
description: The precompiled-vignette workflow for the FINN package. Use whenever creating, editing, rebuilding, or debugging a vignette under vignettes/ (the *.Rmd.orig sources and their committed *.Rmd output), or when a vignette change needs to appear in the built package or pkgdown site.
---

# FINN vignettes are PRECOMPILED — do not treat them as normal vignettes

The fits in FINN's vignettes need a `torch` backend, which is absent on CRAN and
most CI. So the vignettes are **precompiled**: the real code runs once locally, and
the knitted output is committed. They do **not** knit during `R CMD check`.

## The two-file pattern

- `vignettes/<X>.Rmd.orig` — the **source** with live code. Edit this.
- `vignettes/<X>.Rmd` — the **committed, knitted output** that ships. Generated;
  do not hand-edit.

## Workflow to change a vignette

1. Edit the `.Rmd.orig` source.
2. Re-knit with the build script, which needs a torch backend and a UTF-8 locale:
   ```bash
   Rscript vignettes/build.R          # re-knits the .orig -> committed .Rmd (+ figures)
   ```
   If characters like `Käber` come out mangled, ensure a UTF-8 locale
   (`LANG`/`LC_ALL` = `*.UTF-8`) before knitting.
3. Stage **both** the `.Rmd.orig` and the regenerated `.Rmd` (and any figures) in
   the same commit — specific paths only, never `git add -A`.

## Checks

- Do not rely on the vignette knitting during `R CMD check` — it uses the
  committed `.Rmd`. Verify your change by rendering locally, not by the check.
- Vignette code must use the **shipped API** (currently 0.2.0). If you changed the
  package API, update the `.Rmd.orig` and re-knit, or the committed output goes
  stale silently.
- pkgdown site: rebuilt/redeployed to the `gh-pages` branch separately from `main`
  (`docs/` is gitignored on `main`).
