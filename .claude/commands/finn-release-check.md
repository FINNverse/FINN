---
description: Check FINN is ready for a CRAN release (consistency + R CMD check --as-cran)
---

Assess whether FINN is ready for a CRAN submission and report a checklist with
pass/fail for each item. This is a **read-only assessment** — do NOT bump the
version, edit files, commit, push, or submit to CRAN. Those are the maintainer's
explicit actions.

Check and report:

1. **Version.** `DESCRIPTION` `Version:` — state it, and whether it is greater than
   the current CRAN version (0.1.0). A release needs a bumped version.
2. **Changelog.** `NEWS.md` has a section for the DESCRIPTION version, and it
   describes the actual changes (not a placeholder).
3. **cran-comments.md** exists and describes this release (R CMD check result line,
   any expected NOTEs, downstream-dependency statement).
4. **Docs are current.** `Rscript -e 'roxygen2::roxygenise()'` produces **no**
   changes to `man/` or `NAMESPACE` (if it does, docs are stale — report which).
5. **Full check.** `R CMD build . && R CMD check --as-cran FINN_*.tar.gz`. Report
   the exact `errors | warnings | notes` line. CRAN wants 0 errors / 0 warnings;
   the misspelled-words NOTE and the torch/precompiled-vignette reviewer notes are
   expected.
6. **Tests.** `Rscript -e 'devtools::test()'` passes (or skips only on missing torch).

Finish with: **ready / not ready**, and if not ready, the specific blockers. Then
remind the maintainer that submission (`devtools::submit_cran()` or the web form)
is theirs to do.
