---
description: Regenerate docs and run the FINN test suite (optional full CRAN check)
argument-hint: "[full]"
---

Run FINN's standard pre-commit / pre-PR verification and report results
concisely. Do NOT commit, push, or submit anything — verify and report only, and
never use `git add -A`.

1. **Docs.** `Rscript -e 'roxygen2::roxygenise()'`. Then `git status --short` — if
   any `man/*.Rd` or `NAMESPACE` changed, list them and remind that they must be
   staged together with the code change that caused them.
2. **Tests.** `Rscript -e 'devtools::test()'`. Report the pass / fail / skip
   counts. Torch-gated tests skip without a backend — say so if they skip.
3. **Full CRAN check** — only if the argument is `full`:
   `R CMD build . && R CMD check --as-cran FINN_*.tar.gz`, then report the
   errors / warnings / notes (the misspelled-words NOTE is expected and benign).

End with a one-line verdict: is the working tree ready to commit?
