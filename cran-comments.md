## Resubmission

This is a resubmission addressing the points raised by Konstanze Lauseker:

* Replaced all `T`/`F` with `TRUE`/`FALSE` in the R code (and thus in the
  generated `\usage`), and confirmed `T`/`F` are not used as names.
* Removed the example from the unexported `extract_env()` (now internal, `@noRd`).
* Replaced `\dontrun{}`: examples that run in < 5 s are now unwrapped and
  executable; examples that require the (optional, runtime-downloaded) 'torch'
  backend use `@examplesIf torch::torch_is_installed()`.
* No function writes to the user's home/package/working directory by default:
  the internal figure-generation helper that wrote to `man/figures/` has been
  moved out of the package (to `data-raw/`, which is build-ignored).
* `options(na.action=)` in `R/finn-class.R` is now restored with an immediate
  `on.exit()` (the earlier graphical-parameter capture already used `on.exit()`).

## Submission

This is a new submission of FINN, a differentiable forest gap model.

## Test environments

* local macOS, R 4.5.0
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Possibly misspelled words in DESCRIPTION:
    DNNs, FINN, Kaeber (Käber), Pichler, torch

  These are correct: 'DNNs' (deep neural networks), 'FINN' (the package /
  method acronym), 'torch' (the R package FINN is built on), and the two author
  surnames. The DESCRIPTION also cites the accompanying methods paper as
  Pichler and Käber (2026) <doi:10.1111/2041-210x.70347>, which resolves.

  On win-builder the same NOTE also reported a timeout reaching
  <https://www.gnu.org/licenses/gpl-3.0> (the "License: GPL v3" badge target in
  README.md). The URL is valid and reachable; this was a transient network
  timeout on the check host.

## Notes for the reviewer

* FINN builds on the 'torch' package, whose backend (libtorch) is downloaded at
  runtime and is not present on the check machines. All examples that need it are
  guarded with `@examplesIf torch::torch_is_installed()`, and the tests that need
  it skip on CRAN, so the check runs without the backend.
* The vignettes are precompiled (the model-fitting code is run once by the
  authors and the results are committed), so building them does not require the
  torch backend or a long run time.

## Downstream dependencies

There are none, as this is a new submission.
