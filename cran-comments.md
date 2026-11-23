## Submission

This is a bug-fix update of FINN (0.2.0 -> 0.2.1), submitted soon after 0.2.0.
Since the release we found defects that either stop long model fits or silently
affect results:

* `fit()` never ran R's garbage collector, so the torch autograd graphs of
  several epochs piled up until R happened to collect. Long fits with hybrid
  (neural-network) processes on large data sets ran out of memory hours in
  (37 to 85 GB in five epochs on one data set; with the fix flat at 38-39 GB).
  `fit()` now collects once per epoch.
* `makeObsData()` with a species cap dropped the mortality counts of the
  species lumped into "other".
* Optional `fit()` targets (diameter quantiles, size-class counts) with a
  missing value scored later values against the wrong site.
* On the per-site aggregation path a census step without an observation window
  was skipped, so its stand-state terms were not scored.
* A second `fit()` call after a learning-rate schedule kept the rate the
  schedule had left instead of the `lr` it was given.

Each fix is a few lines with a regression test. No argument, default or
default result changes; the fixes change numbers only on the narrow paths
named above (listed in NEWS.md).

## Test environments

* local macOS 26.5 (aarch64), R 4.5.0, `R CMD check --as-cran` on the built
  tarball with the torch backend (libtorch) installed, so every example and
  test ran; the PDF manual built separately with `R CMD Rd2pdf`.
* win-builder R-devel (2026-11-17 r90926 ucrt) and R-release (4.6.2 ucrt):
  1 NOTE each, "Days since last update: 4".

## R CMD check results

0 errors | 0 warnings | 1 note (besides the local machine's notes below)

* "Days since last update: 4": this is the bug-fix release described above.

On the local machine, as for 0.2.0: a WARNING from R 4.5.0's own header
`R_ext/Boolean.h` with Apple clang 21 (`unknown warning group
'-Wfixed-enum-extension'`), not from FINN's code; and NOTEs that the time
server is unreachable and that HTML Tidy is too old. With the backend the
examples of `?FINN` and `?evaluate` use more than 5 s of CPU time (8.8 s and
7.0 s, under 5 s elapsed); they need the backend and do not run on the check
machines (`@examplesIf torch::torch_is_installed()`).

The test suite ran in full with the backend: 1269 expectations passed, 0
failed, 1 skipped (a thread-count test that torch's native macOS backend
cannot repeat in one session). The tarball is 2.7 MB.

## Notes for the reviewer

* FINN builds on the 'torch' package, whose backend (libtorch) is downloaded at
  runtime and is not present on the check machines. Every example that needs
  the backend is guarded with `@examplesIf torch::torch_is_installed()`; the
  test files that need it skip when it is absent.
* The vignettes are precompiled, as in 0.2.0.

## Downstream dependencies

There are no reverse dependencies on CRAN at this time.
