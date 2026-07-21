## Submission

This is a new submission of FINN, a differentiable forest gap model.

## Test environments

* local macOS, R 4.5.0
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Possibly misspelled words in DESCRIPTION:
    Modularity (17:69)

  "Modularity" is spelled correctly.

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
