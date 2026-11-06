## Submission

This is a minor-version update of FINN (0.1.0 -> 0.2.0). It adds model features
developed for two forest-inventory calibration projects: per-response
likelihoods and loss aggregations including per-tree growth and mortality,
multiple shooting and teacher forcing as training regimes, `predictTrees()`,
saturating and adult-limited regeneration, per-site aggregation windows for
mixed remeasurement intervals, and the Oregon FIA calibration sample as
lazy-loaded datasets in place of the CSVs under `inst/extdata`.

Three defaults do not reproduce 0.1.0, each with an argument that restores it:
`mortality(growth_scale = 0.01)` (`growth_scale = 1`), `createProcess(light_mask
= TRUE)` (`light_mask = FALSE`) and `finn(reg_floor = 1e-3)` (`reg_floor =
0.2`). They are corrections rather than preferences: the unmasked height
quantiles made a stand's light depend on which other sites shared its training
batch, the old recruitment floor alone exceeded observed ingrowth by more than
an order of magnitude on a 65-species pool, and the raw growth scale left the
mortality growth coefficient at its bound for most species. Each is documented
in NEWS.md with the measurement behind it. Every other new argument defaults to
the previous behaviour.

A fourth change is locale-dependent rather than argument-controlled:
`resolveSiteIDs()` numbered sites, patches and species through `LC_COLLATE`, so
the same input gave different integer ids on machines with different locales,
and the species id indexes the fitted parameter matrices. Ordering is now
radix. Results are unchanged under `C` and reproducible everywhere else.

New in 0.2.0 (see NEWS.md for details):

* `env_autoscale` gains per-predictor scaling ("auto"/"identity"/"0to1" or a
  user function); the previous TRUE/FALSE/NULL forms are unchanged.
* `finn()` gains `regeneration_saturation` (a Beverton-Holt cap on recruitment,
  read back with the new exported `reg_saturation_K()`), `recruit_obs_weight`
  and `growth_period_scale`.
* `createProcess()` gains a `custom_parameters` interface.
* A regeneration observation operator (internal hook).
* Bug fix: `createProcess(NN = ...)` is no longer dropped in `create_nn()`.

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
  surnames. The DESCRIPTION cites the accompanying methods paper as Pichler and
  Käber (2026) <doi:10.1111/2041-210x.70347>, which resolves.

## Notes for the reviewer

* FINN builds on the 'torch' package, whose backend (libtorch) is downloaded at
  runtime and is not present on the check machines. All examples that need it are
  guarded with `@examplesIf torch::torch_is_installed()`, and the tests that need
  it skip on CRAN, so the check runs without the backend.
* The vignettes are precompiled (the model-fitting code is run once by the
  authors and the results are committed), so building them does not require the
  torch backend or a long run time.

## Downstream dependencies

There are no reverse dependencies on CRAN at this time.
