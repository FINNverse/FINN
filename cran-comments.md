## Submission

This is a minor-version update of FINN (0.1.0 -> 0.2.0). It adds the model
features developed for two forest-inventory calibration projects and the
interface agreed for them:

* `finn()` writes its default processes in the signature, `growth_rate` says
  what the observed growth column is, and the recruitment floor is a setting of
  the regeneration process (`createProcess(reg_floor = )`).
* `fit()` groups its loss arguments (`loss_family`, `loss_aggregation`,
  `loss_resolution`, `loss_weights`, `weight_scale`, `loss_recruit`), collects
  the optimizer settings in `finnControl()`, checks its input tables before
  anything runs, and has the help pages `?finn_data` and `?finn_loss`. Every
  renamed or moved argument of 0.1.0 is still accepted with a warning that
  names the new form.
* Per-response likelihoods and aggregations including per-tree growth and
  mortality (`tree_data`, `predictTrees()`), multiple shooting and teacher
  forcing as training regimes, per-site aggregation windows for mixed
  remeasurement intervals, the regeneration processes `regeneration_capped()`
  (the new default), `regeneration_adult()` and `regeneration_adult_effect()`,
  a `createProcess(custom_parameters = )` interface, per-predictor
  environmental scaling, `fit(threads = )`, and `device = "cuda"` as a synonym
  for `"gpu"`.
* Data preparation on annual simulation steps: `makeObsData()` returns yearly
  rates and `period_length`, `resolveSiteIDs()` an annual time axis; the Oregon
  FIA calibration sample ships as lazy-loaded datasets (`fia_obs_dt` and eight
  more, including the per-tree tables `fia_tree_dt` / `fia_tree_test`) on that
  axis, in place of CSV files under `inst/extdata`.
* Evaluation and plotting of a fitted model: `evaluate()` compares simulations
  with the observations of the training data and a holdout and returns the
  result (also kept in the model, also run by `fit(holdout = )`);
  `comparePlot()` draws observed against simulated (scatter page, heatmap per
  species, time series); `summary()` of a model reports the skill, the
  convergence of the fit and the environmental effects; `plot()` methods for a
  model (the fitted parameters), a simulation and an evaluation; and
  `simulateForest(bare_ground = TRUE)` runs a model from an empty stand. The
  base package 'grid' is added to Imports for arranging the figure panels.
* `fit()` compounds the simulated yearly rates of an observation interval into
  the per-year rate the observation carries (the functional `makeObsData()`
  uses), and the mortality response starts at about 2 % per year instead of
  torch's random start.

Defaults that do not reproduce 0.1.0, each with the value that restores it, are
listed in the two tables at the top of the 0.2.0 section of NEWS.md with the
measurement behind each: `mortality(growth_scale = 0.01)` (`1`),
`createProcess(light_mask = TRUE)` (`FALSE`), `createProcess(reg_floor = 1e-3)`
(`0.2`), `fit(trees_offset = 0.01)` (`NULL`), `fit(rb_gradient = TRUE)`
(`FALSE`), the binomial death draw and the refined recruit draw (`"gamma"`,
`"gaussian"`), the mortality limits `0.001`/`0.999` (`0.01`/`0.99`), the recruit
dispersion floor `theta_floor = 0.5` (`0`), and `regeneration_capped()` as the
default regeneration process (`FINN::regeneration`). `resolveSiteIDs()` numbers
sites, patches and species by radix order instead of the collation locale, so
results are unchanged under `C` and reproducible everywhere else; on its census
axis it now builds the initial cohorts from the first inventory.

## Test environments

* local macOS 26.5 (aarch64), R 4.5.0, `R CMD check --as-cran` on the built
  tarball, twice: with the torch backend (libtorch) installed, so every example
  and test ran, and with the backend hidden, as on the check machines (examples
  and backend tests skip). Both: 0 errors, the one warning below, and a NOTE
  "unable to verify current time" (no time server reachable from the machine).
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macos-latest and windows-latest (R-release), without the backend: OK.

* win-builder (R-devel and R-release) on the submitted state: see the results
  added below once they are in.

## R CMD check results

Local run with the torch backend installed, `R CMD check --as-cran` (PDF
manual included):

0 errors | 1 warning | 0 notes (besides the time-server and HTML Tidy notes of
the local machine)

* checking whether package 'FINN' can be installed ... WARNING
  `R_ext/Boolean.h:62:36: warning: unknown warning group
  '-Wfixed-enum-extension'`. This comes from R 4.5.0's own header with Apple
  clang 21 on the local machine, not from FINN's code; it is absent on the
  check services' compilers.
* With the backend the check notes one example above 5 s elapsed (the
  package help page `?FINN`, which fits a small model on four sites of the
  bundled data, evaluates it and draws the results, 7.6 s). It needs the torch
  backend and does not run on the check machines (`@examplesIf
  torch::torch_is_installed()`).

The test suite ran in full with the backend: 1061 tests passed, 0 failed, 1
skipped (a thread-count test that torch's native macOS backend cannot repeat
in one session). Without the backend the check is clean apart from the local
compiler warning above. The tarball is 2.7 MB.

## Notes for the reviewer

* FINN builds on the 'torch' package, whose backend (libtorch) is downloaded at
  runtime and is not present on the check machines. Every example that needs
  the backend is guarded with `@examplesIf torch::torch_is_installed()`; there
  is no `\dontrun{}`. The test files that need the backend skip when it is
  absent (`skip_if_no_torch()`), so the check runs without it.
* The vignettes are precompiled: the model-fitting code is run once by the
  authors and the knitted output is committed, so building them needs neither
  the torch backend nor a long run time.
* DESCRIPTION: 'DNNs' (deep neural networks), 'FINN' (the package and method
  acronym), 'torch' (the R package FINN is built on) and the two author
  surnames are spelled as intended. The methods paper cited, Pichler and Käber
  (2026) <doi:10.1111/2041-210x.70347>, resolves.

## Downstream dependencies

There are no reverse dependencies on CRAN at this time.
