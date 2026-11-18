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

* win-builder, R-devel (2026-11-17 r90926 ucrt) and R-release (4.6.2 ucrt).
  The run before the final fixes gave one NOTE on both, "height/width
  attributes should be in pixels" for the figure in `?finn`, now given in
  pixels; and on R-devel a URL NOTE for https://www.gnu.org/licenses/gpl-3.0
  (README.md), a connection timeout to www.gnu.org from the check machine; the
  URL resolves.

## R CMD check results

Local run with the torch backend installed, `R CMD check --as-cran --no-manual`
(the machine has no pdflatex; the Rd files pass every Rd check):

0 errors | 1 warning | 0 notes

* checking whether package 'FINN' can be installed ... WARNING
  `R_ext/Boolean.h:62:36: warning: unknown warning group
  '-Wfixed-enum-extension'`. This comes from R 4.5.0's own header with Apple
  clang 21 on the local machine, not from FINN's code; it is absent on the
  check services' compilers.
* On a loaded machine the check notes examples above 5 s elapsed
  (`feature_importance`, `ALE` and a few more, each fitting a small model on
  two or three sites of the bundled data first). They need the torch backend
  and do not run on the check machines (`@examplesIf
  torch::torch_is_installed()`).

The test suite ran in full with the backend: 909 tests passed, 0 failed, 1
skipped (a thread-count test that torch's native macOS backend cannot repeat
in one session). The tarball is 4.2 MB, the installed package under 5 MB.

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
