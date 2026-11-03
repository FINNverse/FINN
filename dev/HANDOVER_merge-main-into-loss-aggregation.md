# Handover: branch `merge-main-into-loss-aggregation`

For Max (and the Claude session he works with). Written 2026-11-03 by Yannek's
session. Read this before touching the branch; it explains what the merge
decided, what was verified, and what is still open.

## TL;DR

- **Branch:** `merge-main-into-loss-aggregation` @ `93f8cdd`, one merge commit
  with parents `6497d01` (`loss-aggregation`, Yannek) and `d4b3907` (`main`, incl.
  PR #21 multishoot-fitting). `main` and `loss-aggregation` are untouched.
  Merging this into `main` is planned but **not done**; we want your review first.
- **Tests:** `devtools::test()` gives 384 passed, 0 failed, 0 skipped (10 warnings,
  all the pre-existing `loss` deprecation warning in `test-per-site-period.R`).
  `R CMD check` has **not** been run on the merge.
- **Reproducibility:** bit-for-bit identical to `6497d01` on the default fit, the
  per-site path, Yannek's BWI stage-2 call and `predictTrees()`; bit-for-bit
  identical to `main` (`d4b3907`) with `nonfinite = "skip_step"`, incl.
  `loss_resolution = "tree"` with `reanchor` / `teacher_forcing`. Details below.
- **The one interface change you will notice:** `loss_resolution = "tree"` now
  takes `tree_data` as a **table** (with a `cohortID` column), not as a list of
  per-timestep arrays. `fit()` builds your arrays from it internally.

## Quick start

```bash
git fetch origin
git switch merge-main-into-loss-aggregation
```

```r
devtools::document()
devtools::test()          # expect 384 pass
?fit                      # section "Scoring individual trees"
```

Install into a separate library if you have fits running against your default one:

```bash
mkdir -p ~/R/lib-finn-merge
R CMD INSTALL -l ~/R/lib-finn-merge .
```

## What overlapped and how it was resolved (Yannek's decisions)

Three features had been built on both branches; one shared code path behaved
differently. The rule was "main's version wins where it is the same feature";
where the two turned out to be different methods, both were kept as explicit
options.

| # | Feature | Resolution | Why |
|---|---|---|---|
| 1 | Multiple shooting: `main` `reanchor`/`anchor_cohorts`/`anchor_mode`/`anchor_min_dbh`/`teacher_forcing` vs `loss-aggregation` `shooting = "multiple"` | **main's kept**, `shooting` removed completely (argument, `private$anchor_state`, `acc_g`/`acc_s`, all `multiple` branches, `test-multiple-shooting.R`) | same idea; main's is more complete ("merge" mode, teacher forcing, more tests); no BWI script used `shooting` |
| 2 | Per-tree growth/mortality loss: `loss_aggregation = c(growth = "tree", mortality = "tree")` vs `loss_resolution = "tree"` | **both kept**, chosen explicitly, documented side by side in `?fit` | different estimators: the first scores the kernels on each OBSERVED tree (measured dbh, light from the observed stand; responses separately; works with binned cohorts), the second follows each SIMULATED cohort along the model's trajectory (growth and mortality together; one cohort per tree; combines with reanchor/teacher forcing). Yannek's final BWI model was trained with `loss_aggregation = c(mortality = "tree")` |
| 3 | `tree_data` argument (same name, two formats) | **always a table**; the option is never inferred from the input type. For `loss_resolution = "tree"` a helper converts the table to main's per-timestep arrays | one consistent user-facing format; FINN's other inputs are tables too; a misaligned hand-built array was scored silently before |
| 4 | Cohort ids: main/base flat ids (new ids from `max(cohort_ids) + 1` after pruning) vs per-(site, patch) ids from a monotone counter | **loss-aggregation's ids kept**; main's tree-loss slot index remapped onto them | main's scheme can re-issue a pruned tree's id to a recruit, which `loss_resolution = "tree"` would then score against the dead tree's target; ids never enter the default loss |
| 5 | Non-finite loss terms: drop only that term (loss-aggregation) vs skip the whole timestep (main) | **option** `fit(nonfinite = c("drop_term", "skip_step"))`, default `"drop_term"` | `"drop_term"` reproduces all existing BWI fits; `"skip_step"` reproduces main exactly |
| 6–9 | main only: FIA datasets as package data, `update_step` "scores no observation" check, CITATION/README/NEWS/.Rbuildignore, annual-rate doc wording | kept | – |
| 10–16 | loss-aggregation only: `loss_family` (`loss` = deprecated alias), `size_classes`, `quantiles`, `predictTrees()`, `regeneration_adult`, optimizer rebuild on resume (efe9512), `createProcess` upper/lower fix (6497d01) | kept | – |

## The interface now

```r
fit(model, data, env, ...,
    loss_family      = c(dbh = "mse", ba = "mse", trees = "poisson",
                         growth = "mse", mortality = "binomial", regeneration = "nbinom"),
    loss_aggregation = NULL,          # e.g. c(mortality = "tree"), list(dbh = list(type = "quantiles", ...))
    tree_data        = NULL,          # ALWAYS a data.frame/data.table
    nonfinite        = c("drop_term", "skip_step"),
    reanchor = FALSE, anchor_cohorts = NULL, anchor_mode = "replace", anchor_min_dbh = NULL,
    teacher_forcing  = FALSE,
    loss_resolution  = c("site", "tree"),
    ...)
```

`tree_data`: one row per tree and census; `year` = the census that CLOSES the
tree's interval. Always `siteID, year, patchID, species, growth, died`; then
- `loss_aggregation = "tree"`: additionally `dbh` at the interval start (optionally `trees`);
- `loss_resolution = "tree"`: additionally `cohortID`, the initial cohort the tree is
  (from `makeInitCohorts(..., treeID_table = TRUE)$init_trees`, unbinned, one row = one
  cohort); `patchID` is the layout patch index; `species` and `dbh` are checked against
  the cohort. Unknown, duplicate or mismatching cohorts, binned initial cohorts, missing
  columns and rows outside the scored years stop with an informative error.
- Combining `loss_resolution = "tree"` with `loss_aggregation` `"tree"` on growth or
  mortality is an error.

**Your existing code** that passes per-timestep `growth`/`mort` arrays to
`loss_resolution = "tree"` must switch to the table. Your tests were converted
this way (see `test-tree-loss.R`); padding slots now carry NaN targets because a
table can only name real trees.

## Where things are (R/finn-class.R on the branch)

| what | where |
|---|---|
| `?fit` section *Scoring individual trees* (comparison table) | roxygen, ~line 473 |
| exported `fit()` signature (`nonfinite`, `loss_resolution`) | ~line 520–540 |
| `predictTrees()` | ~line 608 |
| `forward()`: `observed_state()` (anchor ids from the counter), tree-loss slot remap (`tree_sp_off`, `tree_sp_n`, `tree_gen_nc`) | ~line 940–1010, 1220–1235 |
| `private$scored_terms()`, `private$loss_at_step()` (the `nonfinite` switch) | ~line 2713–2733 |
| `private$resolve_loss_aggregation()` | ~line 2734 |
| `private$tree_predictions()` (observed-tree scoring) | ~line 2880 |
| `private$build_tree_targets()` (table → main's arrays) | ~line 2940 |
| `private$build_tree_obs()` (table → observed-tree tensors) | ~line 3022 |

## Verification that was done

Two temporary libraries (6497d01 / d4b3907 vs the merge), separate Rscript
processes, same seeds, `shuffle = FALSE`, `identical()` on `$history`,
`loss_baseline`, `loss_weights` and every parameter after a few epochs with lr > 0:

| reference | scenario | result |
|---|---|---|
| 6497d01 | default constant-period fit (FIA, 4 sites, 3 epochs) | identical |
| 6497d01 | `period_mode = "per_site"` cell path (mixed periods) | identical |
| 6497d01 | regeneration-only weights | identical |
| 6497d01 | `loss_aggregation` tree + size classes | identical |
| 6497d01 | `predictTrees()` on that model | identical |
| 6497d01 | Yannek's BWI stage-2 call on the BWI pilot set (1999 sites; probe fit for weights, custom Adam with L2 pull, `loss_aggregation = c(mortality = "tree")`, per-site) | identical with `torch_set_num_threads(1)`; with default threads parameters differ by ≤ 3e-8, the same as two runs of 6497d01 itself (torch CPU reduction order) |
| d4b3907 with `nonfinite = "skip_step"` | default fit | identical |
| d4b3907 | `loss_resolution = "tree"`: single shooting, reanchor replace, reanchor merge, teacher forcing (main's arrays vs the equivalent table) | identical |
| d4b3907 | reanchor merge / teacher forcing without tree terms | identical |

## Deviations from main's code you should know about

- Main's id-based assertions in `test-reanchor.R` (merge tests) and
  `test-teacher-forcing.R` (`peak_id`) assumed globally unique flat ids; they now
  key on `(site, patch, id)`.
- One assertion in `test-reanchor.R` ("every cohort the anchor can see is gone")
  already failed on `main` itself (emptied padding slots counted as merged); it
  now checks the cohorts handed over alive.
- Five tests that still read `system.file("extdata", "fia_*.csv")` (fit-resume,
  regeneration-adult, tree-loss, reanchor, teacher-forcing) now use the bundled
  `fia_*` datasets; on `main` those CSVs were already gone.
- A bug in the first merge draft (`tree_gen_base` became double after an anchor →
  `scatter_add` dtype error) was fixed before the merge commit.

## Open points for you

1. **Review the merge**, in particular the tree-loss slot remap in `forward()` and
   `build_tree_targets()`, then decide when to merge into `main` (PR).
2. **Naming.** The two per-tree options live on two arguments. Unifying them on
   `loss_aggregation` (e.g. `"tree"` for the observed-tree scoring and
   `"tree_tracked"` for yours) would be easier to explain, but renames your
   argument. Not done; your call.
3. **Per-site anchoring.** Yannek's removed `shooting = "multiple"` anchored each
   site at the start of its own interval and kept the gradient over the whole
   interval; `reanchor` anchors globally by timestep and needs `update_step` to
   divide the interval. Worth porting if inventories with heterogeneous
   intervals use `reanchor`.
4. **`loss_aggregation = "tree"` only runs on the per-site path**
   (`period_mode = "per_site"`); lifting it into the constant-period branch would
   remove a usability trap.
5. **`anchor_cohorts` as a table** (converted with `makeInitCohorts()`), for
   consistency with `tree_data`. Not done.
6. `R CMD check --as-cran` and a pkgdown build have not been run on the merge;
   the precompiled vignettes were not rebuilt.

## Context: what depends on this (FINN-bwi)

- Yannek's BWI paper model (2018–2020 drought from two inventory decades) was
  trained with `6497d01` and uses `loss_aggregation = c(mortality = "tree")`,
  `period_mode = "per_site"`, `quantiles`, `predictTrees()`, a custom optimizer
  passed via `fit(optimizer = )` and `torch_load()`-resumed fits. All of these are
  reproduced bit-for-bit by the merge (table above).
- Two BWI cluster fits compare the two per-tree options on real data
  (`FINN-bwi/fits/2026-11-03-stage2-pertree-{observed,simulated}`); they run
  against this merge installed in a separate library (`~/R/lib-finn-merge`).

## For the Claude session picking this up

- Work on `merge-main-into-loss-aggregation` (or a branch from it); do not push to
  `main` or `loss-aggregation` without Max's go-ahead.
- Before and after any change in `R/finn-class.R`: `devtools::test()` must stay
  at 384+ passing, and the default path must stay bit-identical to `6497d01`
  (and to `d4b3907` with `nonfinite = "skip_step"`). Compare in separate R
  processes with separate libraries; never reinstall into the default library
  while someone's fits run against it.
- The repo's `CLAUDE.md` rules apply (no `git add -A`, roxygen docs, precompiled
  vignettes, CONTRIBUTING via the skills).
