---
name: finn-process
description: How to add or modify a FINN demographic process function (growth, mortality, regeneration, competition, or a custom_parameters variant) in R/Processes.R. Use whenever writing, editing, or debugging a process function, adding parameters to a process, or touching the forward pass in R/finn-class.R that calls one.
---

# Authoring a FINN process function

FINN's demographic processes are plain R functions in `R/Processes.R`, exported
with roxygen `@export`. They are the mechanistic core; `createProcess()` wires one
into a model. Follow these rules — several are non-obvious and have caused bugs.

## Signatures (match the existing ones exactly)

- `growth(dbh, species, parGrowth, pred, light, ...)`
- `mortality(dbh, species, trees, parMort, pred, light, ..., growth = NULL)`
- `regeneration(species, parReg, pred, light, debug = FALSE)`
- `competition(dbh, species, trees, parComp, ..., patch_size_ha, ...)`

A new alternative reuses the base signature of the process it replaces (e.g.
`regeneration_saturation` mirrors `regeneration`). Keep `debug = FALSE` returning
a list of intermediates if the sibling does.

## The function is bound to the model (`self`)

Before use, `setup_species_parameters()` does `environment(func) = self$.__enclos_env__`,
so **inside the function you can read `self$anything`** — buffers, flags
(`self$record_raws`), and any custom parameters. This is how a process reaches
state beyond its formal arguments.

## Two ways to get parameters

1. **Built-in species table** `par_<type>` — sized and bounded by the *type string*
   ("growth"/"mortality"/"regeneration"/"competition") in
   `init_species_parameters()` and `get_par_boundary()` (`R/Processes_utils.R`).
   The forward pass passes a column, e.g. `parReg = self$par_regeneration[,1]`.
   Changing the column count here affects **every** function of that type, so do
   not repurpose it for one variant.
2. **`custom_parameters`** (preferred for anything extra) — declare in
   `createProcess(custom_parameters = list(myparam = init))`. FINN registers it as
   a trainable `nn_parameter` (or buffer if `trainable = FALSE`), `fit()` optimises
   it, `torch_save/torch_load` round-trip it, and the function reads it as
   `self$myparam`. Parameterise in the **unconstrained space** and apply the link
   inside the function (store `log K`, use `exp(self$myparam)`). Give a clear
   `stop()` if a required custom parameter is missing.

   The init **length sets the resolution**: length 1 → one shared value; length
   `N_species` → per-species. This is how `regeneration_saturation` offers shared
   vs per-species `K` without a flag.

## CRITICAL gotcha: species is the LAST tensor dimension

Per-cohort tensors (`r_mean_ha`, `light`, etc.) carry **species in the last
dimension** (see the `sp = N_species` recruit construction in `forward()`, and the
"n species in last dim" contract on `AL_reg`). When broadcasting a per-species
vector against such a tensor, reshape it to `c(1L, 1L, -1L)` — **not**
`c(1L, -1L, 1L)`. Putting it on the wrong axis works by accident for length-1
(shared) values and crashes for per-species ones (this was the
`regeneration_saturation(shared=FALSE)` bug: `11 vs 4` tensor mismatch).

## Design principle

Prefer a **new process function** over a `finn()`-level option. Density dependence,
alternative links, etc. should compose like any other process
(`FINN::regeneration_saturation` is a process, not a `finn()` argument), so they
work with mechanistic *and* hybrid variants and stay self-contained.

## After editing — always

1. `Rscript -e 'roxygen2::roxygenise()'`; stage the regenerated `man/*.Rd` +
   `NAMESPACE` in the same commit (never hand-edit them).
2. Add a **torch-gated** test in `tests/testthat/` (`skip_if_no_torch()` first;
   reuse the tiny FIA subset like `fit_toy_finn()` in `helper.R`). For a new
   parameterised process, test both the shared and per-species init lengths and
   the missing-parameter error — a real `fit()` for a few epochs, asserting it
   runs and the parameter was optimised.
3. Verify with a quick `pkgload::load_all(".")` fit before committing.
