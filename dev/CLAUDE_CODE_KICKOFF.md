# Claude Code kickoff

How to run this build in Claude Code, and a paste-ready first prompt.

## Setup
1. Open Claude Code in the `FINN-main/` package root (so it reads `CLAUDE.md`).
2. Make sure these are installed: the package itself (`devtools::load_all()` /
   `install`), `torch` with backend (`torch::install_torch()`), `quarto` R
   package + Quarto CLI ≥ 1.4.
3. Read order for the agent: `CLAUDE.md` → `dev/vignette-plan.md` → the three
   `vignettes/*.qmd` stubs.

## Suggested working order
Do them in this sequence — each unblocks the next:
1. **succession-demo** (no data dependency; validates the simulate API end to end).
2. **make_extdata** (`dev/make_extdata.R`): subsample `vignettes/fit-fia-data/`
   into small `inst/extdata/` datasets + `inst/extdata/example_tree_dt.csv`.
3. **data-preparation** (consumes the tiny raw sample).
4. **train_fia_model** (`dev/train_fia_model.R`): fit once, save
   `inst/extdata/fia_process_finn.pt` (or precomputed predictions).
5. **fit-to-fia** (loads the precomputed model).
6. Wire `_pkgdown.yml`, run the verification loop, `R CMD check`.

## Paste-ready first prompt

> You are completing three Quarto vignettes for the FINN R package. Read
> `CLAUDE.md` and `dev/vignette-plan.md` first — they define the API, the data
> contracts, and hard constraints (package vignettes: self-contained, fast, no
> external downloads, no API changes).
>
> Start with `vignettes/succession-demo.qmd`. Fill in the species-parameter
> matrices to encode the trade-off table in the plan (pioneer → climax), run the
> simulation, and produce the headline succession figure. Verify by rendering
> with `quarto render vignettes/succession-demo.qmd` and confirming the figure
> shows pioneers peaking early and shade-tolerant species dominating late. Tune
> parameters until the dynamics are ecologically sensible — do not fake the plot.
> Then stop and show me the rendered result before moving to the next vignette.
>
> After I approve, create `dev/make_extdata.R` to subsample
> `vignettes/fit-fia-data/` into `inst/extdata/` (each file < ~300 KB; ~30–50
> sites), then complete `data-preparation.qmd` and `fit-to-fia.qmd` per the plan.

## Per-vignette acceptance gates (must pass before "done")
- `succession-demo`: renders < ~30 s; figure shows clear successional turnover;
  prose explains it via the parameters.
- `data-preparation`: `makeObsData`/`resolveSiteIDs`/`makeInitCohorts` all run on
  the shipped sample and print output; acquisition steps `eval=FALSE`; column
  tables match `R/createInputs.R`.
- `fit-to-fia`: renders without training (loads precomputed objects); shows ≥1
  fit metric and ≥1 niche/parameter figure; all inputs via `system.file()`.
- Global: `R CMD build .` + `R CMD check --no-manual` clean with vignettes built;
  no absolute paths; nothing read outside the package at knit time.

## Notes / risks
- `torch` model objects can be fragile across versions. If shipping
  `fia_process_finn.pt` causes load errors, fall back to shipping a precomputed
  predictions/parameters data.frame and adjust `fit-to-fia.qmd` accordingly
  (documented as an option in the plan).
- Keep dev scripts in `dev/` (already `.Rbuildignore`d) so they never ship.
