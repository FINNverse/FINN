# `data-raw/` — source data behind the bundled FIA datasets

This directory holds the **source** Oregon FIA + climate CSVs that the bundled
datasets in `inst/extdata/` are derived from. It is **`.Rbuildignore`d**, so none
of it ships to CRAN — it exists only to make the shipped data reproducible.

## Where this comes from (upstream)

The CSVs here are exported by the **FINN-fia analysis repo** (a separate project,
not this package):

| FINN-fia script | produces |
|---|---|
| `scripts/03_attach_environment.R` | attaches the bioclim **climate** variables (`temp`, `tempmax`, `tempmin`, `prec`, `precseas`, `precwarmq`) to each FIA plot |
| `scripts/07_prepare_finn_inputs.R` | writes the source CSVs copied into this `data-raw/` |

So the **climate data** ultimately originates from `03_attach_environment.R`; the
raw (untransformed) values live in `env_unscaled_dt.csv`.

## The pipeline

```
data-raw/*.csv
   │  dev/make_extdata.R        ← STEP 1: subsample + re-index
   ▼
inst/extdata/  fia_obs_dt.csv, fia_env_dt.csv (RAW climate),
               fia_init_trees.csv, fia_species_dt.csv,
               example_tree_dt.csv, example_env_dt.csv
   │
   └─ vignettes/build.R         ← STEP 2: knit *.Rmd.orig -> *.Rmd
                                   (trains the models, bakes in results + figures)
```

Run from the package root; STEP 2 needs a torch backend. `make_extdata.R` is
seeded (`set.seed(42)`), so re-running reproduces the same sample as long as
these source files are unchanged.

## How the vignettes are built (precompiled)

The vignettes train torch models, which CRAN can neither run (no libtorch
backend) nor afford (time). So we use the **precompile** pattern:

| file | role |
|---|---|
| `vignettes/<name>.Rmd.orig` | the **source** — real, live code |
| `vignettes/build.R` | knits `.Rmd.orig` → `.Rmd`, running everything once |
| `vignettes/<name>.Rmd` | **generated & committed** — static, results baked in |
| `vignettes/<name>/` | **generated & committed** — the figure PNGs |

`R CMD build` then only runs pandoc over static markdown, so the package builds
anywhere in seconds without torch.

```sh
# from the package root, with a torch backend available
Rscript vignettes/build.R                 # knit all (fit-to-fia trains 2 models, ~10 min)
Rscript vignettes/build.R succession-demo # knit just one
```

Then **commit the regenerated `.Rmd` and its figure folder** — those are what
ship.

Why this matters:

* **Nothing heavy ships.** The models are trained inside the vignette at knit
  time and never leave the developer machine. There are no `.pt` files and no
  result caches in the package, so vignette size no longer limits what we can
  demonstrate.
* **No drift.** The code shown in a vignette is exactly the code that produced
  the output next to it, because one source generated both. (The previous scheme
  — `eval: false` code plus a separate `.rds` cache — could silently disagree,
  and did.)
* **No quarto dependency.** The vignettes are `knitr::rmarkdown`
  (`VignetteBuilder: knitr`), so no quarto binary is needed to build the package.

Never edit a `.Rmd` by hand — it is overwritten on the next `build.R` run. Edit
the `.Rmd.orig`.

## Files here

**Used by `dev/make_extdata.R`:**

| file | role |
|---|---|
| `full_tree_dt.csv` | raw tree list (for the data-preparation example + site/year keys) |
| `obs_dt.csv` | site-level observations (years 1–2) |
| `env_unscaled_dt.csv` | **RAW climate** per site × year |
| `init_trees.csv` | year-0 tree slice (initial cohorts) |
| `species_dt.csv` | species code ↔ name lookup |

**Not used (kept for reference only):**

| file | why unused |
|---|---|
| `env_dt.csv` | z-scaled climate — the model now standardizes internally (`env_autoscale = TRUE`) |
| `env_scales_dt.csv` | the matching scaling constants — likewise superseded by `m$env_scaling` |
| `full_obs_dt.csv` | superset of `obs_dt.csv`; not read by the pipeline |
