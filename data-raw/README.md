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

## The pipeline (each step is a script)

```
data-raw/*.csv
   │  dev/make_extdata.R        ← STEP 1: subsample + re-index
   ▼
inst/extdata/  fia_obs_dt.csv, fia_env_dt.csv (RAW climate),
               fia_init_trees.csv, fia_species_dt.csv,
               example_tree_dt.csv, example_env_dt.csv
   │
   ├─ dev/train_fia_model.R      ← STEP 2: fit → fia_process_finn.pt, fia_hybrid_finn.pt
   │
   └─ dev/precompute_vignettes.R ← STEP 3: cache vignette results →
                                    vig_fit_to_fia.rds, vig_succession.rds, vig_intro.rds
```

Run them in order from the package root; steps 2 and 3 need a torch backend.
`make_extdata.R` is seeded (`set.seed(42)`), so re-running reproduces the same
sample as long as these source files are unchanged.

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
