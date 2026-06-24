# dev/make_extdata.R  (STEP 1 of the data pipeline; NOT shipped — dev/ is .Rbuildignore'd)
#
# Provenance of the bundled FIA data
# ----------------------------------
# Upstream (the FINN-fia analysis repo, not part of this package):
#   scripts/03_attach_environment.R   attaches bioclim climate to the FIA plots
#   scripts/07_prepare_finn_inputs.R  exports the source CSVs into data-raw/
# This script (make_extdata.R) then subsamples those into the tiny datasets in
# inst/extdata/. Downstream: dev/train_fia_model.R turns them into the .pt
# pre-fits, and dev/precompute_vignettes.R into the vig_*.rds vignette caches.
#
# Source (data-raw/, ~18 MB, build-ignored) -> inst/extdata/ products:
#   Vignette "data-preparation": a RAW tree list + matching env with
#   siteName/patchName/OrigYear keys, so makeObsData -> resolveSiteIDs ->
#   makeInitCohorts run live on the sample.
#     - example_tree_dt.csv, example_env_dt.csv
#   Vignette "fit-to-fia": ID-resolved tables, subsampled to ~40 sites and
#   re-indexed (siteID 1..N, species 1..K) so they fit/simulate as-is.
#     - fia_obs_dt.csv, fia_env_dt.csv (RAW climate), fia_init_trees.csv,
#       fia_species_dt.csv
#
# Run from the package root:  Rscript dev/make_extdata.R
suppressMessages({library(data.table); library(FINN)})
set.seed(42)

src <- "data-raw"          # source CSVs from FINN-fia scripts/07_prepare_finn_inputs.R
out <- "inst/extdata"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Only the five source files actually used are read here. (data-raw/ also holds
# env_dt.csv (z-scaled) + env_scales_dt.csv — now unused because the model stores
# its own standardization via env_autoscale = TRUE — and full_obs_dt.csv, unused.)
full_tree <- fread(file.path(src, "full_tree_dt.csv"))
obs       <- fread(file.path(src, "obs_dt.csv"))             # site-level, years 1..2
env_raw   <- fread(file.path(src, "env_unscaled_dt.csv"))    # RAW climate, years 0..2
init      <- fread(file.path(src, "init_trees.csv"))         # year 0
species_lkp <- fread(file.path(src, "species_dt.csv"))

env_vars <- c("temp", "tempmax", "tempmin", "prec", "precseas", "precwarmq")

## ---------------------------------------------------------------------------
## Vignette 3 sample: pick N_FIT sites with >=2 patches and decent tree counts.
## ---------------------------------------------------------------------------
N_FIT <- 40
patches_per_site <- init[, .(np = uniqueN(patchID), nt = .N), by = siteID]
good <- patches_per_site[np >= 2 & nt >= 10, siteID]
fit_sites <- sort(sample(good, min(N_FIT, length(good))))

# site re-index map (old siteID -> 1..N), applied to every fit table
site_map <- data.table(siteID = fit_sites, new_siteID = seq_along(fit_sites))

reindex_sites <- function(dt) {
  dt <- merge(dt, site_map, by = "siteID")
  dt[, siteID := new_siteID][, new_siteID := NULL]
  dt[]
}

obs_f  <- reindex_sites(obs[siteID %in% fit_sites])
# ship RAW (untransformed) env: FINN standardizes internally via env_autoscale = TRUE
env_f  <- reindex_sites(env_raw[siteID %in% fit_sites])
init_f <- reindex_sites(init[siteID %in% fit_sites])

## --- lump to the K most abundant species; rest -> "other" -------------------
## 28 species over 40 sites is badly under-constrained for a demo fit, so keep
## the dominant species and merge the long tail (mirrors makeObsData's "other").
K_KEEP <- 10
abund  <- obs_f[, .(trees = sum(trees, na.rm = TRUE)), by = species_name][order(-trees)]
keep   <- head(abund$species_name[abund$species_name != "other"], K_KEEP)
obs_f[!(species_name %in% keep),  species_name := "other"]
init_f[!(species_name %in% keep), species_name := "other"]

# tree-weighted mean that ignores NA / zero-weight rows
wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0
  if (any(ok)) sum(x[ok] * w[ok]) / sum(w[ok]) else NA_real_ }

# re-aggregate obs over the lumped species: extensive vars sum, rates tree-weighted
obs_f <- obs_f[, .(
  ba     = sum(ba,    na.rm = TRUE),
  trees  = sum(trees, na.rm = TRUE),
  reg    = sum(reg,   na.rm = TRUE),
  dbh    = wmean(dbh,    trees),
  growth = wmean(growth, trees),
  mort   = wmean(mort,   trees)
), by = .(siteID, year, species_name)]

# species re-index 1..K ordered by abundance, "other" last
sp_order  <- c(keep, "other")
sp_map    <- data.table(species_name = sp_order, species = seq_along(sp_order))
init_f[, species := NULL]                       # drop stale code before remap
obs_f  <- merge(obs_f,  sp_map, by = "species_name")
init_f <- merge(init_f, sp_map, by = "species_name")
species_f <- copy(sp_map)[, .(species, species_name)][order(species)]

setorder(obs_f, siteID, year, species)
setorder(env_f, siteID, year)
setorder(init_f, siteID, patchID, species)

env_out_cols <- c("siteID", "year", env_vars)
fwrite(obs_f,            file.path(out, "fia_obs_dt.csv"))
fwrite(env_f[, ..env_out_cols], file.path(out, "fia_env_dt.csv"))  # RAW units
fwrite(init_f,           file.path(out, "fia_init_trees.csv"))
fwrite(species_f,        file.path(out, "fia_species_dt.csv"))
# note: no fia_env_scales_dt.csv — the model now stores the standardization
# constants itself (m$env_scaling) when fit with env_autoscale = TRUE.

## ---------------------------------------------------------------------------
## Vignette 2 sample: RAW tree list + env for a few sites, keyed by siteName.
## Use OrigYear as `year` so period_length = 10 and fix_period_length works.
## ---------------------------------------------------------------------------
N_PREP <- 8
prep_sites <- sort(sample(fit_sites, N_PREP))      # reuse fit sites (any work)

raw_cols <- c("siteName","patchName","treeName","species_name","dbh",
              "status","status_before","mort_cause","reg","mort","living","complete")
example_tree <- full_tree[siteID %in% prep_sites]
example_tree[, year := OrigYear]
tree_out_cols <- c("year", raw_cols)
example_tree <- example_tree[, ..tree_out_cols]
setcolorder(example_tree, c("siteName","patchName","treeName","year","species_name","dbh","status"))
setorder(example_tree, siteName, patchName, treeName, year)
fwrite(example_tree, file.path(out, "example_tree_dt.csv"))

# matching env keyed by siteName + OrigYear (back from the resolved env_raw)
site_year_map <- unique(full_tree[siteID %in% prep_sites, .(siteID, siteName, year, OrigYear)])
example_env <- merge(env_raw[siteID %in% prep_sites], site_year_map,
                     by = c("siteID","year"))
env_prep_cols <- c("siteName", "OrigYear", env_vars)
example_env <- example_env[, ..env_prep_cols]
setnames(example_env, "OrigYear", "year")
# a site's patches can be measured across two calendar years, so the
# resolved-year -> OrigYear map is one-to-many; the env values are identical
# across those duplicates -> collapse to one row per (siteName, year).
example_env <- unique(example_env)
setorder(example_env, siteName, year)
fwrite(example_env, file.path(out, "example_env_dt.csv"))

## ---------------------------------------------------------------------------
## Report
## ---------------------------------------------------------------------------
files <- list.files(out, full.names = TRUE)
cat("\n=== inst/extdata written ===\n")
info <- file.info(files)
for (f in files) cat(sprintf("  %-26s %6.1f KB\n", basename(f), info[f, "size"]/1024))
cat(sprintf("\nfit sample: %d sites, %d species | prep sample: %d sites\n",
            uniqueN(obs_f$siteID), uniqueN(obs_f$species), N_PREP))
