# dev/make_extdata.R  (STEP 1 of the data pipeline; NOT shipped — dev/ is .Rbuildignore'd)
#
# Provenance of the bundled FIA data
# ----------------------------------
# Upstream (the FINN-fia analysis repo, not part of this package):
#   scripts/03_attach_environment.R   attaches bioclim climate to the FIA plots
#   scripts/07_prepare_finn_inputs.R  exports the source CSVs into data-raw/
# This script (make_extdata.R) then subsamples those into the datasets in
# inst/extdata/. Downstream: vignettes/build.R knits the vignettes, which train
# the models inline (there are no .pt pre-fits or .rds caches any more).
#
# Source (data-raw/, ~18 MB, build-ignored) -> inst/extdata/ products:
#   Vignette "C-Data_preparation": a RAW tree list + matching env with
#   siteName/patchName/OrigYear keys, so makeObsData -> resolveSiteIDs ->
#   makeInitCohorts run live on the sample.
#     - example_tree_dt.csv, example_env_dt.csv
#   Vignette "D-Fit_to_FIA": ID-resolved tables, re-indexed (siteID 1..N,
#   species 1..K) so they fit/simulate as-is, split into TRAIN and a disjoint
#   HOLDOUT so the vignette can report out-of-sample performance.
#     - train:   fia_obs_dt.csv, fia_env_dt.csv (RAW climate), fia_init_trees.csv
#     - holdout: fia_obs_test.csv, fia_env_test.csv, fia_init_test.csv
#     - shared:  fia_species_dt.csv  (species coding is derived from TRAIN)
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
## Mortality as a closed-cohort COUNT pair (n_at_risk / n_died).
##
## obs_dt.csv predates makeObsData() emitting these columns, and the upstream
## FINN-fia script that wrote data-raw/ is not re-runnable as-is (it does not
## produce full_tree_dt.csv). So rather than recompute every response - which
## would risk changing which sites are in the sample - we derive ONLY the
## mortality counts, from the same resolved tree table, using the package's own
## makeObsData(), and merge them onto obs_dt. Filters are switched off so the
## site set is untouched; we take just the two count columns.
##
## Once FINN-fia/scripts/07_prepare_finn_inputs.R is re-run against the current
## FINN, obs_dt.csv will carry n_at_risk/n_died natively and this block can go.
mort_src <- full_tree[, .(siteName, patchName, treeName, year, species_name,
                          dbh, status, living)]
mort_counts <- suppressMessages(
  makeObsData(mort_src, plotsize = 0.06, aggregate_by_site = TRUE,
              minNyears = NULL, dbh_growth_thresh = NULL)
)$obs_dt[, .(siteName, year, species_name, n_at_risk, n_died, growth_n)]

# obs_dt is keyed by siteID; full_tree carries both keys, so map across.
site_key <- unique(full_tree[, .(siteID, siteName)])
mort_counts <- merge(mort_counts, site_key, by = "siteName")[, siteName := NULL][]

obs[, mort := NULL]                                   # replaced by the derived rate
obs <- merge(obs, mort_counts, by = c("siteID", "year", "species_name"), all.x = TRUE)
obs[is.na(n_at_risk), n_at_risk := 0L][is.na(n_died), n_died := 0L]
obs[, mort := fifelse(n_at_risk > 0, n_died / n_at_risk, NA_real_)]
stopifnot(all(obs$n_died <= obs$n_at_risk))           # the closed-cohort invariant

## ---------------------------------------------------------------------------
## Vignette "fit-to-fia": a TRAIN sample and a disjoint HOLDOUT sample.
##
## The vignette fits on the train sites and evaluates on sites the model has
## never seen. That is only meaningful because FINN's parameters are per-species
## x environment (not per-site), so a fitted model transfers to new sites.
##
## Sample size matters: growth/mort/reg are per-tree rates and are only weakly
## constrained by two inventories. At 40 sites the rarest species had 2-3 finite
## growth observations, which made the growth Spearman swing by ~0.13 across
## seeds (sd ~0.05) - pure sampling noise. 200 sites gives ~1190 growth
## observations with >=22 per species and no species under 10.
## ---------------------------------------------------------------------------
N_FIT  <- 200   # sites used for fitting
N_TEST <- 200   # disjoint sites used ONLY for evaluation
patches_per_site <- init[, .(np = uniqueN(patchID), nt = .N), by = siteID]
good <- patches_per_site[np >= 2 & nt >= 10, siteID]
stopifnot(length(good) >= N_FIT + N_TEST)

picked     <- sample(good, N_FIT + N_TEST)          # disjoint by construction
fit_sites  <- sort(picked[seq_len(N_FIT)])
test_sites <- sort(picked[N_FIT + seq_len(N_TEST)])
stopifnot(length(intersect(fit_sites, test_sites)) == 0)

## --- species keep-set: decided on TRAIN, applied to BOTH --------------------
## The holdout must use the identical species coding, so the keep-set and the
## species integer map are derived from the training sites only.
K_KEEP <- 10
abund  <- obs[siteID %in% fit_sites][species_name != "other",
              .(trees = sum(trees, na.rm = TRUE)), by = species_name][order(-trees)]
keep   <- head(abund$species_name, K_KEEP)
sp_map <- data.table(species_name = c(keep, "other"), species = seq_along(c(keep, "other")))

# tree-weighted mean that ignores NA / zero-weight rows
wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0
  if (any(ok)) sum(x[ok] * w[ok]) / sum(w[ok]) else NA_real_ }

build_split <- function(sites) {
  smap <- data.table(siteID = sites, new_siteID = seq_along(sites))
  ri <- function(dt) { dt <- merge(dt, smap, by = "siteID")
                       dt[, siteID := new_siteID][, new_siteID := NULL][] }
  o <- ri(obs[siteID %in% sites])
  # ship RAW (untransformed) env: FINN standardizes internally via env_autoscale = TRUE
  e <- ri(env_raw[siteID %in% sites])
  i <- ri(init[siteID %in% sites])

  o[!(species_name %in% keep), species_name := "other"]
  i[!(species_name %in% keep), species_name := "other"]

  # re-aggregate over the lumped species: extensive vars sum, rates tree-weighted.
  # Mortality is the exception: its counts sum exactly, and the rate is derived
  # from the pooled counts afterwards - no weighting scheme to get wrong.
  o <- o[, .(
    ba        = sum(ba,    na.rm = TRUE),
    trees     = sum(trees, na.rm = TRUE),
    reg       = sum(reg,   na.rm = TRUE),
    dbh       = wmean(dbh,    trees),
    # growth means are weighted by the trees BEHIND THE MEAN (growth_n), not by
    # the standing tree count - they are different numbers, and growth_n is the
    # one that sets the mean's precision.
    growth     = wmean(growth, growth_n),
    growth_n   = sum(growth_n, na.rm = TRUE),
    n_at_risk = sum(n_at_risk, na.rm = TRUE),
    n_died    = sum(n_died,    na.rm = TRUE)
  ), by = .(siteID, year, species_name)]
  o[, mort := fifelse(n_at_risk > 0, n_died / n_at_risk, NA_real_)]

  i[, species := NULL]                       # drop stale code before remap
  o <- merge(o, sp_map, by = "species_name")
  i <- merge(i, sp_map, by = "species_name")
  setorder(o, siteID, year, species); setorder(e, siteID, year)
  setorder(i, siteID, patchID, species)
  list(obs = o, env = e, init = i)
}

tr <- build_split(fit_sites)
te <- build_split(test_sites)
obs_f <- tr$obs                                   # reused by the report below
species_f <- copy(sp_map)[, .(species, species_name)][order(species)]

env_out_cols <- c("siteID", "year", env_vars)
# --- train ---
fwrite(tr$obs,                     file.path(out, "fia_obs_dt.csv"))
fwrite(tr$env[, ..env_out_cols],   file.path(out, "fia_env_dt.csv"))   # RAW units
fwrite(tr$init,                    file.path(out, "fia_init_trees.csv"))
# --- holdout (never seen during fitting) ---
fwrite(te$obs,                     file.path(out, "fia_obs_test.csv"))
fwrite(te$env[, ..env_out_cols],   file.path(out, "fia_env_test.csv"))
fwrite(te$init,                    file.path(out, "fia_init_test.csv"))
# --- shared species coding ---
fwrite(species_f,                  file.path(out, "fia_species_dt.csv"))
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
cat(sprintf("\ntrain: %d sites, %d species | holdout: %d sites (disjoint) | prep sample: %d sites\n",
            uniqueN(tr$obs$siteID), uniqueN(tr$obs$species), uniqueN(te$obs$siteID), N_PREP))
gobs <- function(d) sum(is.finite(d$growth))
cat(sprintf("finite growth obs -- train: %d, holdout: %d (was 226 at 40 sites)\n",
            gobs(tr$obs), gobs(te$obs)))
cat(sprintf("min growth obs per species (train): %d\n",
            min(tr$obs[, .(n = sum(is.finite(growth))), by = species]$n)))
