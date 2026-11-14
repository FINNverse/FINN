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
# Source (data-raw/, ~18 MB, build-ignored) -> inst/extdata/ + data/ products:
#   Vignette "C-Data_preparation": a RAW tree list + matching env with
#   siteName/patchName/OrigYear keys, so makeObsData -> resolveSiteIDs ->
#   makeInitCohorts run live on the sample.
#     - example_tree_dt.csv, example_env_dt.csv
#   Vignettes "D-Fit_to_FIA" / "E-Mortality" (and the torch-gated tests):
#   ID-resolved tables, re-indexed (siteID 1..N, species 1..K) so they
#   fit/simulate as-is, split into TRAIN and a disjoint HOLDOUT so the vignettes
#   can report out-of-sample performance. These ship as package DATASETS
#   (data/*.rda, documented in R/data.R), not as CSVs:
#     - train:   fia_obs_dt, fia_env_dt (RAW climate), fia_init_trees
#     - holdout: fia_obs_test, fia_env_test, fia_init_test
#     - shared:  fia_species_dt  (species coding is derived from TRAIN)
#   All on the ANNUAL time axis of FINN 0.2.0 (one simulation step = one year):
#   `year` = years since the site's first inventory (10 and 20, for two sites
#   21 or 22), `period_length` = years since the previous inventory, growth and
#   mortality as yearly rates, the counts and `reg` over the interval, env rows
#   1..Tmax. Before 0.2.0 the tables used one step per ten-year interval
#   (year 1, 2; env 0..2) with the rates over the interval.
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
## The ANNUAL time axis (FINN 0.2.0): one simulation step is one year.
##
## The source tables index the three inventories as year 0, 1, 2 (OrigYear is
## the calendar year, 2001-02, 2011-12, 2021-24). The datasets put every
## inventory at its CALENDAR OFFSET from the site's first one (10 and 20 for
## nearly every site; 21 or 22 where the third inventory came late), with
## `period_length` = the years since the previous inventory. env gets one row
## per site and year 1..Tmax: env row k drives the year that ends k years after
## the first inventory, so the inventory at offset k is scored against row k.
cal <- unique(full_tree[, .(siteID, year, OrigYear)])
stopifnot(cal[, .N, by = .(siteID, year)][, all(N == 1)])   # one calendar year per site and inventory
setorder(cal, siteID, year)
cal[, `:=`(offset = OrigYear - OrigYear[1L],
           period_length = OrigYear - shift(OrigYear)), by = siteID]

## ---------------------------------------------------------------------------
## Growth and mortality as ANNUAL rates, mortality as a closed-cohort COUNT pair.
##
## obs_dt.csv predates makeObsData() emitting the count columns and annual
## rates, and the upstream FINN-fia script that wrote data-raw/ is not
## re-runnable as-is (it does not produce full_tree_dt.csv). So rather than
## recompute every response - which would risk changing which sites are in the
## sample - we derive ONLY growth, growth_n and the mortality counts, from the
## same resolved tree table with the package's own makeObsData() on CALENDAR
## years (so that period_length is in years and the rates come out per year),
## and merge them onto obs_dt. Filters are switched off so the site set is
## untouched. ba, dbh, trees and reg (recruits per ha over the interval) stay
## those of obs_dt.csv.
rate_src <- full_tree[, .(siteName, patchName, treeName, year = OrigYear, species_name,
                          dbh, status, living)]
rate_cols <- c("siteName", "year", "species_name", "growth", "growth_n", "n_at_risk", "n_died")
rates_annual <- suppressMessages(
  makeObsData(rate_src, plotsize = 0.06, aggregate_by_site = TRUE,
              minNyears = NULL, dbh_growth_thresh = NULL, rates = "annual")
)$obs_dt[, ..rate_cols]
# The growth of obs_dt.csv is the UNWEIGHTED mean of the patch means of the
# interval increment (the site aggregation before FINN 0.2.0), while makeObsData()
# now weights the patches by their trees (growth_n), which is what the model
# predicts and what the shipped growth_n already assumed. Checked here, so the
# replacement below is known to change exactly two things: the unit (per year)
# and the patch weighting.
rates_interval <- suppressMessages(
  makeObsData(rate_src, plotsize = 0.06, aggregate_by_site = FALSE,
              minNyears = NULL, dbh_growth_thresh = NULL, rates = "interval")
)$obs_dt[, .(growth_interval = mean(growth, na.rm = TRUE)), by = .(siteName, year, species_name)]

# obs_dt is keyed by siteID and inventory index; full_tree carries siteName and
# the calendar year, so map across.
site_key <- unique(full_tree[, .(siteID, siteName)])
rates_annual <- merge(rates_annual, rates_interval, by = c("siteName", "year", "species_name"))
rates_annual <- merge(rates_annual, site_key, by = "siteName")[, siteName := NULL][]
rates_annual <- merge(rates_annual, cal[, .(siteID, OrigYear, census = year)],
                      by.x = c("siteID", "year"), by.y = c("siteID", "OrigYear"))
rates_annual[, year := census][, census := NULL]

# the unweighted patch mean of makeObsData()'s interval growth is the growth of obs_dt.csv
chk <- merge(obs[is.finite(growth), .(siteID, year, species_name, growth)],
             rates_annual[, .(siteID, year, species_name, growth_interval)],
             by = c("siteID", "year", "species_name"))
stopifnot(nrow(chk) == obs[is.finite(growth), .N],
          isTRUE(all.equal(chk$growth, chk$growth_interval, tolerance = 1e-6)))
rates_annual[, growth_interval := NULL]

obs[, c("mort", "growth") := NULL]                    # replaced by the annual rates
obs <- merge(obs, rates_annual, by = c("siteID", "year", "species_name"), all.x = TRUE)
obs[is.na(n_at_risk), n_at_risk := 0L][is.na(n_died), n_died := 0L][is.na(growth_n), growth_n := 0L]
stopifnot(all(obs$n_died <= obs$n_at_risk))           # the closed-cohort invariant
obs <- merge(obs, cal[, .(siteID, year, offset, period_length)], by = c("siteID", "year"))
obs[, year := offset][, offset := NULL]

# the yearly mortality rate of a closed cohort: the constant annual probability
# that compounds to n_died / n_at_risk over period_length years
annual_mort <- function(n_died, n_at_risk, p)
  fifelse(n_at_risk > 0, 1 - (1 - n_died / n_at_risk)^(1 / p), NA_real_)

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

# the longest first-to-last-inventory span of the sample: env runs 1..TMAX for
# every site of both splits, so train and holdout share one time axis
TMAX <- cal[siteID %in% picked, max(offset)]

build_split <- function(sites) {
  smap <- data.table(siteID = sites, new_siteID = seq_along(sites))
  ri <- function(dt) { dt <- merge(dt, smap, by = "siteID")
                       dt[, siteID := new_siteID][, new_siteID := NULL][] }
  o <- ri(obs[siteID %in% sites])
  # ship RAW (untransformed) env: FINN standardizes internally via env_autoscale = TRUE.
  # The WorldClim normals do not change between the source's years, so one row
  # per site is repeated for every simulated year 1..TMAX.
  e_site <- unique(env_raw[siteID %in% sites, c("siteID", env_vars), with = FALSE])
  stopifnot(!anyDuplicated(e_site$siteID))
  e <- ri(e_site[, .(year = seq_len(TMAX)), by = c("siteID", env_vars)])
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
  ), by = .(siteID, year, period_length, species_name)]
  o[, mort := annual_mort(n_died, n_at_risk, period_length)]
  setcolorder(o, setdiff(names(o), "period_length"))   # last column

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
# Shipped as lazy-loaded package datasets (data/*.rda; docs in R/data.R). One
# object per file, named exactly like the object, xz-compressed: the seven
# tables are ~120 KB this way versus ~2.2 MB as CSVs. After re-running this,
# R CMD INSTALL . before knitting the vignettes (vignettes/build.R checks).
dir.create("data", showWarnings = FALSE)
save_dataset <- function(x, name) {
  assign(name, x)
  save(list = name, file = file.path("data", paste0(name, ".rda")), compress = "xz")
}
# --- train ---
save_dataset(tr$obs,                   "fia_obs_dt")
save_dataset(unique(tr$env[, ..env_out_cols]), "fia_env_dt")   # RAW units; one row per site x year (the source has one per patch, identical)
save_dataset(tr$init,                  "fia_init_trees")
# --- holdout (never seen during fitting) ---
save_dataset(te$obs,                   "fia_obs_test")
save_dataset(unique(te$env[, ..env_out_cols]), "fia_env_test")
save_dataset(te$init,                  "fia_init_test")
# --- shared species coding ---
save_dataset(species_f,                "fia_species_dt")
# note: no env-scales table — the model stores the standardization constants
# itself (m$env_scaling) when fit with env_autoscale = TRUE.

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

# matching env keyed by siteName + CALENDAR year, one row per year: the annual
# time axis of resolveSiteIDs() needs every year from the year after a site's
# first inventory to the longest first-to-last span of the sample after it.
# The climate is the static WorldClim normal of the plot, repeated; the first
# inventory year is included too, so time_axis = "census" (which keeps the
# inventory years) also finds its rows.
prep_cal  <- cal[siteID %in% prep_sites]
span_prep <- prep_cal[, max(offset)]
prep_first <- merge(prep_cal[year == 0, .(siteID, first = OrigYear)],
                    unique(full_tree[siteID %in% prep_sites, .(siteID, siteName)]), by = "siteID")
env_prep_site <- unique(env_raw[siteID %in% prep_sites, c("siteID", env_vars), with = FALSE])
stopifnot(!anyDuplicated(env_prep_site$siteID))
example_env <- merge(prep_first, env_prep_site, by = "siteID")[
  , .(year = first + 0:span_prep), by = c("siteName", env_vars)]
env_prep_cols <- c("siteName", "year", env_vars)
example_env <- example_env[, ..env_prep_cols]
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
iv <- function(d) unique(d[, .(siteID, year, period_length)])[, .N, by = .(year, period_length)][order(year, period_length)]
cat(sprintf("annual axis: env years 1..%d; inventories (year: period_length x sites) -- train %s | holdout %s\n", TMAX,
            iv(tr$obs)[, paste(sprintf("%d: %d x %d", year, period_length, N), collapse = ", ")],
            iv(te$obs)[, paste(sprintf("%d: %d x %d", year, period_length, N), collapse = ", ")]))
cat(sprintf("annual rates (train) -- growth %.4f..%.4f, mort %.4f..%.4f\n",
            min(tr$obs$growth, na.rm = TRUE), max(tr$obs$growth, na.rm = TRUE),
            min(tr$obs$mort, na.rm = TRUE), max(tr$obs$mort, na.rm = TRUE)))
