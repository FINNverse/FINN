# ============================================================================
# BCI species-level obs build  (shared prerequisite for Path A and Path B)
# ============================================================================
# Re-aggregates the cleaned per-tree BCI file at SPECIES resolution (vs the
# published 5-PFT / genus aggregations), replicating FINNetAl's
# calculate_stand_vars + period35/25patches packaging. Env is identical across
# resolutions, so it is reused from the published pft path.
#
# Config mirrors the reference pft run: period35, 25 patches.
# Outputs (into dev/pft-bci/data/, nothing written into the FINNetAl repo):
#   obs_species.csv   siteID, year, species, ba, dbh, trees, growth, mort, reg
#   env.csv           copy of the published env_dt (siteID, year, 6 predictors)
#   species_pft.csv   final species code -> Rueger PFT_2axes (ground truth)
# ============================================================================

suppressMessages(library(data.table))
FE  <- "/Users/yannekkaber/working-directory/FINNetAl"
OUT <- "/Users/yannekkaber/working-directory/FINN/dev/pft-bci/data"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dbh_cmTOba_m2 <- function(dbh) { dbh <- dbh/100; pi*dbh^2/4 }
area_of_square_ha <- (1000*500/500)/10000    # 0.1 ha per square (as in FINNetAl)
Npatches <- 25L

tr <- fread(file.path(FE, "data/BCI/data-cleaning/pft/all_trees.csv"))
tr[, species := as.integer(as.factor(sp))]                 # 281 species (was PFT before)
sp_pft <- unique(tr[, .(species, sp, PFT_2axes)])          # keep ground-truth PFT per sp

# ---- stand aggregation per patch x census x species  (calculate_stand_vars) --
sd1 <- tr[, .(
  ba       = sum(dbh_cmTOba_m2(dbh_cm)*nostems*(status=="A"), na.rm = TRUE),
  trees    = sum((status=="A" & !is.na(dbh_cm))*nostems, na.rm = TRUE),
  dbh_mean = sum(dbh_cm*(status=="A" & !is.na(dbh_cm))*nostems, na.rm = TRUE) /
             sum((status=="A" & !is.na(dbh_cm))*nostems, na.rm = TRUE),
  g          = mean(relative_growth, na.rm = TRUE),
  fresh_dead = sum(status=="D" & status_before=="A", na.rm = TRUE),
  r          = sum(census != 1985 & status=="A" & status_before=="P", na.rm = TRUE)/area_of_square_ha
), by = .(census, uniquePatch, species, period_length)]

sd1[census == 1985, c("r","g") := NA]
setorder(sd1, uniquePatch, species, census)
sd1[, trees_before := shift(trees, 1, type = "lag"), by = .(uniquePatch, species)]
sd1[, m := 1 - ((trees_before - fresh_dead)/trees_before)^(1/period_length)]
sd1[is.infinite(m), m := 1]

# ---- rare-species lumping (<=5 growth&mort obs -> "other" = code 0) ----------
Nobs <- sd1[census > 1985, .(growth_mort_N = sum(!is.na(g) & !is.na(m))), by = species]
rare_sp <- Nobs[growth_mort_N <= 5]$species
sd1[, species2 := species][species %in% rare_sp, species2 := 0L]
lvl <- sort(unique(sd1$species2))
sd1[, species := match(species2, lvl)]                      # contiguous 1..N ("other"=1 if present)

# map final code -> PFT (non-other only; "other" is a mix, excluded from matching)
sp_pft[, species2 := ifelse(species %in% rare_sp, 0L, species)]
sp_pft[, code := match(species2, lvl)]
final_pft <- unique(sp_pft[species2 != 0, .(species = code, PFT_2axes)])

# ---- package obs_dt: period35, 25 patches  (as in 01-bci-data-preparation) ---
grid <- fread(file.path(FE, "data/BCI/data-cleaning/site/grid_dt_25patches.csv"))
sd2  <- merge(sd1, grid[, .(siteID, patch, uniquePatch)], by = "uniquePatch")
o1 <- sd2[, .(year = census - min(sd2$census) + 1L, siteID, patch, species,
              ba, dbh = dbh_mean, trees, growth = g, mort = m, reg = r)]
obs <- o1[, .(ba = sum(ba, na.rm=TRUE)/Npatches, dbh = sum(dbh, na.rm=TRUE)/Npatches,
              trees = sum(trees, na.rm=TRUE)/Npatches, growth = sum(growth, na.rm=TRUE)/Npatches,
              mort = sum(mort, na.rm=TRUE)/Npatches, reg = sum(reg, na.rm=TRUE)/Npatches),
          by = .(siteID, year, species)]

# fill absent species x site x year with structural zeros (growth/mort stay NA)
full <- CJ(siteID = sort(unique(obs$siteID)), year = sort(unique(obs$year)),
           species = sort(unique(obs$species)))
obs <- merge(full, obs, by = c("siteID","year","species"), all.x = TRUE)
obs[is.na(ba), ba := 0][is.na(trees), trees := 0][is.na(reg) & year != 1, reg := 0]
obs <- obs[year != 1][, year := year - 1L]                  # drop initial year, reindex
# ANNUAL timestep (period35): censuses are 5 yr apart but FINN steps yearly, so
# each observation spans period_length=5 sub-steps. Required so 0.2.0 sums
# regeneration over the 5 steps (else 5x over-count -> the annual blow-up).
obs[, period_length := 5L]

env <- fread(file.path(FE, "data/BCI/noSplits/pft-period35-25patches/env_dt.csv"))

# ---- initial cohorts (1985 alive trees) at species resolution --------------
# same species remap as obs; dbh binned to 0.1 cm as in create_init_cohorts().
tr[, species2 := ifelse(species %in% rare_sp, 0L, species)]
tr[, sp_code  := match(species2, lvl)]
trg <- merge(tr, grid[, .(siteID, patch, uniquePatch)], by = "uniquePatch")
coh <- trg[census == 1985 & status == "A" & !is.na(dbh_cm),
           .(trees = sum(nostems, na.rm = TRUE)),
           by = .(dbh = as.numeric(as.character(cut(dbh_cm,
                    breaks = seq(1, 351, 0.1),
                    labels = seq(1, 351 - 0.1, 0.1), include.lowest = TRUE))),
                  species = sp_code, siteID, patchID = patch, census)]
setorder(coh, siteID, patchID)
coh[, cohortID := seq_len(.N), by = .(siteID, patchID, census)]

fwrite(obs, file.path(OUT, "obs_species.csv"))
fwrite(env, file.path(OUT, "env.csv"))
fwrite(coh, file.path(OUT, "initial_cohorts1985.csv"))
fwrite(final_pft[order(species)], file.path(OUT, "species_pft.csv"))

# ---- diagnostics -----------------------------------------------------------
present <- obs[trees > 0 | growth != 0]
gm <- obs[!is.na(growth), .N, by = species]
cat("\n==================== species-level BCI build ====================\n")
cat("final species (incl. 'other'):", uniqueN(obs$species),
    "  | lumped into 'other':", length(rare_sp), "of 281\n")
cat("sites:", uniqueN(obs$siteID), " years:", paste(range(obs$year), collapse="-"),
    " obs rows:", nrow(obs), "\n")
cat("growth observations per species: min", min(gm$N), " median", median(gm$N),
    " max", max(gm$N), "\n")
cat("species with < 20 growth obs (the 'rare' regime):",
    nrow(gm[N < 20]), "\n")
cat("Rueger PFT sizes among non-other species:\n"); print(final_pft[, .N, by = PFT_2axes][order(PFT_2axes)])
cat("\nWritten to", OUT, "\n")
