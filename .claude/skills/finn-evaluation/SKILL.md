---
name: finn-evaluation
description: The standard protocol for validating or comparing FINN model fits. Use WHENEVER the user asks to validate, evaluate, or assess a fit, asks "is this good / plausible", or asks "which is better / did X help" comparing two models or configurations. Enforces holdout metrics (R2 + RMSE, species-level AND total) plus a bare-ground equilibrium plausibility check.
---

# Evaluating a FINN fit — do BOTH, always

Never judge a FINN model by its training loss alone. Any validation / "is this
good?" / "which is better?" request triggers **two** mandatory checks. When
comparing models, run the identical protocol on each and compare like for like.

## 1. Holdout accuracy — R² and RMSE, species-level AND total

Predict on **held-out** data (never the training set) and score predicted vs
observed for the state variables `dbh`, `ba`, `trees` (add rates
`growth`/`mort`/`reg` if relevant). Report **per species** and **total**
(stand-level) — both, always.

```r
sim <- predict(model, env = test_env, init_cohort = test_init,
               patches = 4, patch_size = 0.06, device = "cpu")   # match training patch_size
pred_dt <- sim$long$site   # siteID, year, species, variable, value

obs_long <- data.table::melt(
  test_obs, id.vars = c("siteID","year","species"),
  measure.vars = c("dbh","ba","trees"), variable.name = "variable", value.name = "obs")
m <- merge(obs_long, pred_dt, by = c("siteID","year","species","variable"))

metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); o <- obs[ok]; p <- pred[ok]
  list(n = length(o), rmse = sqrt(mean((o - p)^2)),
       r2 = 1 - sum((o - p)^2) / sum((o - mean(o))^2))   # can be < 0 for a bad fit
}

# per species
species_scores <- m[, metrics(obs, value), by = .(species, variable)]
# total: stand-level sum across species per site-year (ba, trees), then score
tot <- m[variable %in% c("ba","trees"),
         .(obs = sum(obs), value = sum(value)), by = .(siteID, year, variable)]
total_scores <- tot[, metrics(obs, value), by = variable]
```

Report both tables. When comparing models, put their species-level and total
R²/RMSE side by side and say which wins **on the holdout**, not on training loss.
A negative R² means the model is worse than predicting the mean — call that out.

## 2. Bare-ground equilibrium — does it converge to something plausible?

Run the model forward from **bare ground** (`init_cohort = NULL`) over a long,
constant (or climatological-mean) environment and check that the stand settles to
a plausible steady state rather than exploding or collapsing.

```r
sitesel <- unique(test_env$siteID)[1]          # one or a few representative sites
base    <- test_env[siteID == sitesel][1]      # a single env row to hold constant
eq_env  <- data.table::CJ(siteID = sitesel, year = 1:800)[
             , `:=`(temp = base$temp, prec = base$prec)]     # use the real env var names
sim_eq  <- predict(model, env = eq_env, init_cohort = NULL,
                   patches = 20, patch_size = 0.06, device = "cpu")
eq <- sim_eq$long$site
```

Then check, on the stand totals (remember **`ba`/`trees` are per `patch_size`** —
divide by `patch_size` for per-ha before judging plausibility):

- **Convergence:** the mean of the last ~100 years vs the preceding ~100 years
  changes by only a few percent (state has stabilised, not still drifting).
- **Plausibility (per ha):** basal area in a sensible forest range (roughly
  10–80 m² ha⁻¹, ecosystem-dependent), stem density positive and finite, a
  non-degenerate species composition. Flag anything that **blows up** (BA →
  huge / `Inf`/`NaN`) or **collapses** to ~0.
- A bare-ground run is exactly where **unbounded regeneration explodes** — if it
  does, that is a real defect, and density-dependent recruitment
  (`FINN::regeneration_saturation`) is the intended fix, not a tuning hack.

## Verdict

Finish with an explicit verdict: holdout R²/RMSE (species + total) **and** whether
the equilibrium is convergent and plausible. A model that fits the holdout well but
diverges from bare ground is **not** validated — say so. For "which is better?",
the winner must be better (or at least not worse) on the holdout metrics *and*
remain equilibrium-plausible; if one model trades holdout accuracy for a stable
equilibrium, present that trade-off rather than declaring a single winner.
