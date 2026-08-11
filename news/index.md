# Changelog

## FINN 0.1.0

CRAN release: 2026-08-09

First public release.

FINN is a differentiable forest gap model: a cohort-based dynamic
vegetation model whose demographic processes (competition, growth,
mortality, regeneration) can each be a mechanistic function, a neural
network, or a mixture of the two, all calibrated end-to-end by gradient
descent.

Highlights of the current interface:

- [`finn()`](https://finnverse.github.io/FINN/reference/finn.md)
  assembles a model from one process per demographic component, each
  built with
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  (mechanistic) or
  [`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md)
  (neural network).
- [`simulateForest()`](https://finnverse.github.io/FINN/reference/simulateForest.md)
  runs a model forward;
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)
  calibrates one to data and
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  scores it, returning patch- and site-level results.
- Per-response likelihoods: Gaussian (`"mse"`/`"gaussian"`), Poisson,
  negative binomial, and a `"binomial"` likelihood for mortality that
  takes a closed-cohort count pair (`n_at_risk`, `n_died`).
- `weights = "auto"` (the default in
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)) scales
  each loss by its intercept-only baseline, so the six responses are
  commensurable and every term reads as a fraction of its own null
  deviance.
- Data helpers
  [`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md),
  [`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
  and
  [`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md)
  turn a raw tree list into FINN’s input tables.
- Model interpretation with
  [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md),
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html),
  [`feature_importance()`](https://finnverse.github.io/FINN/reference/feature_importance.md)
  and
  [`conditionalEffects()`](https://finnverse.github.io/FINN/reference/conditionalEffects.md).

See the vignettes for a guided tour: *Introduction to FINN*, *Plausible
succession from a handful of species*, *Preparing your data for FINN*,
*Fitting FINN to forest inventory data*, and *Mortality*.
