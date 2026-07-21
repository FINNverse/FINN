# Package index

## Model setup

Define a FINN model and how each demographic process is represented.

- [`finn()`](https://finnverse.github.io/FINN/reference/finn.md) :
  Forest Informed Neural Network
- [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  : Define a demographic process for FINN
- [`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md)
  : Define a hybrid (deep-neural-network) demographic process for FINN
- [`CohortMat()`](https://finnverse.github.io/FINN/reference/CohortMat.md)
  : Cohort Matrix Class

## Fitting and prediction

Calibrate every component jointly, end-to-end, then simulate.

- [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) : Fit
  FINN
- [`predict(`*`<finn_class>`*`)`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  : Predict from a FINN model
- [`simulateForest()`](https://finnverse.github.io/FINN/reference/simulateForest.md)
  : Simulate
- [`FINN.seed()`](https://finnverse.github.io/FINN/reference/FINN.seed.md)
  : Set Seed for Reproducibility in R and Torch

## Processes

The built-in mechanistic demographic functions and allometries.

- [`growth()`](https://finnverse.github.io/FINN/reference/growth.md) :
  Calculate growth
- [`mortality()`](https://finnverse.github.io/FINN/reference/mortality.md)
  : Mortality
- [`regeneration()`](https://finnverse.github.io/FINN/reference/regeneration.md)
  : Calculate the regeneration of forest patches based on the input
  parameters.
- [`competition()`](https://finnverse.github.io/FINN/reference/competition.md)
  : Compute the fraction of available light (light) for each cohort
  based on the given parameters.
- [`height()`](https://finnverse.github.io/FINN/reference/height.md) :
  Calculate the height of a tree based on its diameter at breast height
  and an allometry parameter.
- [`BA_stand()`](https://finnverse.github.io/FINN/reference/BA_stand.md)
  : Calculate the Basal Area of a Stand
- [`BA_stem()`](https://finnverse.github.io/FINN/reference/BA_stem.md) :
  Calculate the basal area of a tree given the diameter at breast height
  (dbh).
- [`dbh2ba()`](https://finnverse.github.io/FINN/reference/dbh2ba.md) :
  Convert DBH to basal area

## Data preparation

Turn a raw tree list and environment table into FINN’s inputs.

- [`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md)
  : Create observation data from trees
- [`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
  : Resolve site, patch, and year indices for FINN inputs
- [`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md)
  : Make initial cohorts for FINN
- [`climateDF2array()`](https://finnverse.github.io/FINN/reference/climateDF2array.md)
  : Convert a climate data frame to a FINN environment array
- [`obsDF2arrays()`](https://finnverse.github.io/FINN/reference/obsDF2arrays.md)
  : Convert observation data frame to arrays
- [`array2obsDF()`](https://finnverse.github.io/FINN/reference/array2obsDF.md)
  : Transform Arrays to Observation Data Table
- [`pred2DF()`](https://finnverse.github.io/FINN/reference/pred2DF.md) :
  Convert Prediction Arrays to Data Frames
- [`extract_env()`](https://finnverse.github.io/FINN/reference/extract_env.md)
  : Extract Environmental Data for a Process

## Environmental scaling

FINN standardises environmental predictors internally
(`env_autoscale = TRUE`); these expose the stored constants.

- [`apply_env_scaling()`](https://finnverse.github.io/FINN/reference/apply_env_scaling.md)
  : Apply stored z-standardization to environmental predictors
- [`compute_env_scaling()`](https://finnverse.github.io/FINN/reference/compute_env_scaling.md)
  : Learn z-standardization for environmental predictors

## Interpretation (xAI)

Read back what a fitted - or hybrid - model actually learned.

- [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) :
  Accumulated local effect plots
- [`plot(`*`<FINNale>`*`)`](https://finnverse.github.io/FINN/reference/plot.FINNale.md)
  : Plot ALE curves of a FINN model
- [`summary(`*`<finn_class>`*`)`](https://finnverse.github.io/FINN/reference/summary.finn_class.md)
  : Summarise a fitted FINN model
- [`conditionalEffects()`](https://finnverse.github.io/FINN/reference/conditionalEffects.md)
  : Conditional effects of a FINN model
- [`averageConditionalEffects()`](https://finnverse.github.io/FINN/reference/averageConditionalEffects.md)
  : Average conditional effects of a FINN model
- [`feature_importance()`](https://finnverse.github.io/FINN/reference/feature_importance.md)
  : Permutation feature importance for FINN demographic rates

## Low-level utilities

Building blocks used internally and exported for advanced use.

- [`aggregate_results_old()`](https://finnverse.github.io/FINN/reference/aggregate_results_old.md)
  : Aggregate function
- [`binomial_from_bernoulli()`](https://finnverse.github.io/FINN/reference/binomial_from_bernoulli.md)
  : Draw binomial counts from per-trial Bernoulli probabilities
- [`binomial_from_gamma()`](https://finnverse.github.io/FINN/reference/binomial_from_gamma.md)
  : Sample from binomial with gradient
- [`groupby_mean()`](https://finnverse.github.io/FINN/reference/groupby_mean.md)
  : group by mean
- [`index_species()`](https://finnverse.github.io/FINN/reference/index_species.md)
  : Index species
- [`np_runif()`](https://finnverse.github.io/FINN/reference/np_runif.md)
  : Generate random numbers from a uniform distribution
- [`rweibull_cohorts()`](https://finnverse.github.io/FINN/reference/rweibull_cohorts.md)
  : Generate Cohorts Using Weibull Distribution
- [`sample_poisson_relaxed()`](https://finnverse.github.io/FINN/reference/sample_poisson_relaxed.md)
  : Sample poisson relaxed
