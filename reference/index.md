# Package index

## Model setup

Define a FINN model and how each demographic process is represented.

- [`finn()`](https://finnverse.github.io/FINN/reference/finn.md) : Build
  a FINN forest model
- [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  : Define a demographic process for FINN
- [`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md)
  : Define a hybrid (deep-neural-network) demographic process for FINN
- [`CohortMat()`](https://finnverse.github.io/FINN/reference/CohortMat.md)
  : Cohort Matrix Class
- [`recruitment`](https://finnverse.github.io/FINN/reference/recruitment.md)
  [`FINN-recruitment`](https://finnverse.github.io/FINN/reference/recruitment.md)
  : Recruitment in FINN: the options and how to choose them

## Fitting and prediction

Calibrate every component jointly, end-to-end, then simulate.

- [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) :
  Estimate the parameters of a FINN model
- [`finn_data`](https://finnverse.github.io/FINN/reference/finn_data.md)
  [`finn-data`](https://finnverse.github.io/FINN/reference/finn_data.md)
  : The input tables of FINN: observations, environment, initial trees
- [`finn_loss`](https://finnverse.github.io/FINN/reference/finn_loss.md)
  [`finn-loss`](https://finnverse.github.io/FINN/reference/finn_loss.md)
  : How fit() compares the simulation with the observations: the loss
- [`finnControl()`](https://finnverse.github.io/FINN/reference/finnControl.md)
  [`print(`*`<finnControl>`*`)`](https://finnverse.github.io/FINN/reference/finnControl.md)
  : Settings of the optimizer used by fit()
- [`predict(`*`<finn_class>`*`)`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  [`simulateForest()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  : Simulate a forest with a FINN model
- [`plot(`*`<finn_simulation>`*`)`](https://finnverse.github.io/FINN/reference/plot.finn_simulation.md)
  [`print(`*`<finn_simulation>`*`)`](https://finnverse.github.io/FINN/reference/plot.finn_simulation.md)
  : Plot a FINN simulation
- [`predictTrees()`](https://finnverse.github.io/FINN/reference/predictTrees.md)
  : Predict growth and mortality for individual observed trees
- [`FINN.seed()`](https://finnverse.github.io/FINN/reference/FINN.seed.md)
  : Set Seed for Reproducibility in R and Torch

## Evaluation

Judge a fitted model against the observations, on the training data and
a holdout.

- [`evaluate()`](https://finnverse.github.io/FINN/reference/evaluate.md)
  [`print(`*`<finn_evaluation>`*`)`](https://finnverse.github.io/FINN/reference/evaluate.md)
  [`plot(`*`<finn_evaluation>`*`)`](https://finnverse.github.io/FINN/reference/evaluate.md)
  : Evaluate a fitted model
- [`comparePlot()`](https://finnverse.github.io/FINN/reference/comparePlot.md)
  [`print(`*`<finn_heatmap>`*`)`](https://finnverse.github.io/FINN/reference/comparePlot.md)
  [`print(`*`<finn_comparison>`*`)`](https://finnverse.github.io/FINN/reference/comparePlot.md)
  : Compare a simulation with observations
- [`summary(`*`<finn_class>`*`)`](https://finnverse.github.io/FINN/reference/summary.finn_class.md)
  [`print(`*`<summary.finn_class>`*`)`](https://finnverse.github.io/FINN/reference/summary.finn_class.md)
  : Summarise a fitted FINN model

## Processes

The built-in mechanistic demographic functions and allometries.

- [`growth()`](https://finnverse.github.io/FINN/reference/growth.md) :
  Growth
- [`mortality()`](https://finnverse.github.io/FINN/reference/mortality.md)
  : Mortality
- [`regeneration()`](https://finnverse.github.io/FINN/reference/regeneration.md)
  : Regeneration without a recruitment cap
- [`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md)
  [`regeneration_saturation()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md)
  : Regeneration with a Beverton-Holt recruitment cap
- [`regeneration_adult()`](https://finnverse.github.io/FINN/reference/regeneration_adult.md)
  : Regeneration limited by conspecific adult density
- [`regeneration_adult_effect()`](https://finnverse.github.io/FINN/reference/regeneration_adult_effect.md)
  : Regeneration with an effect of conspecific adult density
- [`competition()`](https://finnverse.github.io/FINN/reference/competition.md)
  : Compute the fraction of available light (light) for each cohort
  based on the given parameters
- [`height()`](https://finnverse.github.io/FINN/reference/height.md) :
  Calculate the height of a tree based on its diameter at breast height
  and an allometry parameter
- [`BA_stand()`](https://finnverse.github.io/FINN/reference/BA_stand.md)
  : Calculate the Basal Area of a Stand
- [`BA_stem()`](https://finnverse.github.io/FINN/reference/BA_stem.md) :
  Calculate the basal area of a tree given the diameter at breast height
  (dbh)
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
  : Convert a cohort table to arrays
- [`array2obsDF()`](https://finnverse.github.io/FINN/reference/array2obsDF.md)
  : Transform Arrays to Observation Data Table
- [`pred2DF()`](https://finnverse.github.io/FINN/reference/pred2DF.md) :
  Convert Prediction Arrays to Data Frames

## Data

The bundled Oregon FIA calibration data, in FINN’s input format.

- [`fia_obs_dt`](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
  [`fia_obs_test`](https://finnverse.github.io/FINN/reference/fia_obs_dt.md)
  : Oregon FIA calibration data: training sites
- [`fia_env_dt`](https://finnverse.github.io/FINN/reference/fia_env_dt.md)
  [`fia_env_test`](https://finnverse.github.io/FINN/reference/fia_env_dt.md)
  : Oregon FIA calibration data: site climate
- [`fia_init_trees`](https://finnverse.github.io/FINN/reference/fia_init_trees.md)
  [`fia_init_test`](https://finnverse.github.io/FINN/reference/fia_init_trees.md)
  : Oregon FIA calibration data: initial tree list
- [`fia_tree_dt`](https://finnverse.github.io/FINN/reference/fia_tree_dt.md)
  [`fia_tree_test`](https://finnverse.github.io/FINN/reference/fia_tree_dt.md)
  : Oregon FIA calibration data: the trees behind the mortality counts
- [`fia_species_dt`](https://finnverse.github.io/FINN/reference/fia_species_dt.md)
  : Oregon FIA calibration data: species coding

## Environmental scaling

FINN standardises environmental predictors internally
(`env_autoscale = TRUE`); these expose the stored constants.

- [`apply_env_scaling()`](https://finnverse.github.io/FINN/reference/apply_env_scaling.md)
  : Apply stored z-standardization to environmental predictors
- [`compute_env_scaling()`](https://finnverse.github.io/FINN/reference/compute_env_scaling.md)
  : Compute per-predictor environmental scaling

## Interpretation (xAI)

Read back what a fitted - or hybrid - model actually learned.

- [`plot(`*`<finn_class>`*`)`](https://finnverse.github.io/FINN/reference/plot.finn_class.md)
  : Plot the fitted parameters of a FINN model
- [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) :
  Accumulated local effect plots
- [`plot(`*`<FINNale>`*`)`](https://finnverse.github.io/FINN/reference/plot.FINNale.md)
  : Plot ALE curves of a FINN model
- [`conditionalEffects()`](https://finnverse.github.io/FINN/reference/conditionalEffects.md)
  : Conditional effects of a FINN model
- [`averageConditionalEffects()`](https://finnverse.github.io/FINN/reference/averageConditionalEffects.md)
  : Average conditional effects of a FINN model
- [`feature_importance()`](https://finnverse.github.io/FINN/reference/feature_importance.md)
  : Permutation feature importance for FINN demographic rates

## Low-level utilities

Building blocks used internally and exported for advanced use.

- [`binomial_from_bernoulli()`](https://finnverse.github.io/FINN/reference/binomial_from_bernoulli.md)
  : Draw binomial counts from per-trial Bernoulli probabilities
- [`binomial_from_gamma()`](https://finnverse.github.io/FINN/reference/binomial_from_gamma.md)
  : Draw deaths from a Gamma matched to the binomial
- [`binomial_from_poisson()`](https://finnverse.github.io/FINN/reference/binomial_from_poisson.md)
  : Draw deaths as a Bernoulli/Poisson count
- [`rnbinom_torch_refined()`](https://finnverse.github.io/FINN/reference/rnbinom_torch_refined.md)
  : Negative-binomial recruit draw without the recruitment floor
- [`rweibull_cohorts()`](https://finnverse.github.io/FINN/reference/rweibull_cohorts.md)
  : Generate Cohorts Using Weibull Distribution
