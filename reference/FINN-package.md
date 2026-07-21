# FINN: Forest Informed Neural Networks

FINN is a differentiable forest gap model. A forest is represented as
cohorts of same-species, same-size trees that are updated each timestep
by four demographic processes — competition, growth, mortality and
regeneration. Each process can be a mechanistic function, a neural
network, or a mixture of the two, and the whole model is calibrated
end-to-end by gradient descent through the simulation (implemented in
torch).

## Getting started

- [`finn()`](https://finnverse.github.io/FINN/reference/finn.md)
  assembles a model from one process per component, each built with
  [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)
  (mechanistic) or
  [`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md)
  (neural network).

- [`simulateForest()`](https://finnverse.github.io/FINN/reference/simulateForest.md)
  runs a model forward from known parameters.

- [`fit()`](https://finnverse.github.io/FINN/reference/fit.md)
  calibrates a model to data;
  [`predict.finn_class()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)
  scores it.

- [`makeObsData()`](https://finnverse.github.io/FINN/reference/makeObsData.md),
  [`resolveSiteIDs()`](https://finnverse.github.io/FINN/reference/resolveSiteIDs.md)
  and
  [`makeInitCohorts()`](https://finnverse.github.io/FINN/reference/makeInitCohorts.md)
  turn a raw tree list into FINN's input tables.

- [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md),
  [`summary.finn_class()`](https://finnverse.github.io/FINN/reference/summary.finn_class.md),
  [`feature_importance()`](https://finnverse.github.io/FINN/reference/feature_importance.md)
  and
  [`conditionalEffects()`](https://finnverse.github.io/FINN/reference/conditionalEffects.md)
  interpret a fitted model.

The vignettes give a guided tour: `browseVignettes("FINN")`.

## See also

Useful links:

- <https://github.com/FINNverse/FINN>

- <https://finnverse.github.io/FINN/>

- Report bugs at <https://github.com/FINNverse/FINN/issues>

## Author

**Maintainer**: Yannek K\<c3\>\<a4\>ber <y.kaeber@posteo.de>
([ORCID](https://orcid.org/0000-0002-7041-9849))

Authors:

- Maximilian Pichler <maximilian.pichler@biologie.uni-regensburg.de>
  ([ORCID](https://orcid.org/0000-0003-2252-8327))
