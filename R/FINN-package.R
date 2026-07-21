#' FINN: Forest Informed Neural Networks
#'
#' FINN is a differentiable forest gap model. A forest is represented as cohorts
#' of same-species, same-size trees that are updated each timestep by four
#' demographic processes — competition, growth, mortality and regeneration. Each
#' process can be a mechanistic function, a neural network, or a mixture of the
#' two, and the whole model is calibrated end-to-end by gradient descent through
#' the simulation (implemented in \pkg{torch}).
#'
#' @section Getting started:
#' \itemize{
#'   \item [finn()] assembles a model from one process per component, each built
#'     with [createProcess()] (mechanistic) or [createHybrid()] (neural network).
#'   \item [simulateForest()] runs a model forward from known parameters.
#'   \item [fit()] calibrates a model to data; [predict.finn_class()] scores it.
#'   \item [makeObsData()], [resolveSiteIDs()] and [makeInitCohorts()] turn a raw
#'     tree list into FINN's input tables.
#'   \item [ALE()], [summary.finn_class()], [feature_importance()] and
#'     [conditionalEffects()] interpret a fitted model.
#' }
#'
#' The vignettes give a guided tour: `browseVignettes("FINN")`.
#'
#' @keywords internal
"_PACKAGE"
