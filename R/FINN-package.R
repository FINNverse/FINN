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
#' @importFrom abind abind
#' @importFrom coro loop
#' @importFrom utils modifyList str capture.output
#' @importFrom stats model.matrix runif predict
#' @importFrom graphics matplot
"_PACKAGE"

# data.table uses non-standard evaluation, so column names referenced inside
# `dt[...]` read to R CMD check as undefined global variables. Declaring them
# here silences the (harmless) "no visible binding" NOTEs. Also lists a couple of
# package-internal objects created dynamically via torch::nn_module().
utils::globalVariables(c(
  ".", "..cols_env", "..out_cols_obs_dt", "..out_cols_obs_dt_aggr",
  "..out_cols_siteID_dt", "..spc", "NpatchesPerSite", "NyearsPerPatch",
  "OrigYear", "ale", "ba", "cohortID", "common", "complete", "day", "dbh",
  "dbh_before", "dbh_growth", "gPSize", "growth_n", "i.is_complete",
  "importance", "living", "living_before", "mort", "n_at_risk", "n_died",
  "parHeight", "patchID", "patchName", "period", "period_length", "process",
  "rate_sd", "reg", "rel_growth", "sameYearsAllPatches", "shade", "siteID",
  "siteName", "species", "species_before", "species_name", "status", "td",
  "treeName", "trees", "trees_ha", "value", "var", "variable", "year_before",
  "years", "TransformerEncoderLayer"
))
