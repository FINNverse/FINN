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
#' The functions in the order of a typical analysis:
#' \enumerate{
#'   \item **Prepare the data.** [makeObsData()], [resolveSiteIDs()] and
#'     [makeInitCohorts()] turn a raw tree list into FINN's input tables: the
#'     observations, the environment of every year and the initial cohorts
#'     ([finn_data]).
#'   \item **Build the model.** [finn()] assembles one process per component
#'     (growth, mortality, regeneration, competition), each built with
#'     [createProcess()] (mechanistic, with a linear or neural-network response
#'     to the environment) or [createHybrid()] (a neural network).
#'   \item **Fit it.** [fit()] estimates all processes jointly against all
#'     observed responses; its argument `holdout` evaluates the model right
#'     after the fit.
#'   \item **Evaluate it.** [evaluate()] compares simulations with the
#'     observations of the training data and a holdout and returns the result
#'     (also kept in the model); [comparePlot()] draws observed against
#'     simulated, as a scatter page, a heatmap per species or a time series;
#'     [summary.finn_class()] shows the skill, the convergence of the fit and
#'     the effects of the environment.
#'   \item **Look inside.** [plot.finn_class()] draws the fitted parameters,
#'     [ALE()] the effects of the environment, also of a neural network.
#'   \item **Simulate.** [simulateForest()] runs the model forward, also from
#'     bare ground (`bare_ground = TRUE`); [plot.finn_simulation()] draws the
#'     result.
#' }
#'
#' The vignettes give a guided tour: `browseVignettes("FINN")`.
#'
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # the bundled Oregon FIA data: observations, environment and initial trees
#' Nsp   <- max(fia_obs_dt$species)
#' train <- 1:4; hold <- 5:6
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% train], Nspecies = Nsp)
#'
#' # build the model: the environmental responses as linear functions of
#' # temperature and precipitation
#' FINN.seed(1)
#' m <- finn(N_species = Nsp, recruits_dbh = 12.9,
#'           growth_process       = createProcess(~ temp + prec, growth),
#'           mortality_process    = createProcess(~ temp + prec, mortality),
#'           regeneration_process = createProcess(~ temp + prec, regeneration_capped))
#'
#' # fit all processes jointly (a real calibration runs hundreds of epochs)
#' fit(m, data = fia_obs_dt[siteID %in% train], env = fia_env_dt[siteID %in% train],
#'     init_cohort = init, epochs = 3L, patch_size = 0.06, plot_progress = FALSE)
#'
#' # evaluate on two sites the model has not seen, with a bare-ground run
#' ev <- evaluate(m, holdout = list(
#'   data = fia_obs_dt[siteID %in% hold], env = fia_env_dt[siteID %in% hold],
#'   init_cohort = makeInitCohorts(fia_init_trees[siteID %in% hold], Nspecies = Nsp)),
#'   bare_ground = TRUE, years = 100)
#' plot(ev)                         # observed against simulated
#' plot(ev, type = "bare_ground")   # succession from an empty stand
#'
#' # skill, convergence and environmental effects; the fitted parameters
#' summary(m)
#' plot(m)
#'
#' @keywords internal
#' @import ggplot2
#' @importFrom abind abind
#' @importFrom coro loop
#' @importFrom utils modifyList str capture.output getFromNamespace
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
  "OrigYear", "ale", "ba", "cohortID", "day", "dbh",
  "dbh_before", "dbh_growth", "gPSize", "growth_n",
  "importance", "living", "living_before", "mort", "n_at_risk", "n_died",
  "parHeight", "patchID", "patchName", "period", "period_length", "process",
  "rate_sd", "reg", "rel_growth", "sameYearsAllPatches", "shade", "siteID",
  "siteName", "species", "species_before", "species_name", "status", "td",
  "treeName", "trees", "trees_ha", "value", "var", "variable", "year_before",
  "TransformerEncoderLayer", ".env_row", "first", "growth",
  "period_length_in", "pl_rate", "step",
  # plot.finn_class() and plot.finn_simulation()
  "at_bound", "hi", "lo", "lower", "n_sites", "panel", "predictor", "upper", "x", "year",
  # comparePlot(), evaluate(), summary()
  "auc", "census", "deaths_ratio", "label", "max_year", "obs", "p_interval", "pred",
  "r2", "species_lab", "total_ba", "years", "ratio", "ratio_w",
  "fpr", "tpr", "name", "share", "bound", "total", "colour_group", "site_lab", "value_max",
  "dbh_end", "dbh_max", "died", "any_obs", "n", "step", "obs_sum", "tile_label"
))
