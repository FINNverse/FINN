#' @title The input tables of FINN: observations, environment, initial trees
#'
#' @description
#' [FINN::fit()] and [FINN::simulateForest()] take their input as a few plain
#' tables (`data.frame` or `data.table`), all keyed by an integer site index
#' `siteID` and a `year`. This page says what each table must contain, what
#' may be missing, and what FINN fills in by itself. The bundled Oregon FIA
#' data ([fia_obs_dt], [fia_env_dt], [fia_init_trees]) are a complete example
#' and the vignette *Preparing your data for FINN*
#' (`vignette("C-Data_preparation", package = "FINN")`) builds them from a raw
#' tree list.
#'
#' The short version: build `data` and `env` with [FINN::makeObsData()] and
#' [FINN::resolveSiteIDs()], the initial trees with [FINN::makeInitCohorts()],
#' and pass them to `fit(data = , env = , init_cohort = )`. The three helpers
#' produce exactly the columns below.
#'
#' @section How the tables fit together:
#' * **`env`** is the time axis. Every row is one site in one simulated time
#'   step (usually one year), and every site has the same years. The first row
#'   of a site drives the first step, starting from the site's initial trees.
#' * **`init_cohort`** is the starting state: the trees of every site at the
#'   first census, grouped into cohorts. The number of patches per site is read
#'   from it.
#' * **`data`** is what the simulation is compared with: for every site,
#'   census year and species, the observed stand variables and demographic
#'   rates. An observation with a given `year` is compared with the simulated
#'   state at the end of the step with that `year` in `env`, so every `year` in
#'   `data` must also be a `year` of `env`.
#' * **`disturbance`** (optional) removes trees in given years.
#' * **`tree_data`** (optional) holds the individual trees, for the two ways of
#'   scoring growth and mortality per tree ([finn_loss]).
#'
#' The sites of `data`, `env` and `init_cohort` must be the same set; the ids
#' need not be `1..N`. Species are coded `1..N_species` in every table.
#' [FINN::resolveSiteIDs()] turns site names and species names into these
#' codes and returns the lookup tables.
#'
#' @section `data` -- the observations:
#' One row per `siteID`, `year` and `species`. Columns:
#'
#' | column | what it is |
#' |---|---|
#' | `siteID` | integer site index |
#' | `year` | the census year; must be a `year` of `env` |
#' | `species` | integer species code `1..N_species` |
#' | `dbh` | mean diameter at breast height of the living trees of that species, cm |
#' | `ba` | basal area of the living trees, m2 per patch (the mean over the site's patches) |
#' | `trees` | number of living trees per patch (the mean over the patches; may be fractional) |
#' | `growth` | relative diameter growth of the survivors per time step, `dbh / dbh_before - 1` (the increment between two censuses when one step is one census interval; for intervals of several steps see `finn(growth_rate = )`) |
#' | `mort` | share of the trees alive at the previous census that died by this one (0 to 1) |
#' | `reg` | recruits (trees that grew past `recruits_dbh`) since the previous census, per hectare |
#'
#' `dbh`, `ba` and `trees` describe the stand at the census; `growth`, `mort`
#' and `reg` describe what happened since the previous census. Optional
#' columns, filled in by FINN when absent:
#'
#' | column | what it is | when absent |
#' |---|---|---|
#' | `period_length` | the number of time steps since the previous census (the interval the rates refer to) | one step: each rate is compared with a single simulated step |
#' | `n_at_risk` | trees alive at the previous census behind `mort` | 1: every mortality observation weighs the same |
#' | `growth_n` | trees behind the `growth` mean | 1: every growth observation weighs the same |
#'
#' With `period_length`, the simulated growth and mortality of the interval's
#' steps are averaged and the simulated recruitment is summed before they are
#' compared with the observation. `n_at_risk` and `growth_n` weight the
#' mortality and growth terms by the number of trees behind each observation
#' (a mortality rate from 100 trees says more than one from a single tree).
#' [FINN::makeObsData()] produces all of them. Further optional columns exist
#' for the size-class and quantile aggregations (`n_class_1`, ..., `dbh_q10`,
#' ...; see [finn_loss]).
#'
#' **Missing values.** Every response column must exist, but any value may be
#' `NA`: a response that was not observed in a cell (site, year, species) is
#' simply left out of the comparison for that cell, and the other responses
#' of the cell are still used. A species that does not occur on a site needs
#' a row with 0 for `ba`, `trees` and `reg` (so that FINN learns it is
#' absent) and `NA` for `dbh`, `growth` and `mort` (there is no tree to
#' measure). [FINN::makeObsData()] writes these rows.
#'
#' **Incomplete series.** Sites need not share census years, and a site may
#' have fewer censuses than another: a (site, year) with no row is treated as
#' unobserved in that year, and `period_length` may differ between sites and
#' censuses. When it does, every site is scored over its own interval. A
#' census with `NA` `period_length` next to one with a value counts as a
#' one-step interval on the usual path; on the per-site path (interval
#' lengths that differ between sites, `loss_aggregation = "tree"`) every
#' window is read from the column and a site without a finite value stops
#' the fit with a message. Years of `env` before a site's first census and
#' after its last are simulated but not scored.
#'
#' **Units.** Use the units of your inventory consistently: `dbh` in cm,
#' areas in m2 per patch of `patch_size` ha, recruits per hectare. The
#' per-hectare quantities depend on `patch_size`, so pass the plot size of
#' the inventory to `fit(patch_size = )`.
#'
#' ```r
#' head(fia_obs_dt[, .(siteID, year, species, dbh, ba, trees, growth, mort, reg, n_at_risk, growth_n)], 2)
#' #>    siteID  year species    dbh    ba trees growth  mort   reg n_at_risk growth_n
#' #> 1:      1     1       1 92.837 0.342   0.5  0.036     0     0         2        2
#' #> 2:      1     1       2     NA 0.000   0.0     NA    NA     0         0        0
#' ```
#'
#' @section `env` -- the environment of every simulated step:
#' One row per `siteID` and `year`, with one column per predictor and nothing
#' else: no patch column, no response. Every site must have the same years,
#' every (site, year) a value of every predictor, and the years must run from
#' the first step to simulate to the last observation in `data`; a year that
#' is missing or `NA` cannot be simulated and stops the fit with a message.
#' At least one predictor is needed; a process whose `formula` names none has
#' nothing to respond to.
#'
#' The predictors are given in their natural units. `fit(env_autoscale =
#' TRUE)` (the default) standardises every column to mean 0 and standard
#' deviation 1 and stores the constants on the model, so the same raw table
#' layout works for [FINN::simulateForest()] afterwards. For the bundled FIA
#' data the predictors are the bioclimatic variables of each plot:
#'
#' ```r
#' head(fia_env_dt, 3)
#' #>    siteID  year  temp tempmax tempmin  prec precseas precwarmq
#' #> 1:      1     0  3.96    24.7   -10.4   547     30.0       101
#' #> 2:      1     1  3.96    24.7   -10.4   547     30.0       101
#' #> 3:      1     2  3.96    24.7   -10.4   547     30.0       101
#' ```
#'
#' A climate that does not change between years is repeated for every year,
#' as here; a year-by-year climate gives every step its own row. Which
#' columns a process uses is set by `createProcess(formula = )`, e.g.
#' `~ temp + prec`; the default `~.` uses every predictor column.
#'
#' @section `disturbance` -- optional:
#' One row per `siteID` and `year`, over exactly the sites and years of `env`,
#' with the column `intensity`: the probability (0 to 1) that a patch of the
#' site loses all its trees in that year. Every patch is drawn separately, so
#' an intensity of 0.1 clears on average one patch in ten; use 0 where nothing
#' happens. `NULL` (the default) means no disturbance.
#'
#' @section `init_cohort` -- the initial trees:
#' A `CohortMat`, the trees of every site at the first census grouped into
#' cohorts of one species and one diameter per patch, built by
#' [FINN::makeInitCohorts()] from a table with `siteID`, `patchID`, `species`,
#' `dbh` and optionally `trees` (stems the row stands for, default 1).
#' `makeInitCohorts(dbh_binsize = )` bins the diameters into fewer, cheaper
#' cohorts; without it every tree is a cohort, which the per-tree scoring
#' needs. The number of patches per site is read from `init_cohort`, so
#' `fit(patches = )` is ignored when it is given. Without `init_cohort` every
#' site starts from bare ground with `patches` empty patches -- a succession
#' run, not a calibration.
#'
#' @section `tree_data` -- optional, the individual trees:
#' One row per tree and census, needed only to score growth or mortality per
#' tree. Always: `siteID`, `year` (the census that closes the tree's
#' interval, the `year` the same interval has in `data`), `patchID`,
#' `species`, `growth` (its relative growth rate per time step, `NA` if not
#' observed) and `died` (0/1 over the interval, `NA` if not observed). For
#' `loss_aggregation = "tree"` also `dbh` at the start of the interval, and
#' optionally `trees`; for `loss_resolution = "tree"` also `cohortID`, the
#' cohort the tree started as (`makeInitCohorts(treeID_table = TRUE)` returns
#' the table that carries it). [finn_loss] explains the two routes.
#'
#' @section The helper functions:
#' * [FINN::makeObsData()] turns a tree list (one row per tree and inventory,
#'   with site, patch, year, species name, dbh and status) into the
#'   observation table, per patch or aggregated per site, with the units above.
#' * [FINN::resolveSiteIDs()] replaces site, patch and species names by the
#'   integer codes, aligns the observation, environment and tree tables on
#'   them, and can build the initial cohorts in the same call.
#' * [FINN::makeInitCohorts()] builds `init_cohort` from the trees at the first
#'   census.
#' * [FINN::climateDF2array()], [FINN::obsDF2arrays()] and [FINN::CohortMat()]
#'   are the lower-level converters behind them.
#'
#' @seealso [FINN::fit()], [finn_loss], [fia_obs_dt], [fia_env_dt],
#'   [fia_init_trees], `vignette("C-Data_preparation", package = "FINN")`
#' @name finn_data
#' @aliases finn-data
NULL
