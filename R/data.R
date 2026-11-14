#' Oregon FIA calibration data: training sites
#'
#' A ready-to-fit sample of the US Forest Inventory and Analysis (FIA) network
#' in Oregon: 200 training sites (`fia_*_dt` / `fia_init_trees`) and 200
#' disjoint holdout sites (`fia_*_test`) that are never seen during fitting, so
#' a model fit on the training tables can be scored genuinely out of sample.
#' The tables are already in FINN's input format -- IDs resolved, species
#' re-indexed `1..K` -- so they go straight into [fit()] / [simulateForest()]
#' without [makeObsData()] or [resolveSiteIDs()]. The `D-Fit_to_FIA` and
#' `E-Mortality` vignettes are the worked examples.
#'
#' All seven objects are `data.table`s. Attach `data.table` to use its `[`
#' semantics on them, and `data.table::copy()` one before modifying it by
#' reference (`:=`), as with any package dataset.
#'
#' The tables are on an **annual** time axis: one simulation step is one year.
#' The first inventory (2001-02) is the initial state, year 0; the two
#' remeasurements are observations at the number of years since it, 10 and 20
#' (21 for one training site and 22 for one holdout site, whose third
#' inventory came late), each with `period_length`, the years since the
#' previous inventory. Growth and mortality are yearly rates; [fit()] compares
#' them with the mean of the simulated yearly rates over the `period_length`
#' years before each inventory, and compares `reg` with the recruits simulated
#' over those years. Before FINN 0.2.0 the tables had one step per ten-year
#' interval (`year` 1 and 2, rates over the interval); `makeObsData(rates =
#' "interval")` and `resolveSiteIDs(time_axis = "census")` still build such
#' tables.
#'
#' @format `fia_obs_dt` and `fia_obs_test`: one row per site x year x species
#'   (4400 rows, 14 columns). The year-0 initial state is in `fia_init_trees` /
#'   `fia_init_test`.
#'   \describe{
#'     \item{`siteID`, `year`}{integer site index (`1..200`) and the years since
#'       the site's first inventory (`10`, `20`; `21` or `22` for one site
#'       each).}
#'     \item{`period_length`}{years since the previous inventory (`10`; `11`
#'       or `12` for the late inventory): the window over which the simulation
#'       is aggregated before it meets the observation.}
#'     \item{`species`, `species_name`}{species index `1..11` and its name;
#'       `"other"` pools the species outside the 98 % most abundant stems.}
#'     \item{`ba`}{living basal area, m2 -- the mean over the site's patches
#'       (0.06 ha each; `makeObsData(aggregate_by_site = TRUE)`).}
#'     \item{`trees`}{living stems, mean per patch (hence fractional).}
#'     \item{`dbh`}{mean diameter at breast height of living trees, cm (`NA`
#'       when there are none).}
#'     \item{`growth`, `growth_n`}{mean annual relative dbh growth rate of the
#'       survivors, `(dbh / dbh_before)^(1 / period_length) - 1` per tree
#'       (tree-weighted over patches), and the number of trees behind it.}
#'     \item{`n_at_risk`, `n_died`, `mort`}{trees alive at the start of the
#'       interval and how many of them had died by its end, summed over
#'       patches, and the annual mortality rate
#'       `1 - (1 - n_died / n_at_risk)^(1 / period_length)`. The binomial
#'       mortality loss scores `mort` weighted by `n_at_risk`.}
#'     \item{`reg`}{recruits over the interval, trees per hectare (mean over
#'       patches).}
#'   }
#' @source Subsampled by `dev/make_extdata.R` (`set.seed(42)`) from the Oregon
#'   FIA plots prepared in the (private) `FINNverse/FINN-fia` analysis repo:
#'   tree records via the FIA
#'   database, climate from WorldClim bioclim variables attached per plot. The
#'   full provenance chain is in `data-raw/README.md`.
#' @seealso [fia_env_dt], [fia_init_trees], [fia_species_dt]
#' @family fia datasets
#' @examples
#' library(data.table)
#' fia_obs_dt[year == 10, .(sites = uniqueN(siteID), ba = mean(ba)), by = species_name]
#' # the observation windows: inventory year x interval length, sites
#' unique(fia_obs_dt[, .(siteID, year, period_length)])[, .N, by = .(year, period_length)]
"fia_obs_dt"

#' @rdname fia_obs_dt
"fia_obs_test"

#' Oregon FIA calibration data: site climate
#'
#' Climate per site and year for the [fia_obs_dt] training sites (`fia_env_dt`)
#' and the [fia_obs_test] holdout sites (`fia_env_test`), in **natural units**
#' (not standardised -- [fit()] z-scales internally with `env_autoscale`).
#' One row per simulated year `1..22`: the row of year `k` drives the year
#' that ends `k` years after the first inventory, so every site x year in the
#' observation tables is covered, and the years after a site's last
#' inventory are simulated but not scored. The climate is the plot's
#' WorldClim normal, the same in every year.
#'
#' @format A `data.table` with 4400 rows and 8 columns: 200 sites x 22 years,
#'   one row per site x year.
#'   \describe{
#'     \item{`siteID`, `year`}{site index and simulated year (`1..22`).}
#'     \item{`temp`, `tempmax`, `tempmin`}{annual mean, warmest-month maximum
#'       and coldest-month minimum temperature, degrees C (bioclim 1, 5, 6).}
#'     \item{`prec`, `precwarmq`}{annual and warmest-quarter precipitation, mm
#'       (bioclim 12, 18).}
#'     \item{`precseas`}{precipitation seasonality, coefficient of variation
#'       (bioclim 15).}
#'   }
#' @inherit fia_obs_dt source
#' @family fia datasets
"fia_env_dt"

#' @rdname fia_env_dt
"fia_env_test"

#' Oregon FIA calibration data: initial tree list
#'
#' The year-0 tree list of the [fia_obs_dt] training sites (`fia_init_trees`)
#' and the [fia_obs_test] holdout sites (`fia_init_test`): what a site looked
#' like at the first census, from which [makeInitCohorts()] builds the initial
#' cohorts a fit or simulation starts from
#' (`makeInitCohorts(fia_init_trees, Nspecies = 11)`).
#'
#' @format A `data.table` with one row per tree (6039 / 6083 rows) and 28
#'   columns. The ones FINN uses:
#'   \describe{
#'     \item{`siteID`, `patchID`, `year`}{site and patch index; `year` is 0.}
#'     \item{`species`, `species_name`}{species index `1..11` and name, coded
#'       as in [fia_species_dt].}
#'     \item{`dbh`}{diameter at breast height, cm.}
#'     \item{`trees`}{trees represented by the record (`1`; `NA` for records
#'       that are not living stems).}
#'     \item{`treeName`}{unique tree identifier.}
#'   }
#'   The remaining columns (`siteName`, `patchName`, `OrigYear`, FIA's `SPCD`
#'   and `STATUSCD`, `status`, `living`, and the per-tree change fields such as
#'   `dbh_before`, `dbh_growth`, `mort`, `reg`) are carried over from the
#'   remeasurement tree list produced by [makeObsData()]; at year 0 most of
#'   them are `NA` by construction.
#' @inherit fia_obs_dt source
#' @family fia datasets
"fia_init_trees"

#' @rdname fia_init_trees
"fia_init_test"

#' Oregon FIA calibration data: species coding
#'
#' The species index shared by all `fia_*` tables. Derived from the training
#' sites, so `species` is `1..11` in every table.
#'
#' @format A `data.table` with 11 rows and 2 columns, `species` (integer index)
#'   and `species_name` (Latin name, or `"other"` for the pooled rare species).
#' @inherit fia_obs_dt source
#' @family fia datasets
"fia_species_dt"
