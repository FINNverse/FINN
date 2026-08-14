#' NW-FVA growth-and-yield tables (bundled reference data)
#'
#' Loads the growth-and-yield reference table shipped with FINN: the NW-FVA
#' "gestaffelte Durchforstung 2021" yield tables for the five main commercial
#' species of northwest Germany (oak, beech, spruce, Douglas-fir, pine), tidied
#' into one long data frame with a few derived growth/thinning quantities. These
#' give FINN real, regionally valid growth-and-yield information out of the box
#' (e.g. to parameterise or benchmark the growth process and the management layer).
#'
#' @details
#' One row per species (`Art`) x relative yield class (`Ekl`, -1..4) x age
#' (`Alter`, 5-year steps). Remaining stand after thinning: `N` (stems/ha), `Hg`
#' (mean height, m), `H100` (top height, m), `G` (basal area, m2/ha), `Dg`
#' (quadratic mean dbh, cm), `Dw` (dbh of the mean-basal-area tree, cm), `V`
#' (volume, m3/ha). Removed (thinning) stand: `N_aus`, `G_aus`, `Dg_aus`,
#' `V_aus`. Performance: `iV`, `GWL`, `dGZ`. Derived here: `dDg_dt` (cm/yr),
#' `thin_N_frac`, `thin_G_frac`, `Dg_aus_rel` (= `Dg_aus/Dg`; < 1 thinning from
#' below, > 1 from above).
#'
#' @section Source and license:
#' Derived from Albert, M., Nagel, J., Schmidt, M., Nagel, R.-V., Spellmann, H.
#' (2021) *Eine neue Generation von Ertragstafeln fuer Eiche, Buche, Fichte,
#' Douglasie und Kiefer* \[Dataset\], Zenodo,
#' \doi{10.5281/zenodo.6343906}, licensed **CC-BY-4.0** (copyright NW-FVA);
#' extracted with the R package *et.nwfva* (Nuske, Staupendahl, Albert 2022,
#' \doi{10.5281/zenodo.7207597}, GPL >= 2). The CC-BY-4.0 data is redistributed
#' under its own license, separate from FINN's GPL-3. Changes made (tidying,
#' derived columns) are listed in `inst/extdata/nwfva_gy_tables_SOURCE.md`. Please
#' cite Albert et al. (2021) when using these tables.
#'
#' @return A `data.frame` of the yield-table rows described above.
#' @examples
#' gy <- nwfva_yield_tables()
#' subset(gy, Art == "Fichte" & Ekl == 1, c("Alter", "Dg", "dDg_dt"))
#' @export
nwfva_yield_tables <- function() {
  path <- system.file("extdata", "nwfva_gy_tables.csv", package = "FINN")
  if (!nzchar(path))
    stop("nwfva_gy_tables.csv not found in the installed FINN package.",
         call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}
