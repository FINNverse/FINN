
#' Initialize a CohortMat Instance
#'
#' A convenience function that calls `CohortMat$new` to initialize a new `CohortMat` object. This function has the same parameters and behavior as the `initialize` method of the `CohortMat` class.
#'
#' @param obs_df A data frame containing columns "siteID", "patchID", "species", "dbh", and "trees". If provided, it will be used to initialize the tensors.
#' @param dbh A tensor or array representing the diameter at breast height. Defaults to `NULL`.
#' @param trees A tensor or array representing the number of trees. Defaults to `NULL`.
#' @param species A tensor or array representing the species. Defaults to `NULL`.
#' @param dims A numeric vector representing the dimensions of the arrays (sites, patches, cohorts). Defaults to `c(50, 30, 10)`.
#' @param sp An integer representing the number of species. Defaults to `10`.
#' @param device A character string specifying the device to use ('cpu' or 'cuda'). Defaults to `'cpu'`.
#'
#' @return A new instance of the `CohortMat` class.
#' @examples
#' \dontrun{
#' # Initialize a CohortMat instance
#' cohort <- initCohort(obs_df = your_data_frame)
#' }
#' @export
initCohort <- function(obs_df = NULL, dbh = NULL, trees = NULL, species = NULL, dims = c(50, 30, 10), sp = 10, device = 'cpu') {
  CohortMat$new(obs_df = obs_df, dbh = dbh, trees = trees, species = species, dims = dims, sp = sp, device = device)
}
