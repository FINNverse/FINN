#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List obsDF2arraysCpp(DataFrame obs_dt, CharacterVector extra_cols) {
  // Extract columns from the data frame
  IntegerVector siteID = obs_dt["siteID"];
  IntegerVector patchID = obs_dt["patchID"];
  IntegerVector cohortID = obs_dt["cohortID"];
  IntegerVector species = obs_dt["species"]; // Change to IntegerVector
  NumericVector dbh = obs_dt["dbh"];
  NumericVector trees = obs_dt["trees"];

  // Get the dimensions
  int Nsites = max(siteID);
  int Npatches = max(patchID);
  int maxNcohorts = max(cohortID);
  int n = siteID.size();

  // Initialize vectors to hold the data
  IntegerVector species_vec(Nsites * Npatches * maxNcohorts, NA_INTEGER); // Change to IntegerVector
  NumericVector dbh_vec(Nsites * Npatches * maxNcohorts, NA_REAL);
  NumericVector trees_vec(Nsites * Npatches * maxNcohorts, NA_REAL);

  // Populate the vectors
  for (int i = 0; i < n; ++i) {
    int index = (siteID[i] - 1) + (patchID[i] - 1) * Nsites + (cohortID[i] - 1) * Nsites * Npatches;
    species_vec[index] = species[i];
    dbh_vec[index] = dbh[i];
    trees_vec[index] = trees[i];
  }

  // Create the result list
  List result = List::create(
    Named("species") = species_vec,
    Named("dbh") = dbh_vec,
    Named("trees") = trees_vec,
    Named("dim") = IntegerVector::create(Nsites, Npatches, maxNcohorts)
  );

  // Handle additional columns
  for (int j = 0; j < extra_cols.size(); ++j) {
    String col_name = extra_cols[j];
    SEXP col = obs_dt[col_name];

    if (Rf_isInteger(col)) {
      IntegerVector col_vec(Nsites * Npatches * maxNcohorts, NA_INTEGER);
      IntegerVector col_data = as<IntegerVector>(col);
      for (int i = 0; i < n; ++i) {
        int index = (siteID[i] - 1) + (patchID[i] - 1) * Nsites + (cohortID[i] - 1) * Nsites * Npatches;
        col_vec[index] = col_data[i];
      }
      result[col_name] = col_vec;
    } else if (Rf_isReal(col)) {
      NumericVector col_vec(Nsites * Npatches * maxNcohorts, NA_REAL);
      NumericVector col_data = as<NumericVector>(col);
      for (int i = 0; i < n; ++i) {
        int index = (siteID[i] - 1) + (patchID[i] - 1) * Nsites + (cohortID[i] - 1) * Nsites * Npatches;
        col_vec[index] = col_data[i];
      }
      result[col_name] = col_vec;
    } else if (Rf_isString(col)) {
      CharacterVector col_vec(Nsites * Npatches * maxNcohorts, NA_STRING);
      CharacterVector col_data = as<CharacterVector>(col);
      for (int i = 0; i < n; ++i) {
        int index = (siteID[i] - 1) + (patchID[i] - 1) * Nsites + (cohortID[i] - 1) * Nsites * Npatches;
        col_vec[index] = col_data[i];
      }
      result[col_name] = col_vec;
    }
  }

  return result;
}
