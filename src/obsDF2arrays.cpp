#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List obsDF2arraysCpp(DataFrame obs_dt) {
  IntegerVector siteID = obs_dt["siteID"];
  IntegerVector patchID = obs_dt["patchID"];
  IntegerVector cohortID = obs_dt["cohortID"];
  CharacterVector species = obs_dt["species"];
  NumericVector dbh = obs_dt["dbh"];
  NumericVector trees = obs_dt["trees"];

  int Nsites = max(siteID);
  int Npatches = max(patchID);
  int maxNcohorts = max(cohortID);
  int n = siteID.size();

  CharacterVector species_vec(Nsites * Npatches * maxNcohorts, NA_STRING);
  NumericVector dbh_vec(Nsites * Npatches * maxNcohorts, NA_REAL);
  NumericVector trees_vec(Nsites * Npatches * maxNcohorts, NA_REAL);

  for (int i = 0; i < n; ++i) {
    int index = (siteID[i] - 1) + (patchID[i] - 1) * Nsites + (cohortID[i] - 1) * Nsites * Npatches;
    species_vec[index] = species[i];
    dbh_vec[index] = dbh[i];
    trees_vec[index] = trees[i];
  }

  return List::create(
    Named("species") = species_vec,
    Named("dbh") = dbh_vec,
    Named("trees") = trees_vec,
    Named("dim") = IntegerVector::create(Nsites, Npatches, maxNcohorts)
  );
}
