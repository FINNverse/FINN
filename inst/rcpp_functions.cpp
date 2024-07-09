#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector populate_array(NumericMatrix climate_dt, IntegerVector site_ids, IntegerVector year_ids, bool include_month, int Nsites, int Nyears, int Nmonths, int Nenv) {
  int total_rows = climate_dt.nrow();
  NumericVector env_array;

  if (include_month) {
    env_array = NumericVector(Nsites * Nyears * Nmonths * Nenv, NA_REAL);
  } else {
    env_array = NumericVector(Nsites * Nyears * Nenv, NA_REAL);
  }

  for (int i = 0; i < total_rows; ++i) {
    int site = -1;
    int year = -1;

    for (int j = 0; j < site_ids.size(); ++j) {
      if (site_ids[j] == climate_dt(i, 2)) {
        site = j;
        break;
      }
    }

    for (int j = 0; j < year_ids.size(); ++j) {
      if (year_ids[j] == climate_dt(i, 0)) {
        year = j;
        break;
      }
    }

    for (int e = 0; e < Nenv; ++e) {
      if (include_month) {
        int month = climate_dt(i, 1) - 1;
        env_array(site + Nsites * (year + Nyears * (month + Nmonths * e))) = climate_dt(i, 3 + e);
      } else {
        env_array(site + Nsites * (year + Nyears * e)) = climate_dt(i, 3 + e);
      }
    }
  }

  return env_array;
}
