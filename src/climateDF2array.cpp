#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector climateDF2arrayCpp(NumericMatrix climate_dt, IntegerVector site_ids, IntegerVector year_ids, bool include_month, int Nsites, int Nyears, int Nmonths, int Nenv, CharacterVector env_vars) {
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
      if (site_ids[j] == climate_dt(i, 1)) { // uniquePLOTid is the 2nd column (index 1)
        site = j;
        break;
      }
    }

    for (int j = 0; j < year_ids.size(); ++j) {
      if (year_ids[j] == climate_dt(i, 0)) { // year is the 1st column (index 0)
        year = j;
        break;
      }
    }

    if (site == -1 || year == -1) {
      stop("Site or year index not found. Check your data.");
    }

    for (int e = 0; e < Nenv; ++e) {
      int env_col = include_month ? 3 + e : 2 + e; // Adjusting column index for environmental variables
      if (include_month) {
        int month = climate_dt(i, 2) - 1; // month is the 3rd column (index 2)
        env_array(site + Nsites * (year + Nyears * (month + Nmonths * e))) = climate_dt(i, env_col);
      } else {
        env_array(site + Nsites * (year + Nyears * e)) = climate_dt(i, env_col);
      }
    }
  }

  return env_array;
}
