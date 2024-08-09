#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector climateDF2arrayCpp(NumericMatrix climate_dt,
                                 IntegerVector site_ids,
                                 IntegerVector year_ids,
                                 bool include_month,
                                 bool include_day,
                                 int Nsites,
                                 int Nyears,
                                 int Nmonths,
                                 int Ndays,
                                 int Nenv,
                                 CharacterVector env_vars) {

  int total_rows = climate_dt.nrow();
  NumericVector env_array;

  if (include_day) {
    env_array = NumericVector(Nsites * Nyears * Nmonths * Ndays * Nenv, NA_REAL);
  } else if (include_month) {
    env_array = NumericVector(Nsites * Nyears * Nmonths * Nenv, NA_REAL);
  } else {
    env_array = NumericVector(Nsites * Nyears * Nenv, NA_REAL);
  }

  for (int i = 0; i < total_rows; ++i) {
    int site = -1;
    int year = -1;

    // Match the siteID
    for (int j = 0; j < site_ids.size(); ++j) {
      if (site_ids[j] == climate_dt(i, 0)) { // siteID is now the 1st column (index 0)
        site = j;
        break;
      }
    }

    // Match the year
    for (int j = 0; j < year_ids.size(); ++j) {
      if (year_ids[j] == climate_dt(i, 1)) { // year is now the 2nd column (index 1)
        year = j;
        break;
      }
    }

    if (site == -1 || year == -1) {
      stop("Site or year index not found. Check your data.");
    }

    for (int e = 0; e < Nenv; ++e) {
      int env_col;
      if (include_day) {
        int month = climate_dt(i, 2) - 1; // month is now the 3rd column (index 2)
        int day = climate_dt(i, 3) - 1;   // day is now the 4th column (index 3)
        env_col = 4 + e;                  // Environmental variables start from the 5th column (index 4)
        env_array(site + Nsites * (year + Nyears * (month + Nmonths * (day + Ndays * e)))) = climate_dt(i, env_col);
      } else if (include_month) {
        int month = climate_dt(i, 2) - 1; // month is now the 3rd column (index 2)
        env_col = 3 + e;                  // Environmental variables start from the 4th column (index 3)
        env_array(site + Nsites * (year + Nyears * (month + Nmonths * e))) = climate_dt(i, env_col);
      } else {
        env_col = 2 + e;                  // Environmental variables start from the 3rd column (index 2)
        env_array(site + Nsites * (year + Nyears * e)) = climate_dt(i, env_col);
      }
    }
  }

  return env_array;
}
