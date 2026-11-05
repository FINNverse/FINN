#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

// Index of the first entry of `ids` per value, keyed on the value as a double
// (the climate matrix is numeric). Looking a row up in this map gives the same
// answer as the linear scan `for (j ...) if (ids[j] == x) break;` it replaces:
// first match, exact numeric equality (int promoted to double), and no match
// for NaN/NA (NaN never compares equal). The scan made every row O(Nsites), so
// a call was O(rows x sites) - ~90 s at 32,686 sites; the map makes it O(rows).
static std::unordered_map<double, int> first_index_map(const IntegerVector& ids) {
  std::unordered_map<double, int> m;
  m.reserve(ids.size() * 2 + 1);
  for (int j = 0; j < ids.size(); ++j) m.emplace(static_cast<double>(ids[j]), j);  // emplace keeps the first j
  return m;
}

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

  const std::unordered_map<double, int> site_map = first_index_map(site_ids);
  const std::unordered_map<double, int> year_map = first_index_map(year_ids);

  for (int i = 0; i < total_rows; ++i) {
    // Match the siteID (1st column, index 0) and the year (2nd column, index 1)
    std::unordered_map<double, int>::const_iterator s_it = site_map.find(climate_dt(i, 0));
    std::unordered_map<double, int>::const_iterator y_it = year_map.find(climate_dt(i, 1));

    if (s_it == site_map.end() || y_it == year_map.end()) {
      stop("Site or year index not found. Check your data.");
    }
    int site = s_it->second;
    int year = y_it->second;

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
