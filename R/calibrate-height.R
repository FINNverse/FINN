#' Calibration of Tree Heights
#'
#' This script performs the calibration of tree heights by optimizing parameters for different species.
#' The optimization is done using the `optim` function, and the results include parameter estimates,
#' standard errors, t-values, p-values, and confidence intervals.
#'
#' @param treevars_dt A `data.table` containing tree variables, including species, DBH, and actual heights.
#' @return A `data.table` with the optimized parameters and statistics for each species.
#' @examples
#' # Assuming treevars_dt is already defined with the necessary columns
#' results <- calibrate_height(treevars_dt)
#' print(results)
#' @export

library(data.table)

# Function to calculate standard errors, t-values, p-values, and confidence intervals
#' @param hessian The Hessian matrix obtained from the optimization.
#' @param estimate The estimated parameter value.
#' @return A list containing standard errors, t-values, p-values, and confidence intervals.
calc_ci <- function(hessian, estimate) {
  var_cov_matrix <- solve(hessian)
  se <- sqrt(diag(var_cov_matrix))
  t_values <- estimate / se
  p_values <- 2 * (1 - pnorm(abs(t_values)))
  ci_lower <- estimate - 1.96 * se
  ci_upper <- estimate + 1.96 * se
  list(se = se, t_values = t_values, p_values = p_values, CI_lower = ci_lower, CI_upper = ci_upper)
}

# Function to calculate the sum of squared residuals for height optimization
#' @param par The parameter to optimize.
#' @param data The data containing DBH and actual heights.
#' @return The sum of squared residuals.
optimHeight <- function(par, data) {
  # Calculate the predicted height
  predicted_height <- FINN::height(data$DBH, par)

  # Calculate the residuals
  residuals <- data$ACTUALHTm - predicted_height

  # Calculate the sum of squared residuals
  sum_sq_residuals <- sum(residuals^2, na.rm = T)

  # Return the sum of squared residuals
  sum_sq_residuals
}

# Function to calibrate tree height parameters for different species
#' @param treevars_dt A `data.table` containing tree variables, including species, DBH, and actual heights.
#' @return A `data.table` with the optimized parameters and statistics for each species.
calibrate_height <- function(treevars_dt){
  parHeight_dt <- data.table()
  species_list <- unique(treevars_dt$species)
  total_species <- length(species_list)
  start_time <- Sys.time()

  # Iterate over each species
  for(i in seq_along(species_list)) {
    i_species <- species_list[i]
    data <- treevars_dt[species == i_species]

    # Perform optimization with error handling
    result <- tryCatch({
      out <- optim(par = 0.5, fn = optimHeight, data = data, method = "Brent", lower = 0, upper = 1, hessian = TRUE)
      ci_results <- calc_ci(out$hessian, out$par)
      list(par = out$par, ci_results = ci_results, error = "Success")
    }, error = function(e) {
      list(par = NA, ci_results = list(se = NA, t_values = NA, p_values = NA, CI_lower = NA, CI_upper = NA), error = e$message)
    })

    # Append results to the data.table
    parHeight_dt <- rbind(parHeight_dt, data.table(
      species = i_species,
      Estimate = result$par,
      StdError = result$ci_results$se,
      tValue = result$ci_results$t_values,
      pValue = result$ci_results$p_values,
      CI_lower = result$ci_results$CI_lower,
      CI_upper = result$ci_results$CI_upper,
      status = result$error,
      N = nrow(data[!is.na(ACTUALHTm)])
    ))

    # Calculate elapsed time and estimate remaining time
    elapsed_time <- Sys.time() - start_time
    avg_time_per_iteration <- elapsed_time / i
    remaining_time <- avg_time_per_iteration * (total_species - i)

    # Overwrite previous console output with estimated remaining time
    cat("\r", sprintf("%d of %d - %s - Elapsed time: %.2f mins, Estimated remaining time: %.2f mins",
                      i, total_species, result$error, as.numeric(elapsed_time, units = "mins"), as.numeric(remaining_time, units = "mins")), sep = "")
  }

  cat("\n")
  return(parHeight_dt)
}
