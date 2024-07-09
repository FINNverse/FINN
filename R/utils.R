#' Generate random numbers from a uniform distribution
#'
#' This function generates random numbers from a uniform distribution with specified low and high values and size similar to np.random.uniform in Python.
#'
#' @param low numeric Lower bound of the uniform distribution.
#' @param high numeric Upper bound of the uniform distribution.
#' @param size numeric Size of the output array.
#'
#' @return array A numeric array of random numbers.
#'
#' @examples
#' np_runif(0, 1, c(2, 3))
#'
#' @export
np_runif = function(low, high, size) {
  N = prod(size)  # Calculate total number of values to generate
  array(runif(N, low, high), dim = size)  # Generate random numbers and reshape into desired size
}

