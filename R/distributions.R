#' Sample poisson relaxed
#'
#' @param lmbd lambda
#' @param num_samples number of samples
#' @param temperature temperature
#'
#' @export
sample_poisson_relaxed = function(lmbd, num_samples=50, temperature = 1e-2) {
  t = ((torch::torch_rand(c(num_samples,lmbd$shape))+0.0001)$log()$negative()/lmbd)$cumsum(dim=1L)
  relaxed_indicator = torch_sigmoid((1.0 - t) / temperature)
  N = relaxed_indicator$sum(1)
  return(N)
}


#' Sample from binomial with gradient
#'
#' @param n number of trials
#' @param p probability of success
#' @param sample_size sample size
#'
#' @export
binomial_from_gamma = function(n, p, sample_size=1) {
  mean = n * p
  variance = n * p * (1 - p)
  alpha = mean**2 / variance
  beta = mean / variance
  gamma_dist = torch::distr_gamma(alpha, beta)
  samples_gamma = gamma_dist$rsample(sample_size)$squeeze(1)
  return(samples_gamma)
}
