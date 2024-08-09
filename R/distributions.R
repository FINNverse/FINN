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

sample_poisson_gaussian <- function(lmbd) {
  epsilon <- torch_randn_like(lmbd)
  poisson_samples <- lmbd + epsilon * lmbd$sqrt()
  poisson_samples <- poisson_samples$clamp(min = 0)
  return(poisson_samples)
}

sample_poisson_refined <- function(lmbd, temperature = 1e-2) {
  # Ensure lambda is a tensor

  # Generate Gaussian noise
  epsilon <- torch_randn_like(lmbd)

  # Compute Gaussian approximation
  gaussian_samples <- lmbd + epsilon * lmbd$sqrt()

  # Apply sigmoid function to smooth transition between small and large lambda
  smooth_transition = torch_sigmoid((lmbd - 2) / temperature)

  # Use Poisson approximation for small lambda and Gaussian approximation for large lambda
  poisson_samples_small_lambda = torch_poisson(lmbd)
  poisson_samples_large_lambda = gaussian_samples$clamp(min = 0)$floor()

  # Interpolate between the two approaches
  poisson_samples = smooth_transition * poisson_samples_large_lambda + (1 - smooth_transition) * poisson_samples_small_lambda

  return(poisson_samples)
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


#' @export
binomial_from_bernoulli = function(n, p) {

  max_trials = n$max()$item()
  samples_shape = c(n$shape, max_trials)

  u = torch_rand(samples_shape)
  logistic_samples = torch_log(u / (1 - u))
  logits = torch_log(p / (1 - p))$unsqueeze(-1)
  scaled_samples = logits + logistic_samples
  probabilities = torch_sigmoid(scaled_samples)
  #binomial_samples = torch_sum(torch_bernoulli(probabilities$expand(c(n, p$shape))), dim=1)

  # Generate binomial samples by summing up to the number of trials for each element
  bernoulli_samples = torch_bernoulli(probabilities)
  #mask = torch_arange(1, max_trials)$unsqueeze(1)$unsqueeze(1)$expand_as(bernoulli_samples) < n$unsqueeze(-1)
  mask = torch_sigmoid((torch_arange(1, max_trials)$unsqueeze(1)$unsqueeze(1)$expand_as(bernoulli_samples) - n$unsqueeze(-1) - 0.1)/1e-3 )

  masked_samples = bernoulli_samples * mask

  binomial_samples = masked_samples$sum(dim=-1)
  return(binomial_samples)
}



sample_poisson_gumbel <- function(lmbd, num_samples=50, temperature = 1e-2) {
  shape <- lmbd$shape
  gumbels <- -torch::torch_log(-torch::torch_log(torch::torch_rand(c(num_samples, shape))))
  logits <- lmbd$log()
  poisson_samples <- (logits + gumbels) / temperature
  poisson_samples <- poisson_samples$exp()$cumsum(dim=1L)
  relaxed_indicator <- torch_sigmoid((1.0 - poisson_samples) / temperature)
  N <- relaxed_indicator$sum(1)
  return(N)
}



