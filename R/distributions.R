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

binomial_from_gamma(torch_ones(c(2, 2), requires_grad = TRUE)+8, p = torch_ones(c(2,2))-0.5)



