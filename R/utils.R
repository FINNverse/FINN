build_NN = function(input_shape,
                    output_shape,
                    hidden, # vector
                    bias, # vector
                    activation, # vector
                    dropout,
                    last_activation = "sigmoid") {

  model_list <- list()
  if (length(hidden) != length(activation)) {
    activation <- rep(activation[1], length(hidden))
  }

  if (length(bias) == 1) {
    bias <- rep(bias[1], length(hidden))
  }
  bias <- c(FALSE, bias)

  if (length(hidden) > 0) {
    for (i in 1:length(hidden)) {
      if (i == 1) {
        model_list <- c(model_list, torch::nn_linear(input_shape, hidden[i], bias = TRUE))
      } else {
        model_list <- c(model_list, torch::nn_linear(hidden[i - 1], hidden[i], bias = TRUE))
      }
      if (activation[i] == "relu") model_list <- c(model_list, torch::nn_relu())
      if (activation[i] == "selu") model_list <- c(model_list, torch::nn_selu())
      if (activation[i] == "leakyrelu") model_list <- c(model_list, torch::nn_leaky_relu())
      if (activation[i] == "tanh") model_list <- c(model_list, torch::nn_tanh())
      if (activation[i] == "sigmoid") model_list <- c(model_list, torch::nn_sigmoid())
      if (dropout > 0.0) model_list <- c(model_list, torch::nn_dropout(p = dropout))
    }
  }

  if (length(hidden) > 0) {
    model_list <- c(model_list, torch::nn_linear(hidden[length(hidden)], output_shape, bias = TRUE))
  } else {
    model_list <- c(model_list, torch::nn_linear(input_shape, output_shape, bias = FALSE))
  }
  if (last_activation == "sigmoid") model_list <- c(model_list, torch::nn_sigmoid())
  if (last_activation == "relu") model_list <- c(model_list, torch::nn_relu())
  return(do.call(torch::nn_sequential, model_list))
}




# The progress plot of fit(plot_progress = TRUE): the six loss terms per epoch,
# each as a fraction of its first-epoch value. The whole plot is redrawn every
# epoch with the y-axis fitted to the data, so a term that rises above its
# starting value is not cut off (the previous version fixed ylim at [0, 1] and
# only appended the last segment); a term that is switched off ("none") or has
# no data is NaN here and is simply not drawn.
visualize.training <- function(losses, epoch, max_epoch = 100L, main, new = FALSE, baseline = NULL){
  cols = c("#de324c", "#9656a2", "#95cf92", "#f4895f", "#000000", "#369acc")
  lab  = c("dbh", "ba", "trees", "growth", "mort", "reg")
  y = losses[, 1:6, drop = FALSE]
  y[!is.finite(y)] = NA
  shown = colSums(!is.na(y)) > 0
  top = suppressWarnings(max(y, na.rm = TRUE))
  if (!is.finite(top) || top <= 0) top = 1
  graphics::plot(c(), c(), xlim = c(1, max(max_epoch, 2)), ylim = c(0, top * 1.05),
                 main = main,
                 xlab = "epoch",
                 ylab = "loss relative to epoch 1",
                 las = 1)
  graphics::abline(h = 1, col = "grey80", lty = 2)
  if (any(shown))
    graphics::matplot(seq_len(nrow(y)), y[, shown, drop = FALSE], type = "l", lty = 1,
                      col = cols[shown], lwd = 1.5, add = TRUE)
  graphics::legend("topright", pch = 15, bty = "n", col = cols[shown],
                   legend = if (any(shown)) lab[shown] else character(0))
  invisible(NULL)
}




# Set torch's intra-op thread count for a fit (FINN#48). Returns the count to
# restore on exit, or NULL when there is nothing to restore. With OpenMP the
# count can change at any time; torch's native backend (macOS) accepts one change
# per session, so there the count is set but not restored, and a refused change
# is reported instead of silently running on another count.
set_torch_threads = function(threads) {
  if (is.null(threads)) return(NULL)
  threads = as.integer(threads)
  if (length(threads) != 1L || is.na(threads) || threads < 1L)
    stop("`threads` must be NULL or one positive integer.", call. = FALSE)
  old = torch::torch_get_num_threads()
  if (old == threads) return(NULL)
  torch::torch_set_num_threads(threads)
  now = torch::torch_get_num_threads()
  if (now != threads) {
    message(sprintf(paste0("torch keeps %d intra-op threads: this session's thread count was already set ",
                           "(torch's native backend allows one change per session). Restart R, or call ",
                           "torch::torch_set_num_threads(%d) before any other, to fit on %d."),
                    now, threads, threads))
    return(NULL)
  }
  if (torch::backends_openmp_is_available()) old else NULL
}

#' Index species
#'
#' @param pred predictions
#' @param species species index vector, must be int64
#'
#' @return A `torch` tensor of `pred` gathered along its species dimension by
#'   `species`.
#' @noRd
index_species = function(pred, species) {
  shapes = pred$shape
  X_expanded = pred$unsqueeze(2)$expand(c(shapes[1], species$shape[2], shapes[2]))
  pred = torch_gather(X_expanded, 3, species)
  return(pred)
}
# TODO: gradient necessary for indexing?! Test!



#' Sum cohort values per site and species
#'
#' The hot-path form of `aggregate_results()` (an internal function) used by the time loop of
#' `forward()`: every sample is scattered into a fresh species tensor and summed
#' over the patches in ONE `scatter_add` for all samples (stacked along a new
#' leading dimension), instead of one `torch_zeros` + `scatter_add` + `sum` +
#' `add_` per sample. The values are the same (a `scatter_add` accumulates each
#' output element over the cohorts in index order whether the samples are
#' stacked or not, and the sum over the patches is unchanged); only the number
#' of torch calls differs, which is what costs time on the small tensors of a
#' batch.
#'
#' @param labels torch.Tensor int64 species index, `[sites, patches, cohorts]`.
#' @param samples list of torch.Tensor, each `[sites, patches, cohorts]` and of
#'   the same dtype, the values to sum.
#' @param sp integer Number of species (the size of the last output dimension).
#'
#' @return A list of `torch` tensors `[sites, sp]`, one per element of `samples`,
#'   the per-site, per-species sums over patches and cohorts.
#' @keywords internal
#' @noRd
aggregate_species = function(labels, samples, sp) {
  k = length(samples)
  if(k == 1L) {
    s1 = samples[[1]]
    return(list(torch::torch_zeros(s1$shape[1], s1$shape[2], sp, dtype = s1$dtype, device = s1$device)$scatter_add(3, labels, s1)$sum(dim = 2)))
  }
  st = torch::torch_stack(samples)                                          # [k, sites, patches, cohorts]
  idx = labels$unsqueeze(1)$expand(c(k, -1, -1, -1))
  out = torch::torch_zeros(k, st$shape[2], st$shape[3], sp, dtype = st$dtype, device = st$device)$scatter_add(4, idx, st)$sum(dim = 3)
  lapply(seq_len(k), function(v) out[v, , ])
}


# Per-site, per-species sums (or means) of cohort values: every sample is
# scattered into a species tensor and reduced over the patches. The time loop
# of forward() uses the stacked form aggregate_species(); this one serves the
# adult-density terms of the regeneration processes and the per-tree light.
aggregate_results = function(labels, samples, Results, drop_rows = TRUE, sp_max = NULL, aggregation = "sum") {

  if(is.null(sp_max)) sp_max = as.numeric(Results[[1]]$shape[2])

  if(length(aggregation) == 1) aggregation = rep(aggregation, length(samples))

  for( v in 1:length(samples)) {
    if(aggregation[v] == "sum") {
    Results[[v]]$add_(torch::torch_zeros(samples[[v]]$shape[1], samples[[v]]$shape[2], sp_max, dtype=samples[[v]]$dtype, device = samples[[v]]$device)$scatter_add(3, labels, samples[[v]]$expand(c(samples[[v]]$shape[1], samples[[v]]$shape[2], -1)))$sum(dim = 2))
    } else {
      Results[[v]]$add_(torch::torch_zeros(samples[[v]]$shape[1], samples[[v]]$shape[2], sp_max, dtype=samples[[v]]$dtype, device = samples[[v]]$device)$scatter_add(3, labels, samples[[v]]$expand(c(samples[[v]]$shape[1], samples[[v]]$shape[2], -1)))$mean(dim = 2))
      }

  }
  return(Results)

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
  # floor() loses about 0.5 recruits per draw above lambda 2; round() is planned
  # for 0.2.1, where it changes results (see ?rnbinom_torch_refined)
  poisson_samples_large_lambda = gaussian_samples$clamp(min = 0)$floor()

  # Interpolate between the two approaches
  poisson_samples = smooth_transition * poisson_samples_large_lambda + (1 - smooth_transition) * poisson_samples_small_lambda

  return(poisson_samples)
}

# `mortality_floor` / `mortality_ceiling` of createProcess()/createHybrid(): one probability in (0, 1).
check_mortality_limit = function(value, what) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0 || value >= 1)
    stop("`", what, "` must be one number in (0, 1).", call. = FALSE)
  as.numeric(value)
}

#' Draw deaths from a Gamma matched to the binomial
#'
#' @description
#' The mortality draw of FINN before 0.2.0, selected with
#' [`createProcess(mortality_draw = "gamma")`][createProcess]: a Gamma distribution with the mean
#' `n * p` and the variance `n * p * (1 - p)` of the binomial, sampled with
#' `rsample()` so that the draw carries a pathwise gradient to `n` and `p`. The
#' draw is continuous, not an integer count; the simulation rounds it afterwards.
#' See [binomial_from_poisson()] for what that rounding costs and
#' [FINN::createProcess()] for the other draws.
#'
#' @param n torch.Tensor Number of trees at risk per cohort.
#' @param p torch.Tensor Mortality probability per cohort, strictly between 0
#'   and 1 (the Gamma's shape `n * p / (1 - p)` degenerates at both ends).
#' @param sample_size integer Number of draws; only the default `1` is used by
#'   FINN (the leading sample dimension is dropped).
#'
#' @return A `torch` tensor of the shape of `n` with continuous, positive draws
#'   whose mean is `n * p` and whose variance is `n * p * (1 - p)`.
#' @seealso [binomial_from_poisson()], [binomial_from_bernoulli()],
#'   [FINN::createProcess()] (`mortality_draw`)
#' @examplesIf torch::torch_is_installed()
#' FINN.seed(1)
#' n <- torch::torch_full(c(10000), 40)   # 10 000 cohorts of 40 trees
#' p <- torch::torch_full(c(10000), 0.3)
#' d <- as.numeric(binomial_from_gamma(n, p))
#' # the binomial has mean 12 and variance 8.4; the draws are not integers
#' c(mean = mean(d), var = var(d))
#' head(d, 3)
#' @keywords internal
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


#' Draw deaths as a Bernoulli/Poisson count
#'
#' @description
#' The alternative mortality draw, selected with
#' [`createProcess(mortality_draw = "poisson")`][createProcess]. The default since 0.2.0 is
#' `mortality_draw = "binomial"`, the exact binomial; [binomial_from_gamma()]
#' (`mortality_draw = "gamma"`) reproduces results published before 0.2.0.
#' A cohort of one tree or fewer (`n <= 1.5`) is a Bernoulli trial and is
#' drawn as one; every larger cohort is drawn from a Poisson capped at `n`.
#'
#' @details
#' [binomial_from_gamma()] matches a Gamma to the binomial's first two moments and
#' `rsample()`s it, which buys a pathwise gradient -- the hazard it guards against
#' is a discrete draw having no gradient at all. The forward pass then rounds that
#' draw straight-through (`finn-class.R`), so the gradient flows through a quantity
#' the simulation discards, and the rounded count delivers a measured 0.68-0.72 of
#' the calibrated deaths for one-tree cohorts at annual probabilities of 0.01-0.10
#' (5e5 draws per probability; 0.76 at `p = 0.26`). This draw is integer-valued, exact
#' in expectation for small `n*p`, and carries no pathwise gradient, which is how
#' Lin et al. (2026) treat discrete events in a differentiable crop model: resolve
#' them before differentiation and hold them fixed during gradient evaluation.
#' The Bernoulli branch is not cosmetic -- capping a Poisson at a one-tree cohort
#' delivers only 88\% of the calibrated probability at the `p = 0.26` of a
#' suppressed stem, which is the stem that decides thinning.
#'
#' The choice shows up in the simulated stand rather than in a holdout score: a
#' long simulation to equilibrium separates the two draws by a factor of several
#' in basal area and stem number, while the fitted loss moves far less, because
#' the rounded deficit is a bias on the trajectory and not on the parameters.
#'
#' @param n torch.Tensor Number of trees at risk per cohort.
#' @param p torch.Tensor Annual mortality probability per cohort.
#'
#' @return A `torch` tensor of integer death counts, never exceeding `n`.
#' @references
#' Lin et al. (2026). \doi{10.48550/arXiv.2610.36827}
#' @seealso [binomial_from_gamma()], [binomial_from_bernoulli()],
#'   [FINN::createProcess()] (`mortality_draw`)
#' @examplesIf torch::torch_is_installed()
#' FINN.seed(1)
#' # cohorts of one tree are Bernoulli trials, larger cohorts a capped Poisson
#' n <- torch::torch_tensor(rep(c(1, 40), each = 10000))
#' p <- torch::torch_full(c(20000), 0.1)
#' d <- as.numeric(binomial_from_poisson(n, p))
#' # mean deaths per cohort: expected 0.1 (one tree) and 4 (forty trees)
#' tapply(d, rep(c("n = 1", "n = 40"), each = 10000), mean)
#' range(d)   # whole trees, and never more than the cohort holds
#' @keywords internal
#' @export
binomial_from_poisson = function(n, p) {
  one = n$le(1.5)$to(dtype = n$dtype)
  one * torch::torch_bernoulli(p) +
    (1 - one) * torch::torch_minimum(torch::torch_poisson(n * p), n)
}


#' Draw binomial counts from per-trial Bernoulli probabilities
#'
#' @description
#' The mortality draw selected with `createProcess(mortality_draw =
#' "bernoulli")`: every element of `n` gets `max(n)` reparameterised Bernoulli
#' trials with probability `p`, of which the first `n` are kept and summed. The
#' result is integer-valued and matches the binomial in mean and variance, and
#' it keeps a gradient path to `p` (a straight-through sigmoid) and to `n` (a
#' soft mask). The price is memory: the trials are laid out along an extra
#' dimension of length `max(n)`, so the cost grows with the largest cohort.
#'
#' @param n torch.Tensor Number of Bernoulli trials (trees at risk) per
#'   element, with at least two dimensions (FINN passes
#'   `[sites, patches, cohorts]`) and a largest value of at least 1.
#' @param p torch.Tensor Success (death) probability per element, of the shape
#'   of `n`.
#' @return A `torch` tensor of the shape of `n`: the number of successes per
#'   element, a whole number between 0 and `n`.
#' @seealso [binomial_from_gamma()], [binomial_from_poisson()],
#'   [FINN::createProcess()] (`mortality_draw`)
#' @examplesIf torch::torch_is_installed()
#' FINN.seed(1)
#' # one site, 5000 patches, two cohorts of 5 and 40 trees, 30 % mortality
#' n <- torch::torch_tensor(array(rep(c(5, 40), each = 5000), dim = c(1, 5000, 2)))
#' p <- torch::torch_full(c(1, 5000, 2), 0.3)
#' deaths <- torch::as_array(binomial_from_bernoulli(n, p))
#' dim(deaths)
#' # binomial mean n * p = 1.5 and 12, variance n * p * (1 - p) = 1.05 and 8.4
#' sapply(c("n = 5" = 1, "n = 40" = 2),
#'        function(k) c(mean = mean(deaths[1, , k]), var = var(deaths[1, , k])))
#' @export
binomial_from_bernoulli = function(n, p) {

  max_trials = n$max()$item()
  samples_shape = c(n$shape, max_trials)

  u = torch_rand(samples_shape)
  logistic_samples = torch_log(u / (1 - u))
  logits = torch_log(p / (1 - p))$unsqueeze(-1)
  scaled_samples = logits + logistic_samples

  # THRESHOLD the reparameterised draw; do not feed it to a second torch_bernoulli().
  # With L ~ Logistic(0,1), P(logit(p) + L > 0) = p exactly, which is the trial. The
  # old code drew a second Bernoulli with probability sigmoid(logit(p) + L), whose
  # mean is not p -- measured 0.181 against 0.100 at p = 0.1 (FINNverse/FINN#62).
  # The sigmoid surrogate keeps a gradient path to `p`.
  hard = (scaled_samples > 0)$to(dtype = scaled_samples$dtype)
  soft = torch_sigmoid(scaled_samples / 1e-3)
  bernoulli_samples = soft + (hard - soft)$detach()
  # Keep trial k while k <= n. The sign used to be the other way round, which kept
  # exactly the trials that should be dropped and made the sum identically zero for
  # every input (FINNverse/FINN#62). The sigmoid rather than a hard `<=` is what lets
  # a gradient reach `n`.
  mask = torch_sigmoid((n$unsqueeze(-1) + 0.1 - torch_arange(1, max_trials)$unsqueeze(1)$unsqueeze(1)$expand_as(bernoulli_samples))/1e-3)

  masked_samples = bernoulli_samples * mask

  binomial_samples = masked_samples$sum(dim=-1)
  return(binomial_samples)
}





#' Set Seed for Reproducibility in R and Torch
#'
#' This function sets the seed for both R's random number generator and Torch's random number generator,
#' ensuring reproducibility across operations that involve both R and Torch.
#'
#' @param seed An integer value to set as the seed for both R and Torch. This ensures that random operations
#'     in both environments produce consistent results.
#'
#' @return This function does not return a value. It sets the seed internally for both R and Torch.
#'
#' @details The function calls `set.seed()` to set the seed for R's random number generator and
#'     `torch::torch_manual_seed()` to set the seed for Torch's random number generator. This is useful
#'     for ensuring reproducibility in scripts that rely on both R and Torch for random operations.
#'
#' @seealso [base::set.seed()], [torch::torch_manual_seed()]
#'
#' @examplesIf torch::torch_is_installed()
#' # the same seed gives the same draws, in R and in torch
#' FINN.seed(123)
#' a <- c(runif(1), as.numeric(torch::torch_rand(1)))
#' FINN.seed(123)
#' b <- c(runif(1), as.numeric(torch::torch_rand(1)))
#' identical(a, b)
#'
#' # seed before finn() and fit() to make a fit reproducible: the starting
#' # values of the parameters and the simulated deaths and recruits are random
#'
#' @import torch
#' @export
FINN.seed <- function(seed) {
  set.seed(seed)
  torch::torch_manual_seed(seed)
}


# default ranges
default_speciesPars_ranges = list(
  parGrowth = rbind(
    c(0.01, 0.99),
    c(0.0, 1)
  ),
  parMort = rbind(
    c(0.01, 0.99),
    c(0, 4)
  ),
  parReg = c(0.01, 0.99),
  parHeight = c(0.3, 0.7)
)



rnbinom_torch = function(pred, theta) {
  # cat("\nPred: ")
  # print(pred$min())
  # print(pred$max())
  # cat("\nTheta: ")
  # print(theta$min())
  #pred = pred$abs()
  #theta = theta$abs()
  # clamp away from exact 0/1, mirroring dnbinom_torch: pred -> 0 (e.g. a
  # near-zero predicted regeneration rate) otherwise sends probs -> 0 exactly,
  # so log(probs) -> -Inf and the resulting 0*Inf term in the backward pass
  # through distr_gamma()$rsample() shows up as NaN gradients.
  eps = 0.0001
  probs = torch::torch_clamp(1.0 - theta/(theta + pred) + eps, 0.0, 1.0 - eps)
  total_count = theta
  logits = torch::torch_log(probs) - torch::torch_log1p(-probs)
  return(sample_poisson_gaussian((torch::distr_gamma(total_count, (-logits)$exp() )$rsample() +0.001)$abs() ))
  # stats::rgamma(length(pred),total_count, exp(- logits ))
}

#' Negative-binomial recruit draw without the recruitment floor
#'
#' @description
#' The default regeneration sampler since 0.2.0, selected with
#' [`createProcess(regeneration_sampler = "refined")`][createProcess]; `rnbinom_torch()`
#' (`regeneration_sampler = "gaussian"`) reproduces results published before
#' 0.2.0. Three differences: the probability clamp is `1e-10` instead of `1e-4`,
#' the Gamma draw is clamped at that same `eps` instead of taking `+0.001` and an
#' absolute value, and the draw goes through `sample_poisson_refined()` (an exact
#' `torch_poisson` below `lambda = 2`, the Gaussian approximation above) instead of
#' `sample_poisson_gaussian()`.
#'
#' @details
#' The `+0.001` offset and the `1e-4` clamp exist for the Gaussian branch: they
#' keep the Gamma draw off zero and stop `log(0)` sending a NaN gradient back
#' through `distr_gamma()$rsample()`. Neither hazard applies to `torch_poisson`,
#' and together they impose a recruitment FLOOR of 0.0011 per cell-year, which at a
#' FINN patch of 0.0168 ha is 0.066 recruits ha^-1 yr^-1 *per species* -- over a
#' pool of a hundred-odd species, as large as the whole ingrowth a temperate
#' inventory records. `sample_poisson_gaussian()`
#' (`lambda + eps*sqrt(lambda)`, rounded straight-through by the forward pass)
#' rounds the floor away but also rounds the signal away: at a FINN patch of
#' 0.0168 ha the per-cell mean is 0.001-0.008, and over 2e6 draws the delivered
#' fraction of the calibrated mean was 0.000-0.048 at `theta = 1`.
#'
#' One bias is left and is not repaired in 0.2.0: above `lambda = 2` the draw is
#' the Gaussian approximation FLOORED, so it delivers about half a recruit per
#' draw too few (measured over 1e6 draws: -0.43 at `lambda = 2.5`, -0.50 from
#' `lambda = 20`, against about 0 with rounding). It bites only where a cell's
#' mean recruitment is above two, which the per-patch means above are not.
#' Rounding instead of flooring is planned for 0.2.1; it changes results, so it
#' is not in this release.
#'
#' For the measured consequence of both repairs on the simulated stand, and why
#' the holdout cannot separate them, see [binomial_from_poisson()].
#'
#' @param pred torch.Tensor Expected recruits per patch.
#' @param theta torch.Tensor Negative-binomial dispersion.
#'
#' @return A `torch` tensor of recruit counts, of the shape of `pred`.
#' @seealso [recruitment], [FINN::createProcess()] (`regeneration_sampler`,
#'   `theta_floor`)
#' @examplesIf torch::torch_is_installed()
#' FINN.seed(1)
#' # sparse recruitment: an expected 0.05 recruits per patch and year
#' mu    <- torch::torch_full(c(100000), 0.05)
#' theta <- torch::torch_tensor(1)
#' r <- as.numeric(rnbinom_torch_refined(mu, theta))
#' mean(r)    # the draw delivers the expected count ...
#' table(r)   # ... as whole recruits, most patches getting none
#' @keywords internal
#' @export
rnbinom_torch_refined = function(pred, theta) {
  eps = 1e-10
  probs = torch::torch_clamp(1.0 - theta/(theta + pred) + eps, eps, 1.0 - eps)
  logits = torch::torch_log(probs) - torch::torch_log1p(-probs)
  sample_poisson_refined(torch::distr_gamma(theta, (-logits)$exp())$rsample()$clamp(min = eps))
}

dnbinom_torch = function(pred, true, theta) {
  eps = 0.0001
  #pred = torch::torch_exp(pred)
  theta_tmp = theta#$abs() # 1.0/(torch::nnf_softplus(theta)+eps)
  probs = torch::torch_clamp(1.0 - theta_tmp/(theta_tmp+pred)+eps, 0.0, 1.0-eps)
  total_count = theta_tmp
  value = true
  logits = torch::torch_log(probs) - torch::torch_log1p(-probs)
  log_unnormalized_prob <- total_count * torch::torch_log(torch::torch_sigmoid(-logits)) + value * torch::torch_log(torch::torch_sigmoid(logits))
  log_normalization <- -torch::torch_lgamma(total_count + value) + torch::torch_lgamma(1 + value) + torch::torch_lgamma(total_count)
  log_normalization <- torch::torch_where(total_count + value == 0, torch::torch_tensor(0, dtype = log_normalization$dtype, device = pred$device), log_normalization)
  return( - (log_unnormalized_prob - log_normalization))
  }



# The simulated yearly rates of an observation window, compounded into the
# per-year rate makeObsData(rates = "annual") reports: growth as the geometric
# mean of 1 + g over the window, mortality as 1 minus the geometric mean
# survival. `x` is [sites, years, species] over the years 1..i (or the window),
# `w` an optional [sites, years, 1] mask of the years inside each site's window,
# `p` the window length per site ([sites, 1]) or one number. A year the mask
# excludes, and a year in which a species has no trees (its rate is 0),
# contributes growth 0 and survival 1 while the denominator stays the full
# period_length: the observed rate is derived from the counts over the whole
# interval, so the simulated one is spread over the whole interval too.
# Mortality is clamped below 1 - 1e-4 before the log: sigmoid saturates to
# exactly 1 in float32 (logit > 17.3), where log(0) would make the gradient
# infinite; 1e-4 bounds it at 1e4. Growth is not clamped: the package's growth
# is non-negative, and a custom process that goes below -1 should fail visibly.
# The input is cloned first: log1p() and clamp() save their input for the
# backward pass, and the caller passes a view of the Result tensors, which the
# time loop keeps writing into before the per-site path runs its one backward.
window_growth = function(g, w = NULL, p) {
  lg = torch_log1p(g$clone())
  if (!is.null(w)) lg = lg*w
  torch_expm1(lg$sum(2)/p)
}
window_mortality = function(m, w = NULL, p) {
  ls = torch_log1p(-m$clone()$clamp(max = 1 - 1e-4))
  if (!is.null(w)) ls = ls*w
  1 - torch_exp(ls$sum(2)/p)
}

# `"cuda"` is what everyone who has used torch directly reaches for, and
# match.arg()'s "should be one of cpu, gpu" does not say that it is accepted.
resolve_device = function(device) {
  if (is.character(device) && length(device) > 1L) device = device[1]
  if (identical(device, "cuda")) device = "gpu"
  match.arg(device, c("cpu", "gpu"))
}

# torch's binomial is exact in both moments at constant memory but is not exported
# from torch and has no derivative: "derivative for aten::binomial is not
# implemented". Only reachable where the draw is detached anyway (rb_gradient).
torch_binomial_exact = function(n, p) {
  getFromNamespace("torch_binomial", "torch")(n, p$detach())
}
