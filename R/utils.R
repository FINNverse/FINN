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


check_and_recreate = function(tensor, r_obj, dtype=torch::torch_float32(), device="cpu", requires_grad = FALSE) {
  pointer_check <- tryCatch(torch::as_array(tensor), error = function(e) e)
  if(!is.null(attr(r_obj, "requires_grad"))) requires_grad = attr(r_obj, "requires_grad")
  if(inherits(pointer_check,"error")){
    tensor = torch::torch_tensor(r_obj, dtype=dtype, device=device, requires_grad = requires_grad)
  }
  return(tensor)
}


#' Index species
#'
#' @param pred predictions
#' @param species species index vector, must be int64
#'
#' @export
index_species = function(pred, species) {
  shapes = pred$shape
  X_expanded = pred$unsqueeze(2)$expand(c(shapes[1], species$shape[2], shapes[2]))
  pred = torch_gather(X_expanded, 3, species)
  return(pred)
}
# TODO: gradient necessary for indexing?! Test!



#' Aggregate function
#'
#' @param labels labels
#' @param samples samples
#' @param Results Results
#'
#' @export
aggregate_results_old = function(labels, samples, Results, drop_rows = TRUE, sp_max = NULL) {

  if(drop_rows) {
    # Drop rows for speed up - memory intensive
    if(is.null(sp_max)) sp_max = as.numeric(Results[[1]]$shape[2]) # Not good, better to pass sp_max to function!

    old_shape = labels$shape

    labels_new = (labels + torch::torch_tensor(seq(1, by = sp_max, length.out = labels$shape[1])-1,
                                               dtype = torch_int64(),
                                               device = labels$device)$reshape(c(labels$shape[1], 1L, 1L))$expand(labels$shape))$reshape(c(1, labels$shape[1]*labels$shape[2], labels$shape[3]))

    samples_new = lapply(samples, function(sample) {
      return(sample$reshape(c(1, labels$shape[1]*labels$shape[2], labels$shape[3])))
    })
    Results_new = lapply(Results, function(result) {
      return(result$reshape(c(1, labels$shape[1]*sp_max)))
    })
    Results = Results_new
    labels = labels_new
    samples = samples_new
  }

  for(k in 1:labels$shape[1]) {
    #samples_tmp = [sample[k,,,drop=FALSE]$flatten()$view(list(-1,1)) for sample in samples]
    samples_tmp = lapply(samples, function(sample) sample[k,,,drop=FALSE]$flatten()$view(list(-1,1)))

    results = groupby_mean(samples_tmp, labels[k,,,drop=FALSE]$flatten())
    values = results[[1]]
    positions = results[[2]]

    for( v in 1:length(samples)) {
      Results[[v]][k,positions] = Results[[v]][k,positions] + values[[v]]$squeeze()
    }
  }

  if(drop_rows) {
    Results = lapply(Results, function(result) {
      return(result$reshape(c(old_shape[1], sp_max)))
    })
  }

  return(Results)
}



aggregate_results = function(labels, samples, Results, drop_rows = TRUE, sp_max = NULL) {

  # if(is.null(sp_max)) sp_max = as.numeric(Results[[1]]$shape[2])
  # old_shape = labels$shape
  # labels = labels$flatten(start_dim=1, end_dim = 2)
  #
  # for( v in 1:length(samples)) {
  #   da = samples[[v]]$flatten(start_dim=1, end_dim=2)
  #   Results[[v]]$add_(torch::torch_zeros(da$shape[1], sp_max, dtype=da$dtype, device = da$device)$scatter_add(2, labels, da$expand(c(da$shape[1], -1)))$view(c(old_shape[1:2], sp_max))$sum(dim = 2))
  # }
  # return(Results)
  if(is.null(sp_max)) sp_max = as.numeric(Results[[1]]$shape[2])

  for( v in 1:length(samples)) {
    Results[[v]]$add_(torch::torch_zeros(samples[[v]]$shape[1], samples[[v]]$shape[2], sp_max, dtype=samples[[v]]$dtype, device = samples[[v]]$device)$scatter_add(3, labels, samples[[v]]$expand(c(samples[[v]]$shape[1], samples[[v]]$shape[2], -1)))$sum(dim = 2))
  }
  return(Results)

}

#' group by mean
#'
#' @param values list of value tensors
#' @param labels labels
#'
#' @export
groupby_mean = function(values, labels) {
  uniques = unique(as.matrix(labels))
  key_val = 1:length(uniques)
  names(key_val) = val_key = sort(uniques)
  names(val_key) = 1:length(uniques)
  labels2 = as.matrix(labels)[,1]
  labels = torch_tensor(key_val[as.character(labels2)],  dtype = torch_int64(), device = values[[1]]$device)
  labels = labels$view(c(labels$size(1), 1))$expand(c(-1, values[[1]]$size(2)))
  unique_labels_r = matrix(sort(unique(as.matrix(labels))), ncol = 1L)
  unique_labels = torch_tensor(unique_labels_r, dtype = torch_int64())
  results = lapply(values, function(value) torch_zeros_like(unique_labels, dtype = value$dtype, device = values[[1]]$device)$scatter_add_(1, labels, value))
  new_position = torch_tensor(val_key[as.character(unique_labels_r[,1])], device = values[[1]]$device,  dtype = torch_int64())
  return(list(results, new_position))
}


#' Padding tensors to remove dead cohorts

pad_tensors_speed_up = function(value, indices, org_dim, const = 0) {
  # KK = torch::torch_split(value$flatten(start_dim = 1, end_dim = 2), split_size = org_dim[1]*org_dim[2], dim = 1) TODO does not work...bug?!
  KK = value$flatten(start_dim = 1, end_dim = 2)$split(1, 1)
  KK_new = vector("list", length = length(KK))
  for( i in 1:indices$shape[1]) {
    KK_new[[i]] = KK[[i]][1, indices[i, ]$detach()]
  }
  return(torch::nn_utils_rnn_pad_sequence(KK_new, batch_first=TRUE, padding_value = const))
}

# pad_tensors_speed_up2 <- function(value, indices, org_dim, const = 0) {
#   # Flatten the value tensor to 2D
#   flat_value <- value$flatten(start_dim = 1, end_dim = 2)
#
#   # Use the boolean indices to filter out the alive cohorts
#   alive_cohorts <- flat_value$masked_select(indices$unsqueeze(1))$view(c(indices$sum(dim = 2)$max()$item(), -1))
#
#   # Pad the sequences
#   padded_tensor <- torch::nn_utils_rnn_pad_sequence(list(alive_cohorts), batch_first = TRUE, padding_value = const)
#
#   return(padded_tensor)
# }
#
# pad_tensors_speed_up3 <- function(value, indices, org_dim, const = 0) {
#   # Flatten the value tensor to 2D
#   flat_value <- value$flatten(start_dim = 1, end_dim = 2)
#
#   # Get the shape parameters
#   sites_patches <- indices$shape[1]
#   cohorts <- indices$shape[2]
#
#   # Calculate the maximum number of alive cohorts
#   max_alive <- indices$sum(dim = 2)$max()$item()
#
#   # Create a tensor to store the padded cohorts
#   padded_tensor <- torch::torch_full(c(sites_patches, max_alive), const)
#
#   # Expand the indices tensor for broadcasting
#   indices_expanded <- indices$unsqueeze(1)
#
#   # Gather the alive cohorts
#   alive_cohorts <- flat_value$masked_select(indices_expanded)
#
#   # Reshape the alive cohorts tensor
#   alive_cohorts_reshaped <- alive_cohorts$view(c(sites_patches, -1))
#
#   # Fill the padded tensor with alive cohorts
#   padded_tensor[, 1:alive_cohorts_reshaped$shape[2]] <- alive_cohorts_reshaped
#
#   return(padded_tensor)
# }
# value <- torch_tensor(array(sample(0:5, 10*50*9, replace = TRUE), dim = c(10, 50, 9)), requires_grad = TRUE, dtype = torch_float32()) # [sites, patches, cohorts]
# indices = (value > 0)$flatten(start_dim = 1, end_dim = 2) # boolean [sites*patches, cohorts], reporting which cohorts are alive or dead
# org_dim <- c(10, 50, 9)
#
# padded_tensor <- pad_tensors_speed_up(value, indices, org_dim) # output [sites*patches, mininum number of cohort dimension required (usually == highest number of alive cohorts in all sites and patches)]
#
# padd
#
# print(padded_tensor)
# # in this example: [5,5]
#
# indices = (M>0)$flatten(start_dim = 1, end_dim = 2)
#
# as_array(pad_tensors_speed_up(M, indices, org_dim = c(1, 10, 5)))
#
# M = matrix(sample(c(0, 1), 50, replace = TRUE), 10, 5)
# M
#
# remove_zeros_and_pad(M[1,,])
#
#
#
# # Print the result tensor
# print(result_tensor)
#
#
# remove_zeros_and_pad = function(tensor, index) {
#   mask = index
#   non_zero_counts = mask$sum(dim=2)
#   max_non_zeros = non_zero_counts$max()$item()
#   sorted_tensor = torch_sort(mask$float(), dim=2, descending=TRUE)[[2]]
#   sorted_tensor = tensor$gather(2, sorted_tensor)[, 1:max_non_zeros]
#
#   return(sorted_tensor)
# }
#
# pad_tensors_speed_up(value, indices, org_dim)
# tensor = value$flatten(start_dim = 1, end_dim = 2)
# remove_zeros_and_pad(tensor, indices) == pad_tensors_speed_up(value, indices, org_dim)
#
#
# microbenchmark::microbenchmark(remove_zeros_and_pad(tensor, indices) , pad_tensors_speed_up(value, indices, org_dim))


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
#' @examples
#' \dontrun{
#' FINN.seed(123)
#' # Now both R and Torch are seeded with 123, ensuring reproducible results
#' }
#'
#' @import torch
#' @export
FINN.seed <- function(seed) {
  set.seed(seed)
  torch::torch_manual_seed(seed)
}


checkPars = function(inputPar, parRange){
  if(is.list(inputPar)) inputPar = inputPar[[1]]
  if(is.vector(inputPar)) {
    Npar = 1
    Nsp = length(inputPar)
    inputPar = matrix(inputPar, nrow = Nsp, ncol = Npar)
    parRange = matrix(parRange, nrow = Npar, ncol = 2)
  }else if(is.matrix(inputPar)){
    Npar = ncol(inputPar)
    Nsp = nrow(inputPar)
  }else{
    stop("speciesPars and speciesPars_ranges must contain vectors or a matrices")
  }
  if(Npar != nrow(parRange)){
    stop("speciesPars and speciesPars_ranges must have the same number of parameters")
  }
  checkedPar = matrix(nrow = Nsp, ncol = Npar)
  j=1
  for(j in 1:Nsp){
    for(i in 1:Npar){
      lower = parRange[i,1]
      upper = parRange[i,2]
      checkedPar[j,i] = inputPar[j,i] < upper & inputPar[j,i] > lower
    }
  }
  return(list(invalid = any(!checkedPar), checkedPar = checkedPar, inputPar = inputPar, parRange = parRange))
}

checkParInput = function(speciesPars, speciesPars_ranges){
  checked = list()
  valid_pars = T
  for(i in names(speciesPars_ranges)){
    checked[[i]] = checkPars(speciesPars[[i]], speciesPars_ranges[[i]])
    if(checked[[i]]$invalid) valid_pars = F
  }
  if(!valid_pars){
    stop_message = paste("speciesPars must be within the range of speciesPars_ranges", sep = "\n")
    stop_message <- paste(stop_message, "The following parameters are out of range:", sep = "\n")
    for(i in names(speciesPars_ranges)){
      if(any(!checked[[i]]$checkedPar)){
        false_idx = which(!checked[[i]]$checkedPar, arr.ind = TRUE)
        for(ipar in 1:nrow(false_idx)){
          stop_message <- paste(stop_message, paste0("speciesPars$", i, "[", false_idx[ipar,1], ",", false_idx[ipar,2], "] = ", checked[[i]]$inputPar[false_idx[ipar,1], false_idx[ipar,2]]), sep = "\n")
        }
      }
    }
    stop_message <- paste(stop_message, "Adjust the upper and lower range in speciesPars_ranges or adjust the parameters in initSpecies", "", sep = "\n")
    stop(stop_message)
  }
}

# default ranges
default_speciesPars_ranges = list(
  parGrowth = rbind(
    c(0.01, 0.99),
    c(0.01, 4)
  ),
  parMort = rbind(
    c(0.01, 0.99),
    c(0, 4)
  ),
  parReg = c(0.01, 0.99),
  parHeight = c(0.3, 0.7)
)

to_r = function(par, numeric = FALSE) {
  if(numeric) {
    tmp = par |> as.numeric()
  } else {
    tmp = par |> as.matrix()
  }
  attributes(tmp) = append(attributes(tmp), list(requires_grad = par$requires_grad))
  return(tmp)
}
