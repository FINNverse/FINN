#' Sample poisson relaxed
#'
#' @param lmbd lambda
#' @param num_samples number of samples
#' @param temperature temperature
#'
#' @export
sample_poisson_relaxed = function(lmbd, num_samples=50, temperature = 1e-2) {
  t = (torch::torch_rand(c(num_samples,lmbd$shape))$log()$negative()/lmbd)$cumsum(dim=1L)
  relaxed_indicator = torch_sigmoid((1.0 - t) / temperature)
  N = relaxed_indicator$sum(1)
  return(N)
}


#' Basal aera multiplied with number of species
#'
#' @param dbh dbh
#' @param nTree number of Trees
#'
#' @export
BA_T_P = function(dbh, nTree) {
  return(pi*(dbh/100./2.)$pow(2.0)*nTree)
}

#' Index species
#'
#' @param pred predictions
#' @param Species Species index vector, must be int64
#'
#' @export
index_species = function(pred, Species) {
  shapes = pred$shape
  X_expanded = pred$unsqueeze(2)$expand(c(shapes[1], Species$shape[2], shapes[2]))
  pred = torch_gather(X_expanded, 3, Species)
  return(pred)
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

#' Mortality
#'
#' @param dbh dbh
#' @param Species Species
#' @param nTree nTree
#' @param parMort parMort
#' @param pred predictions
#' @param AL available light
#'
#' @export
mortFP = function(dbh, Species, nTree, parMort, pred, AL) {
  shade = 1-torch_sigmoid((AL + (1-parMort[,1][Species]) - 1)/1e-2)
  environment = index_species(pred, Species)
  gPSize = 0.1*(torch_clamp(dbh/(parMort[,2][Species]*100), min = 0.00001) )$pow(2.3) #.reshape([-1,1])
  predM = torch_clamp((shade*0.1+environment)+gPSize*1, min = 1-0.9999, max = 0.9999)
  #mort = torch.distributions.Beta(predM*nTree+0.00001, nTree - predM*nTree+0.00001).rsample()*nTree
  mort = binomial_from_gamma(nTree, predM)
  return( mort + mort$round()$detach() - mort$detach() )
}


#' Aggregate function
#'
#' @param labels labels
#' @param samples samples
#' @param Results Results
#'
#' @export
aggregate_results = function(labels, samples, Results, drop_rows = TRUE, sp_max = NULL) {

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


A = torch_ones(10, 2)
A$split(dim = 1L, split_size = 10L)


pad_tensors_speed_up = function(value, indices, org_dim) {
  # KK = torch::torch_split(value$flatten(start_dim = 1, end_dim = 2), split_size = org_dim[1]*org_dim[2], dim = 1) TODO does not work...bug?!
  KK = value$flatten(start_dim = 1, end_dim = 2)$split(1, 1)
  KK_new = vector("list", length = length(KK))
  for( i in 1:indices$shape[1]) {
    KK_new[[i]] = KK[[i]]$masked_select(indices[i, ]$detach())
  }
  return(torch::nn_utils_rnn_pad_sequence(KK_new, batch_first=TRUE))
}
