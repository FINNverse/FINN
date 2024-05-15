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
