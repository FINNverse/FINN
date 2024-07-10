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

  if(is.null(sp_max)) sp_max = as.numeric(Results[[1]]$shape[2])
  old_shape = labels$shape
  labels = labels$flatten(start_dim=1, end_dim = 2)

  for( v in 1:length(samples)) {
    da = samples[[v]]$flatten(start_dim=1, end_dim=2)
    Results[[v]]$add_(torch::torch_zeros(da$shape[1], sp_max, dtype=da$dtype, device = da$device)$scatter_add(2, labels, da$expand(c(da$shape[1], -1)))$view(c(old_shape[1:2], sp_max))$sum(dim = 2))
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
