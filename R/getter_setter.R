
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


set_weights_nnGrowthEnv = function(weights) {

  self$nnGrowthEnv = set_weights(weights, self$nnGrowthEnv)

}

set_weights_nnMortEnv = function(weights) {

  self$nnMortEnv = set_weights(weights, self$nnMortEnv)

}

set_weights_nnRegEnv = function(weights) {

  self$nnRegEnv  = set_weights(weights, self$nnRegEnv )

}

set_weights = function(weights, NN) {
  torch::with_no_grad({
    counter = 1
    for(i in 1:length(NN)) {
      if(inherits(NN[[i]], "nn_linear")) {
        NN$modules[[i]]$parameters$`0.weight`$set_data(weights[[counter]])
        counter <<- counter + 1
        if(!is.null(NN[[i]]$bias)) {
          NN$modules[[i]]$parameters$`0.bias`$set_data(weights[[counter]])
        }
      }
    }
  })
  return(NN)
}



get_parMort =  function() {
  #return(torch::torch_cat(list(self$parMort[,1,drop=FALSE]$sigmoid(), self$parMort[,2,drop=FALSE]$sigmoid()*4.0), dim = 2L ))
  return(torch::torch_cat(list(self$parMort[,1,drop=FALSE]$clamp(min = 0., max = 1.0), self$parMort[,2,drop=FALSE]$clamp(min = 0.0, max = 4.0)), dim = 2L ))
}
get_parGrowth =  function() {
  #return(torch::torch_cat(list(self$parGrowth[,1,drop=FALSE]$sigmoid(), self$parGrowth[,2,drop=FALSE]$exp()), dim = 2L ))
  return(torch::torch_cat(list(self$parGrowth[,1,drop=FALSE]$clamp(min = 0.0, max = 1.0), self$parGrowth[,2,drop=FALSE]$clamp(min = 0.0)), dim = 2L ))
}
get_parHeight =  function() {
  # return(self$parHeight$sigmoid())
  return(self$parHeight$clamp(min = 0.0, max = 1.0))
}
get_parReg = function() {
  #return(self$parReg$sigmoid())
  return(self$parReg$clamp(min = 0.0, max = 1.0))
}


set_parMort =  function(value) {
  #return(cbind(binomial()$linkfun(value[,1]), binomial()$linkfun(value[,2]/4.0)))
  #self$parMort = torch::torch_tensor(cbind(stats::binomial()$linkfun(value[,1]), stats::binomial()$linkfun(value[,2]/4.0)), requires_grad = TRUE, device = self$device, dtype=self$dtype)
  self$parMort = torch::torch_tensor(value, requires_grad = TRUE, device = self$device, dtype=self$dtype)
}
set_parGrowth =  function(value) {
  #self$parGrowth = torch::torch_tensor(cbind(stats::binomial()$linkfun(value[,1]), log(value[,2])), requires_grad = TRUE, device = self$device, dtype=self$dtype)
  self$parGrowth = torch::torch_tensor(value, requires_grad = TRUE, device = self$device, dtype=self$dtype)
}
set_parHeight =  function(value) {
  #self$parHeight = torch_tensor(stats::binomial()$linkfun(value), requires_grad = TRUE, device = self$device, dtype=self$dtype )
  self$parHeight = torch_tensor(value, requires_grad = TRUE, device = self$device, dtype=self$dtype )
}
set_parReg = function(value) {
  #self$parReg = torch_tensor(stats::binomial()$linkfun(value), requires_grad = TRUE, device = self$device, dtype=self$dtype )
  self$parReg = torch_tensor(value, requires_grad = TRUE, device = self$device, dtype=self$dtype )
}
