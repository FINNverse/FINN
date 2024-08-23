setPars = function(inputPar, parRange){
  # inputPar = torch::torch_tensor(inputPar, requires_grad = TRUE)
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
  out = matrix(nrow = Nsp, ncol = 0)
  # for(j in 1:Nsp){
  for(i in 1:Npar){
    lower = as.numeric(parRange[i,1,drop=FALSE])
    upper = as.numeric(parRange[i,2,drop=FALSE])
    out <- cbind(out, qlogis((inputPar[,i, drop=FALSE] - lower) / (upper - lower)))
  }
  if(Npar == 1) out = as.vector(out) # TODO replace vectors with matrices as species input everywhere
  # }
  # out <- torch::torch_tensor(out, requires_grad = TRUE, device = "cpu", dtype = "float32")
  out <- torch::torch_tensor(out, requires_grad = TRUE, device = self$device, dtype = self$dtype)
  return(out)
},

getPars = function(internalPar, parRange) {
  if (is.vector(parRange)) {
    # Case where internalPar is a 1D tensor and parRange is a vector
    Npar <- length(parRange) / 2
    lower <- parRange[1:Npar]
    upper <- parRange[(Npar + 1):(2 * Npar)]

    out <- internalPar$sigmoid() * (upper - lower) + lower
  } else {
    # Case where internalPar is a matrix and parRange is a matrix
    Npar <- ncol(internalPar)
    out <- list()
    for (i in 1:Npar) {
      lower <- parRange[i, 1, drop = FALSE]
      upper <- parRange[i, 2, drop = FALSE]
      out[[i]] <- internalPar[, i, drop = FALSE]$sigmoid() * (upper - lower) + lower
    }
    out <- torch::torch_cat(out, dim = 2L)
  }
  return(out)
}


self <- list(
  device = "cpu",
  dtype = "float32"
)


# default ranges
speciesPars_ranges = list(
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

sp = 2

#case 1 matrix input
  inputPar1 = cbind(
    runif(sp, min = speciesPars_ranges$parGrowth[1,1], speciesPars_ranges$parGrowth[1,2]),
    runif(sp, min = speciesPars_ranges$parGrowth[2,1], speciesPars_ranges$parGrowth[2,2])
  )

inputPar = inputPar1
parRange = speciesPars_ranges$parGrowth
internalPar = setPars(inputPar1, speciesPars_ranges$parGrowth)

inputPar
getPars(internalPar, speciesPars_ranges$parGrowth)
internalPar = setPars(inputPar1, speciesPars_ranges$parGrowth)


#case 2 vector input
inputPar2 =parHeight = runif(sp, min = speciesPars_ranges$parHeight[1], max = speciesPars_ranges$parHeight[2])
inputPar = inputPar2
parRange = speciesPars_ranges$parHeight

internalPar = setPars(inputPar2, speciesPars_ranges$parHeight)
inputPar
getPars(internalPar, speciesPars_ranges$parHeight)







