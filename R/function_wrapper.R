#' Competition function
#'
#' @param dbh dbh of species, should be of dimensions [sites, patches, species, 1]
#' @param Species species identifier, should be of dimensions [sites, patches, species]
#' @param nTree number of trees, should be of dimensions [sites, patches, species]
#' @param parGlobal species global parameter, should be of dimensions [species] (so it should be a vector)
#' @param h height of species, if not null, should be [1,1]
#' @param minLight minimal light
#' @param dtype dtype
#' @param devie device
#'
#' @details
#' if [length(dim) == 1] this dimension is assumed to be the species dimension (works only for R objects!)
#'
#' @examples
#' competition(rep(50, 5),rep(1, 5), 0:4, rep(1., 5))
#' @export

competition = function(dbh, Species,nTree, parGlobal, h = NULL, minLight = 50.0, dtype = "float32", device= "cpu") {

  torch = pkg.env$torch

  if(dtype == "float32") dtype = torch$float32  else dtype = torch$float64
  dbh = FINN:::to_Tensor(dbh, dtype, device, TRUE, TRUE)
  Species = FINN:::to_Tensor(Species, dtype=torch$int64, device=device, TRUE, FALSE)
  nTree = FINN:::to_Tensor(nTree, dtype, device, TRUE, TRUE)

  parGlobal = FINN:::to_Tensor(parGlobal, dtype = dtype, device = device, FALSE, correct_zero = TRUE)
  if(!is.null(h)){
    if(h == 0) h = torch$zeros(list(1L, 1L)) else h = to_Tensor(h)
    }

  # AL = pkg.env$FINN$compF_P(dbh, Species, nTree, parGlobal, h= torch$zeros(list(1L, 1L)))$cpu()$data$numpy()
  result = pkg.env$FINN$compF_P(dbh, Species, nTree, parGlobal, h = h, minLight = minLight)$cpu()$data$numpy()
  return(force_r(result))

}


#' Growth function
#'
#' @param dbh dbh of species, should be of dimensions [sites, patches, species]
#' @param Species species identifier, should be of dimensions [sites, patches, species]
#' @param nTree number of trees, should be of dimensions [sites, patches, species]
#' @param parGlobal species global parameter, should be of dimensions [species] (so it should be a vector)
#' @param parGrowth species growth parameter, should be of dimensions [species, 2]
#' @param pred predictions for environmental niche, should be of dimensions [sites, species]
#' @param dtype dtype
#' @param devie device
#'
#' @details
#' if [length(dim) == 1] this dimension is assumed to be the species dimension (works only for R objects!)
#'
#' @examples
#' growth(rep(50, 5), 0:4,rep(1., 5), rep(1., 5), matrix(1., 5, 2), rep(0.5, 5))
#' @export
growth = function(dbh, Species, nTree, parGlobal, parGrowth, pred,  dtype = "float32", device= "cpu") {
  torch = pkg.env$torch

  if(dtype == "float32") dtype = torch$float32
  else dtype = torch$float64
  dbh = to_Tensor(dbh, dtype, device, TRUE, TRUE)
  pred = to_Tensor(pred, dtype=torch$int64, device=device, TRUE, FALSE, env = TRUE)
  Species = to_Tensor(Species, dtype=torch$int64, device=device, TRUE, FALSE)
  nTree = to_Tensor(nTree, dtype, device, TRUE, TRUE)

  parGrowth = to_Tensor(parGrowth, dtype = dtype, device = device, FALSE)
  parGlobal = to_Tensor(parGlobal, dtype = dtype, device = device, FALSE, correct_zero = TRUE)

  AL = pkg.env$FINN$compF_P(dbh, Species, nTree, parGlobal)
  result = pkg.env$FINN$growthFP(dbh, Species, parGlobal, parGrowth, pred, AL)$cpu()$data$numpy()
  return(force_r(result))
}


#' Mortality function
#'
#' @param dbh dbh of species, should be of dimensions [sites, patches, species]
#' @param Species species identifier, should be of dimensions [sites, patches, species]
#' @param nTree number of trees, should be of dimensions [sites, patches, species]
#' @param parGlobal species global parameter, should be of dimensions [species] (so it should be a vector)
#' @param parMort species mort parameter, should be of dimensions [species, 2]
#' @param pred predictions for environmental niche, should be of dimensions [sites, species]
#' @param dtype dtype
#' @param devie device
#'
#' @details
#' if [length(dim) == 1] this dimension is assumed to be the species dimension (works only for R objects!)
#'
#' @examples
#' mort(rep(50, 5), 0:4, rep(5, 5), rep(1., 5), matrix(1., 5, 2), rep(0.5, 5))
#'
#' @export
mort = function(dbh, Species, nTree, parGlobal, parMort, pred,  dtype = "float32", device= "cpu") {
  torch = pkg.env$torch

  if(dtype == "float32") dtype = torch$float32
  else dtype = torch$float64
  dbh = to_Tensor(dbh, dtype, device, TRUE, TRUE)
  nTree = to_Tensor(nTree, dtype, device, TRUE, TRUE)
  pred = to_Tensor(pred, dtype=torch$int64, device=device, TRUE, FALSE, env = TRUE)
  Species = to_Tensor(Species, dtype=torch$int64, device=device, TRUE, FALSE)

  parGlobal = to_Tensor(parGlobal, dtype = dtype, device = device, FALSE, correct_zero = TRUE)
  parMort = to_Tensor(parMort, dtype = dtype, device = device, FALSE)

  AL = pkg.env$FINN$compF_P(dbh, Species, nTree, parGlobal)
  result = pkg.env$FINN$mortFP(dbh, Species,nTree, parGlobal, parMort, pred, AL)$cpu()$data$numpy()
  return(force_r(result))
}


#' Reg function
#'
#' @param dbh dbh of species, should be of dimensions [sites, patches, species]
#' @param Species species identifier, should be of dimensions [sites, patches, species]
#' @param nTree number of trees should be of dimensions [sites, patches, species]
#' @param parGlobal species global parameter, should be of dimensions [species] (so it should be a vector)
#' @param parReg species growth parameter, should be of dimensions [species]
#' @param pred predictions for environmental niche, should be of dimensions [sites, species]
#' @param dtype dtype
#' @param devie device
#'
#' @details
#' if [length(dim) == 1] this dimension is assumed to be the species dimension (works only for R objects!)
#'
#' @examples
#' reg(rep(50, 5), 0:4, rep(1, 5), rep(1., 5), rep(1., 5), rep(0.5, 5))
#' @export
reg = function(dbh, Species, nTree, parGlobal, parReg, pred,  dtype = "float32", device= "cpu") {
  torch = pkg.env$torch

  if(dtype == "float32") dtype = torch$float32
  else dtype = torch$float64
  dbh = to_Tensor(dbh, dtype, device, TRUE, TRUE)
  nTree = to_Tensor(nTree, dtype, device, TRUE, TRUE)

  pred = to_Tensor(pred, dtype=torch$int64, device=device, TRUE, FALSE, env = TRUE)
  Species = to_Tensor(Species, dtype=torch$int64, device=device, TRUE, FALSE)


  parGlobal = to_Tensor(parGlobal, dtype = dtype, device = device, FALSE, correct_zero = TRUE)
  parReg = to_Tensor(parReg, dtype = dtype, device = device, FALSE, correct_zero = TRUE)

  AL = pkg.env$FINN$compF_P(dbh, Species, nTree, parGlobal, h= torch$zeros(list(1L, 1L)))
  result = pkg.env$FINN$regFP(dbh, Species, parGlobal, parReg, pred, AL)$cpu()$data$numpy()
  return(force_r(result))
}




to_Tensor = function(m, dtype, device, correct_dim = TRUE, add_last = TRUE, env = FALSE, correct_zero = FALSE) {
  tmp = m
  if(!inherits(m, "torch.tensor")) {
    if(correct_dim) {
      if(is.null(dim(m))) {
        if(add_last) m = array(m, dim = c(1, 1, length(m)))
        else m = array(m, dim = c(1, 1, length(m)))

        if(env) m = array(m, dim = c(1, length(m)))
      }
    }
    m = pkg.env$torch$tensor(m, dtype = dtype, device = device)
  } else  {
    m = m$to(device)
  }

  if(correct_zero) {
    if(is.vector(tmp) && length(tmp) == 1) {
      m = m$view(list(1L))
    }
  }
  return(m)
}
