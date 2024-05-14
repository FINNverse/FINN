library(R6)

CohortMat = R6::R6Class("CohortMat", public = list(
  dbh=NULL,
  nTree=NULL,
  Species=NULL,
  dims=c(50, 30, 10),
  sp = 10,
  device ='cpu',
  initialize = function(dbh=self$dbh, nTree=self$nTree, Species=self$Species, dims=self$dims, sp = self$sp, device = self$device) {
    self$dbh = if(is.null(dbh)) array(0.0, dim = dims) else dbh
    self$nTree = if(is.null(nTree)) array(0.0, dim = dims) else nTree
    self$Species = if(is.null(Species)) array(sample.int(sp, prod(dims),replace = TRUE), dim = dims) else Species
    self$dims = dims
    self$sp = sp
    self$device = torch::torch_device(device)

    if(!inherits(self$dbh, "torch_tensor")) self$dbh = torch::torch_tensor(self$dbh, dtype=torch_float32(), device=self$device)
    if(!inherits(self$nTree, "torch_tensor")) self$nTree = torch::torch_tensor(self$nTree, dtype=torch_float32(), device=self$device)
    if(!inherits(self$Species, "torch_tensor")) self$Species = torch::torch_tensor(self$Species, dtype=torch_int64(), device=self$device)

  }
))

