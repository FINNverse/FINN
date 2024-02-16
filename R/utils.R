force_r = function(x) {
  if(inherits(x, "python.builtin.object")) return(reticulate::py_to_r( x ))
  else return(x)
}

softplus = function(x) log(1+exp(x))

check_installation = function() {
  # check if dependencies are installed
  torch_ = pyro_ = torch_optimizer_ = madgrad_ = c(crayon::red(cli::symbol$cross), 0)
  if(reticulate::py_module_available("torch")) torch_ =  c(crayon::green(cli::symbol$tick), 1)
  if(reticulate::py_module_available("torch_optimizer")) torch_optimizer_ =  c(crayon::green(cli::symbol$tick), 1)
  return(rbind("torch" = torch_,  "torch_optimizer" = torch_optimizer_))
}

force_r = function(x) {
  if(inherits(x, "python.builtin.object")) return(reticulate::py_to_r( x ))
  else return(x)
}



checkModel = function(object) {
  check_module()

  if(!reticulate::py_is_null_xptr(object$model)) return(object)

  object$model = object$get_model()

  if(inherits(object, c("sjSDM", "sjSDM_DNN"))){
    object$model$set_env_weights(lapply(object$weights, function(w) reticulate::r_to_py(w)$copy()))
    if(!is.null(object$spatial)) object$model$set_spatial_weights(lapply(object$spatial_weights, function(w) reticulate::r_to_py(w)$copy()))
    object$model$set_sigma(reticulate::r_to_py(object$sigma)$copy())
  }

  if(object$family$family$family == "nbinom") {
    object$model$set_theta(reticulate::r_to_py(object$theta)$copy())
  }
  return(object)
}


#' check module
#'
#' check if module is loaded
check_module = function(){
  if(is.null(pkg.env$fa)){
    .onLoad()
  }

  if(is.null(pkg.env$fa)) {
    stop("PyTorch not installed", call. = FALSE)
  }

  if(reticulate::py_is_null_xptr(pkg.env$FINN)) .onLoad()
}
