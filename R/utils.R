force_r = function(x) {
  if(inherits(x, "python.builtin.object")) return(reticulate::py_to_r( x ))
  else return(x)
}

softplus = function(x) log(1+exp(x))

check_installation = function() {
  # check if dependencies are installed
  torch_ = pyro_ = torch_optimizer_ = madgrad_ = c(crayon::red(cli::symbol$cross), 0)
  if(reticulate::py_module_available("torch")) torch_ =  c(crayon::green(cli::symbol$tick), 1)
  if(reticulate::py_module_available("pyro")) pyro_ =  c(crayon::green(cli::symbol$tick), 1)
  if(reticulate::py_module_available("torch_optimizer")) torch_optimizer_ =  c(crayon::green(cli::symbol$tick), 1)
  if(reticulate::py_module_available("madgrad")) madgrad_ =  c(crayon::green(cli::symbol$tick), 1)
  return(rbind("torch" = torch_,  "torch_optimizer" = torch_optimizer_, "pyro" = pyro_, "madgrad" = madgrad_))
}

force_r = function(x) {
  if(inherits(x, "python.builtin.object")) return(reticulate::py_to_r( x ))
  else return(x)
}
