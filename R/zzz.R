pkg.env = new.env()
pkg.env$name = "r-FINN"
pkg.env$torch = NULL
pkg.env$FINN = NULL



.onLoad = function(libname, pkgname){
  msg( text_col( cli::rule(left = "Attaching FINN", right = utils::packageVersion("FINN")) ), startup = TRUE)

  # causes problems on macOS systems
  if( is_osx() ) Sys.setenv( KMP_DUPLICATE_LIB_OK=TRUE )

  # load r-sjsdm environment
  success_env = try({
    envs = reticulate::conda_list()
    env_path = envs[which(envs$name %in% "r-FINN", arr.ind = TRUE), 2]
    reticulate::use_python(env_path, required = TRUE)
  }, silent = TRUE)

  # check if dependencies are installed
  check = check_installation()

  # load modules only if dependencies are available
  modules_available = any(check[,2] == "0")
  if(!modules_available) {
    # load torch
    pkg.env$torch = reticulate::import("torch", delay_load = FALSE, convert = FALSE )

    # 'compile' and load sjSDM python package
    path = system.file("python", package = "FINN")
    try({
      compile = reticulate::import("compileall", delay_load = FALSE)
      tmp = compile$compile_dir(paste0(path, "/FINN_python"),quiet = 2L,force=TRUE)
    }, silent = FALSE)
    pkg.env$FINN = reticulate::import_from_path("FINN_python", path, delay_load = FALSE, convert = FALSE)
    check= cbind(check, crayon::black( c(pkg.env$torch$`__version__`, rep("", 1))))
  }

  check[,2] = crayon::black( rownames(check) )
  check = cbind(check, "\n")

  msg(paste0(apply(check, 1, function(d) paste0(d,collapse = " "))), startup = TRUE)

  if(modules_available) {
    msg( crayon::red( "Torch or other dependencies not found:" ), startup = TRUE)
    info =
      c(
        "\t1. Use install_FINN() to install Pytorch and conda automatically \n",
        "\t2. Installation trouble shooting guide: ?installation_help \n")
    msg( info, startup = TRUE )
  }
  invisible()
}

# copied from the tidyverse package
msg <- function(..., startup = FALSE) {
  if (startup) {
    if (!isTRUE(getOption("tidyverse.quiet"))) {
      packageStartupMessage(text_col(...))
    }
  } else {
    message(text_col(...))
  }
}

# copied from the tidyverse package
text_col <- function(x) {
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  theme <- rstudioapi::getThemeInfo()

  if (isTRUE(theme$dark)) crayon::white(x) else crayon::black(x)

}






# FINN = function(env,
#                 reg = ~ Temp + Precp,
#                 growth = ~ Temp,
#                 mort = ~ Temp + Precp,
#                 models = list(reg = NULL,
#                               growth = dnn(c(20L, 20L), lambda = 0.001),
#                               # Neural network with two hidden layers, each with 20 units
#                               mort = c(20L, 20L)
#                 ) # lineares Modell
#
#                 parameters = list(...),
#                 simulate = TRUE
#
# )
