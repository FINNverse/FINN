#' Install FINN and its dependencies
#'
#' @param conda path to conda
#' @param version version = "cpu" for CPU version, or "gpu" for GPU version. (note MacOS users have to install 'cuda' binaries by themselves)
#' @param restart_session Restart R session after installing (note this will
#'   only occur within RStudio).
#' @param ... not supported
#'
#' @return
#'
#' No return value, called for side effects (installation of 'python' dependencies).
#'
#' @export
install_FINN = function(conda = "auto",
                         version = c("cpu", "gpu"),
                         restart_session = TRUE, ...) {

  version = match.arg(version)

  method = "conda"
  envname = "r-FINN"

  # install conda if not installed
  conda = tryCatch(reticulate::conda_binary(), error = function(e) e)
  if(inherits(conda, "error")) {
    reticulate::install_miniconda(update = TRUE)
  }

  # get python dependencies
  pkgs = get_pkgs(version = version)

  # torch will be installed via pip on macOS because of mkl dependencies
  pip = FALSE
  channel = "pytorch"


  # install dependencies
  error = tryCatch({

    reticulate::py_install(
      pkgs$conda,
      envname = envname,
      method = "conda",
      conda = "auto",
      python_version = "3.10",
      channel = channel,
      pip = pip
    )

    reticulate::py_install(
      c("numpy", "pandas", "tqdm", "torch_optimizer"),
      envname = envname,
      method = "conda",
      conda = "auto",
      pip = TRUE
    )
  }, error = function(e) e)

  # check if instllation was successfull
  if(!inherits(error, "error")) {
    cli::cli_alert_success("\nInstallation complete.\n\n")

    if (restart_session && rstudioapi::hasFun("restartSession"))
      rstudioapi::restartSession()

    invisible(NULL)
  } else {
    cli::cli_alert_danger("\nInstallation failed... try to install manually PyTorch (install instructions: https://github.com/FINNverse/FINN\n")
    cli::cli_alert_info("If the installation still fails, please report the following error on https://github.com/FINNverse/FINN/issues\n")
    cli::cli_alert(error$message)
  }
}



is_windows = function() {
  identical(.Platform$OS.type, "windows")
}

is_unix = function() {
  identical(.Platform$OS.type, "unix")
}

is_osx = function() {
  Sys.info()["sysname"] == "Darwin"
}

is_linux = function() {
  identical(tolower(Sys.info()[["sysname"]]), "linux")
}

get_pkgs = function(version="cpu") {

  channel = "pytorch"
  if(is_windows() || is_linux() || is_unix()) {
    package = list()
    package$conda =
      switch(version,
             cpu = "pytorch torchvision torchaudio cpuonly",
             gpu = "pytorch torchvision torchaudio cudatoolkit=11.8")
  }

  if(is_osx()) {
    package = list()
    package$conda = "pytorch::pytorch torchvision torchaudio"
    if(version == "gpu") message("PyTorch does not provide cuda binaries for macOS, installing CPU version...\n")
  }

  packages = strsplit(unlist(package), " ", fixed = TRUE)
  return(packages)
}
