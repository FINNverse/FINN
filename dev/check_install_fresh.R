# Sourced by dev scripts before they touch system.file(). Aborts if the INSTALLED
# package is out of step with this source tree.
#
# Why: dev scripts and vignettes reach their data through system.file(), which
# resolves to the INSTALLED package - so editing R/ or re-running make_extdata.R
# without reinstalling means the code silently runs against stale data or stale
# functions. This has produced confident, wrong results three times; fail loudly.
check_install_fresh <- function(pkg = "FINN", src = ".") {
  inst_dir <- system.file("extdata", package = pkg)
  if (!nzchar(inst_dir)) stop(pkg, " is not installed - run: R CMD INSTALL . --no-docs")
  src_dir <- file.path(src, "inst", "extdata")
  if (dir.exists(src_dir)) {
    f <- list.files(src_dir)
    stale <- f[!vapply(f, function(x)
      file.exists(file.path(inst_dir, x)) &&
      identical(tools::md5sum(file.path(src_dir, x))[[1]], tools::md5sum(file.path(inst_dir, x))[[1]]),
      logical(1))]
    if (length(stale)) stop("INSTALLED inst/extdata is stale vs source:\n  ",
                            paste(stale, collapse = "\n  "),
                            "\nRun: R CMD INSTALL . --no-docs")
  }
  r_src <- list.files(file.path(src, "R"), pattern = "[.]R$", full.names = TRUE)
  desc  <- file.path(inst_dir, "..", "DESCRIPTION")
  if (length(r_src) && file.exists(desc)) {
    newer <- r_src[file.mtime(r_src) > file.mtime(desc)]
    if (length(newer)) stop("R/ is newer than the INSTALLED package:\n  ",
                            paste(basename(newer), collapse = ", "),
                            "\nRun: R CMD INSTALL . --no-docs")
  }
  invisible(TRUE)
}
