# vignettes/build.R — precompile the FINN vignettes (cito-style).
#
# WHY: the vignettes train torch models. CRAN machines have no libtorch backend
# and no time budget, so the vignettes must NOT run at build time. We therefore
# keep the real code in `*.Rmd.orig`, run it ONCE here, and commit the knitted
# `*.Rmd` with all results and figures baked in. R CMD build then only has to
# run pandoc over static markdown.
#
# CONSEQUENCES (this is the whole point):
#   * Nothing heavy ships. The .pt models and .rds caches are gone; only the
#     knitted .Rmd + its figure PNGs are in the package. Vignette size is no
#     longer a constraint on what we can demonstrate.
#   * No drift. The code shown in a vignette is exactly the code that produced
#     the output next to it, because one source generated both. (The previous
#     `eval:false` + cached-.rds scheme could silently disagree.)
#
# USAGE — from the package root, with a torch backend available:
#
#     Rscript vignettes/build.R                 # knit all vignettes, then build the site
#     Rscript vignettes/build.R fit-to-fia      # knit one, then build the site
#     Rscript vignettes/build.R --no-site       # knit only, skip pkgdown
#
# This is SLOW on purpose (fit-to-fia trains two models, ~10 min). Run it after
# changing a vignette's analysis, then commit the regenerated .Rmd + figure dir.
# See data-raw/README.md for the full data pipeline.
#
# The pkgdown step renders the website (https://finnverse.github.io/FINN) into
# docs/. It reuses the knitted .Rmd above, so the site never retrains anything.

stopifnot(requireNamespace("knitr", quietly = TRUE))

# Each vignette knits in vignettes/ so that fig.path (and the figure folder that
# ships next to the .Rmd) resolves relative to the vignette itself.
#
# Naming follows cito: "<LETTER>-<Name>.Rmd", with figures in a matching folder
# via fig.path = "<LETTER>/<LETTER>-". The letter prefix is what fixes the
# READING ORDER — pkgdown auto-discovers articles and sorts them alphabetically
# by filename, so without it succession lands last instead of second.
vignettes <- c("A-Introduction_to_FINN", "B-Succession_demo",
               "C-Data_preparation", "D-Fit_to_FIA")

args     <- commandArgs(trailingOnly = TRUE)
do_site  <- !("--no-site" %in% args)
args     <- setdiff(args, "--no-site")
if (length(args)) {
  unknown <- setdiff(args, vignettes)
  if (length(unknown)) stop("unknown vignette(s): ", paste(unknown, collapse = ", "))
  vignettes <- args
}

owd <- setwd("vignettes")
on.exit(setwd(owd), add = TRUE)

for (v in vignettes) {
  src <- paste0(v, ".Rmd.orig")
  out <- paste0(v, ".Rmd")
  if (!file.exists(src)) { message("[", v, "] no ", src, " - skipping"); next }
  message("[", v, "] knitting ", src, " -> ", out, " ...")
  t0 <- Sys.time()
  knitr::knit(src, output = out, quiet = FALSE)
  message("[", v, "] done in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}

message("\nKnitted: ", paste(vignettes, collapse = ", "),
        "\nCommit the regenerated .Rmd files and their figure folders.")

setwd(owd)

# --- website -----------------------------------------------------------------
# Same final step as cito's build.R: render the pkgdown site from the *knitted*
# vignettes (plus the website-only articles under vignettes/articles/, which are
# .Rbuildignore'd and therefore never ship to CRAN).
if (do_site) {
  if (!requireNamespace("pkgdown", quietly = TRUE)) {
    message("\npkgdown not installed - skipping the site. install.packages('pkgdown')")
  } else {
    message("\nBuilding the pkgdown site into docs/ ...")
    pkgdown::build_site(preview = FALSE)
    message("Site built: docs/index.html")
  }
}
