## Regenerate the reference objects of tests/testthat/test-identity.R.
## Run ONLY when a change to the default objective is intended.
pkgload::load_all(".", quiet = TRUE)
library(testthat)
source("tests/testthat/fixture-identity.R")
d <- "tests/testthat/ref-identity"
dir.create(d, recursive = TRUE, showWarnings = FALSE)
for (cfg in names(identity_configs())) {
  saveRDS(identity_run(cfg), file.path(d, paste0(cfg, ".rds")))
  message("wrote ", cfg)
}
