## The refactor's guard rail: the default objective must not move.
##
## Five configurations are fitted for three epochs with a non-zero learning rate
## and every parameter compared against a stored reference. The references describe
## the objective as of `mortality(growth_scale = 0.01)` and
## `createProcess(light_mask = TRUE)`; regenerate them only when
## a change to the objective is INTENDED, with
## `Rscript dev/make_identity_snapshots.R`, and say so in the commit message.
skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)

source(test_path("fixture-identity.R"))

snap_dir <- test_path("ref-identity")

test_that("the five reference configurations are bit-identical to the snapshots", {
  for (cfg in names(identity_configs())) {
    f <- file.path(snap_dir, paste0(cfg, ".rds"))
    expect_true(file.exists(f),
                info = paste("missing reference", f, "- see dev/make_identity_snapshots.R"))
    # 1e-6, not 1e-10: the same code drifts by ~1e-7 between BLAS implementations
    # (reported on macOS), while every objective change this guards against moves
    # the parameters by orders of magnitude more.
    expect_equal(identity_run(cfg), readRDS(f), tolerance = 1e-6, info = cfg)
  }
})
