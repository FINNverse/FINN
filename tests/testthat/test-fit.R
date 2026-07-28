library(testthat)
library(FINN)
library(data.table)

# End-to-end smoke test: finn() -> fit() -> predict() on a tiny subset of the
# bundled FIA data (the real input format). Torch-gated and deliberately small
# (4 sites, 2 epochs, 2 patches) so it runs in ~1s. Setup lives in the shared
# fit_toy_finn() helper (also used by the xAI tests).

test_that("finn -> fit -> predict runs end to end", {
  skip_if_no_torch()

  f <- fit_toy_finn()          # builds, then fits; must not error
  expect_length(f$model$history, 2)   # one record per epoch

  # the fitted model must predict, returning the long site table
  sim <- predict(f$model, env = f$env, init_cohort = f$init, patches = 2, device = "cpu")
  expect_true(all(c("siteID", "year", "species", "variable", "value")
                  %in% names(sim$long$site)))
  expect_true(any(is.finite(sim$long$site$value)))
})
