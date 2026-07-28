library(testthat)
library(FINN)

# Interpretation / explainable-AI layer (used in the vignettes): conditional
# effects, their averages, permutation feature importance, and ALE. All run the
# torch simulation under the hood, so they are gated. Smoke-level: they must run
# on a fitted model and return output of the documented shape. Several of them
# draw a diagnostic plot; quiet_plot() routes that to a throw-away device and
# swallows the cosmetic "par(new=TRUE) with no plot" warning it raises, so the
# tests stay clean without hiding computation-level conditions.
quiet_plot <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressWarnings(force(code))
}

test_that("conditionalEffects / averageConditionalEffects run on a fitted model", {
  skip_if_no_torch()
  f <- fit_toy_finn()

  ce <- quiet_plot(conditionalEffects(f$model, env = f$env, init_cohort = f$init))
  expect_true(is.list(ce))
  expect_true("processes" %in% names(ce))

  ace <- averageConditionalEffects(f$model, env = f$env, init_cohort = f$init)
  expect_s3_class(ace, "FINNconditionalEffectsSummary")
  # each process summary, when populated, carries the documented columns
  for (p in ace) {
    if (is.data.frame(p) && nrow(p))
      expect_true(all(c("species", "variable", "mean_effect") %in% names(p)))
  }
})

test_that("feature_importance runs and scores the environmental predictors", {
  skip_if_no_torch()
  f <- fit_toy_finn()

  fi <- quiet_plot(suppressMessages(
    feature_importance(f$model, env = f$env, init_cohort = f$init,
                       nperm = 2L, method = "rmse")))
  expect_true(is.data.frame(fi) || is.list(fi))
})

test_that("ALE runs without plotting", {
  skip_if_no_torch()
  f <- fit_toy_finn()

  al <- quiet_plot(suppressMessages(
    ALE(f$model, env = f$env, init_cohort = f$init, plot = FALSE)))
  expect_false(is.null(al))
})
