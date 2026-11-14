library(testthat)
library(FINN)
library(data.table)

# FINN::regeneration_adult_effect: the density of conspecific adults as an
# effect on recruitment, f(A) = f0 + (1 - f0) A / (A + A_half), instead of the
# gate of regeneration_adult (f0 = 0). The unit tests call the process function
# on hand-built stands through the model it is bound to; the last test fits a
# tiny FIA subset. Torch-gated.

# len = 1 -> shared parameters, len = Nsp -> one per species
build_effect_model <- function(Nsp, logitf0 = qlogis(0.3), logx = log(20), logA = log(20),
                               len = Nsp, logitf0_len = len) {
  FINN.seed(1)
  finn(
    N_species            = Nsp, recruits_dbh = 12.9,
    competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    growth_process       = createProcess(~ temp + prec, FINN::growth,    optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = createProcess(~ temp + prec, FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
    regeneration_process = createProcess(
      ~ temp + prec, FINN::regeneration_adult_effect,
      custom_parameters = list(reg_logK          = rep(log(50), len),
                               reg_adult_logx    = rep(logx, len),
                               reg_adult_logA    = rep(logA, len),
                               reg_adult_logitf0 = rep_len(logitf0, logitf0_len)),
      optimizeSpecies = TRUE, optimizeEnv = TRUE)
  )
}

# one site, 2 patches, 2 species; cohort c has species sp[c], dbh d[c], n[c]
# stems. `func` (default: the model's own regeneration function) is bound to the
# model, so the sibling processes can be evaluated with the same parameters.
effect_call <- function(m, d, n, sp, debug = FALSE, stand = TRUE, func = NULL) {
  m$patch_size_ha <- 0.1
  shp <- c(1L, 2L, length(d))
  dbh   <- torch_tensor(array(rep(d,  each = 2), dim = shp))
  trees <- torch_tensor(array(rep(n,  each = 2), dim = shp))
  species <- torch_tensor(array(rep(sp, each = 2), dim = shp), dtype = torch_int64())
  pred  <- torch_zeros(1L, 2L)                        # exp(0) = 1 environment
  light <- torch_ones(1L, 2L, 2L)                     # full light
  if (is.null(func)) func <- m$regeneration_func
  else environment(func) <- m$.__enclos_env__
  args <- list(species = species, parReg = torch_tensor(c(0.1, 0.1)), pred = pred, light = light)
  # dbh/trees only for a signature that has them, as forward() does
  if (stand && all(c("dbh", "trees") %in% names(formals(func))))
    args <- c(args, list(dbh = dbh, trees = trees))
  args$debug <- debug
  do.call(func, args)
}

# the potential recruitment at full light, parReg 0.1 and exp(env) = 1
regP_full <- (plogis(9) - plogis(-1)) / (1 - plogis(-9))
bh <- function(m, K = 50) K * m / (K + m)

test_that("shared and per-species init lengths both work", {
  skip_if_no_torch()
  shared <- build_effect_model(2L, len = 1L)
  expect_equal(shared$reg_adult_logitf0$numel(), 1L)
  expect_equal(shared$reg_adult_logx$numel(), 1L)
  out_s <- effect_call(shared, d = c(30, 14), n = c(3, 2), sp = c(1L, 2L))
  expect_equal(out_s$shape, c(1L, 2L, 2L))
  expect_true(all(is.finite(as.numeric(out_s))))

  per_sp <- build_effect_model(2L, logitf0 = c(qlogis(0.1), qlogis(0.6)), len = 2L)
  expect_equal(per_sp$reg_adult_logitf0$numel(), 2L)
  out_p <- effect_call(per_sp, d = c(30, 30), n = c(0, 0), sp = c(1L, 2L), debug = TRUE)
  expect_equal(out_p$mean$shape, c(1L, 2L, 2L))
  # each species gets its own f0 on the empty stand: species is the LAST dimension
  expect_equal(as.numeric(out_p$seed[1, 1, ]), c(0.1, 0.6), tolerance = 1e-6)

  # a shared f0 next to per-species x, A_half and cap
  mixed <- build_effect_model(2L, len = 2L, logitf0_len = 1L)
  expect_equal(mixed$reg_adult_logitf0$numel(), 1L)
  expect_true(all(is.finite(as.numeric(effect_call(mixed, d = c(30, 14), n = c(3, 2), sp = c(1L, 2L))))))
})

test_that("regeneration_adult_effect without its parameters errors clearly", {
  skip_if_no_torch()
  m <- finn(
    N_species = 2L, recruits_dbh = 12.9,
    competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    regeneration_process = createProcess(~ temp + prec, FINN::regeneration_adult_effect,
                                         custom_parameters = list(reg_adult_logx = rep(log(20), 2),
                                                                  reg_adult_logA = rep(log(20), 2)),
                                         optimizeSpecies = TRUE, optimizeEnv = TRUE))
  expect_false(is.null(m$reg_logK))     # the cap is created by finn(), as for regeneration_adult
  expect_error(effect_call(m, d = c(30, 30), n = c(1, 1), sp = c(1L, 2L)), "reg_adult_logitf0")
})

test_that("f0 -> 1 is regeneration_capped and f0 -> 0 is regeneration_adult", {
  skip_if_no_torch()
  d <- c(30, 14); n <- c(3, 2); sp <- c(1L, 2L)

  one <- build_effect_model(2L, logitf0 = 30)
  capped <- effect_call(one, d, n, sp, func = FINN::regeneration_capped)
  expect_equal(as.numeric(effect_call(one, d, n, sp)), as.numeric(capped), tolerance = 1e-6)

  zero <- build_effect_model(2L, logitf0 = -30)
  adult <- effect_call(zero, d, n, sp, func = FINN::regeneration_adult)
  expect_equal(as.numeric(effect_call(zero, d, n, sp)), as.numeric(adult), tolerance = 1e-6)

  # in between, adults still raise recruitment above the share without them
  mid <- effect_call(build_effect_model(2L), d = c(30, 30), n = c(10, 0), sp = sp, debug = TRUE)
  expect_gt(as.numeric(mid$seed[1, 1, 1]), 0.3)
  expect_equal(as.numeric(mid$seed[1, 1, 2]), 0.3, tolerance = 1e-6)  # species 2 has no adults
})

test_that("an empty stand recruits at potential x f0, far above the floor", {
  skip_if_no_torch()
  m <- build_effect_model(2L, logitf0 = qlogis(0.3))
  expected <- bh(regP_full * 0.3 + 1e-3)

  empty <- effect_call(m, d = c(30, 30), n = c(0, 0), sp = c(1L, 2L), debug = TRUE)
  expect_equal(as.numeric(empty$adults), c(0, 0))
  expect_equal(as.numeric(empty$mean), rep(expected, 4), tolerance = 1e-5)
  expect_true(all(as.numeric(empty$mean) > 100 * 1e-3))

  # a stand with no cohorts at all gives the same
  none <- effect_call(m, d = numeric(0), n = numeric(0), sp = integer(0))
  expect_equal(as.numeric(none), rep(expected, 4), tolerance = 1e-5)

  # regeneration_adult on the same empty stand is stuck at the floor
  gate <- as.numeric(effect_call(m, d = c(30, 30), n = c(0, 0), sp = c(1L, 2L), func = FINN::regeneration_adult))
  expect_equal(gate, rep(bh(1e-3), 4), tolerance = 1e-6)

  # and the gradient reaches f0 there: a higher f0 gives more recruits
  out <- effect_call(m, d = c(30, 30), n = c(0, 0), sp = c(1L, 2L))
  out$sum()$backward()
  expect_true(all(as.numeric(m$reg_adult_logitf0$grad) > 0))
})

test_that("without dbh/trees the seed source counts as present (f = 1)", {
  skip_if_no_torch()
  m <- build_effect_model(2L)
  out <- effect_call(m, d = c(30, 30), n = c(0, 0), sp = c(1L, 2L), debug = TRUE, stand = FALSE)
  expect_equal(as.numeric(out$seed), rep(1, 4))
})

test_that("regeneration_adult_effect fits, and f0 moves", {
  skip_if_no_torch()
  obs <- copy(FINN::fia_obs_dt); env <- copy(FINN::fia_env_dt); it <- copy(FINN::fia_init_trees)
  sites <- head(sort(unique(obs$siteID)), 4L)
  Nsp <- max(obs$species)
  m <- build_effect_model(Nsp)
  before <- as.numeric(m$reg_adult_logitf0)
  suppressMessages(
    fit(m, env = env[siteID %in% sites], data = obs[siteID %in% sites],
        init_cohort = makeInitCohorts(it[siteID %in% sites], Nspecies = Nsp), device = "cpu",
        epochs = 3L, patches = 4, patch_size = 0.06, lr = 0.02, plot_progress = FALSE)
  )
  after <- as.numeric(m$reg_adult_logitf0)
  expect_true(all(is.finite(after)))
  expect_false(isTRUE(all.equal(before, after)))
})
