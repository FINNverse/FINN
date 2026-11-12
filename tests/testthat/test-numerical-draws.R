## The two numerical repairs, default since 0.2.0:
## createProcess(mortality_draw = "poisson") and
## createProcess(regeneration_sampler = "refined"). The defaults are guarded by
## test-identity.R (config `sample_reg` draws recruits); here the repaired paths
## themselves are exercised, and that "gamma"/"gaussian" still restore the
## published behaviour.
skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(torch)
library(data.table)

N <- 2e5

test_that("the defaults are the repaired draws and both arguments reach the model", {
  expect_equal(createProcess(func = FINN::mortality)$mortality_draw, "poisson")
  expect_equal(createProcess(func = FINN::regeneration)$regeneration_sampler, "refined")
  expect_equal(createHybrid()$mortality_draw, "poisson")
  d <- finn(N_species = 2)
  expect_equal(d$mortality_draw, "poisson")
  expect_equal(d$regeneration_sampler, "refined")
  # the published behaviour is still reachable, which is what makes an old fit reproducible
  m <- finn(N_species = 2,
            mortality_process    = createProcess(~1, func = FINN::mortality, mortality_draw = "gamma"),
            regeneration_process = createProcess(~1, func = FINN::regeneration, regeneration_sampler = "gaussian"))
  expect_equal(m$mortality_draw, "gamma")
  expect_equal(m$regeneration_sampler, "gaussian")
  expect_error(createProcess(func = FINN::mortality, mortality_draw = "binomial"))
})

test_that("a one-tree cohort is a Bernoulli trial, not a capped Poisson", {
  FINN.seed(1)
  p0 <- 0.26                                   # the p of a suppressed stem; the stem that decides thinning
  n <- torch_ones(N); p <- torch_full(N, p0)
  d <- binomial_from_poisson(n, p)
  expect_setequal(unique(as.numeric(d)), c(0, 1))
  expect_equal(as.numeric(d$mean()), p0, tolerance = 0.02)
  ## capping a Poisson at n = 1 instead would deliver (1-exp(-p))/p = 88% of p
  capped <- torch_minimum(torch_poisson(n * p), n)
  expect_equal(as.numeric(capped$mean()) / p0, (1 - exp(-p0)) / p0, tolerance = 0.02)
  expect_lt(as.numeric(capped$mean()) / p0, 0.9)
  ## and the default gamma draw, rounded as the forward pass rounds it, delivers ~0.7
  expect_lt(as.numeric(binomial_from_gamma(n, torch_full(N, 0.05))$round()$mean()) / 0.05, 0.8)
})

test_that("larger cohorts draw an integer Poisson capped at the cohort size", {
  FINN.seed(2)
  n <- torch_full(N, 5); p <- torch_full(N, 0.3)
  d <- binomial_from_poisson(n, p)
  expect_equal(as.numeric(d$max()), 5)
  expect_true(all(as.numeric(d) %% 1 == 0))
  expect_equal(as.numeric(d$mean()), 5 * 0.3, tolerance = 0.05)
})

test_that("the refined sampler has no recruitment floor", {
  FINN.seed(3)
  zero <- torch_zeros(5e5); theta <- torch_tensor(1)
  expect_lt(as.numeric(rnbinom_torch_refined(zero, theta)$mean()), 1e-5)
  ## the default's +0.001 offset and 1e-4 probability clamp are a floor of 0.0011
  ## per cell-year: the Gaussian draw rounds it away, an exact sampler delivers it
  ## (0.0011 / 0.0168 ha = 0.066 recruits ha^-1 yr^-1 per species)
  eps <- 1e-4
  probs <- torch_clamp(1 - theta / (theta + zero) + eps, 0, 1 - eps)
  logits <- torch_log(probs) - torch_log1p(-probs)
  lambda_floor <- (distr_gamma(theta, (-logits)$exp())$rsample() + 0.001)$abs()
  expect_equal(as.numeric(lambda_floor$mean()), 0.0011, tolerance = 0.1)
  expect_gt(as.numeric(FINN:::sample_poisson_refined(lambda_floor)$mean()), 5e-4)
})

test_that("both repaired paths run end to end and move the simulated stand", {
  Nsp <- 2
  mk <- function(draw, sampler) finn(N_species = Nsp, recruits_dbh = 5,
    competition_process  = createProcess(~0, func = FINN::competition),
    growth_process       = createProcess(~1, initEnv = matrix(c(0.5, 0.5), Nsp, 1),
                                         initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp, 2), func = FINN::growth),
    mortality_process    = createProcess(~1, initEnv = matrix(c(-2.5, -2.5), Nsp, 1),
                                         initSpecies = matrix(c(0.2, -0.2, 0.1, 0.1, 0, 0), Nsp, 3),
                                         func = FINN::mortality, mortality_draw = draw),
    regeneration_process = createProcess(~1, initEnv = matrix(c(1.5, 1.0), Nsp, 1),
                                         initSpecies = c(0.1, 0.1), func = FINN::regeneration,
                                         regeneration_sampler = sampler))
  it <- data.table(siteID = rep(1:2, each = 6), patchID = rep(rep(1:2, each = 3), 2),
                   species = rep(1:2, 6), dbh = rep(c(10, 25, 40), 4),
                   treeName = paste0("t", 1:12), trees = 4)
  ic <- makeInitCohorts(it, Nspecies = Nsp)
  env <- data.table(expand.grid(siteID = 1:2, year = 1:30))[, env1 := 0]
  sim <- function(draw, sampler) {
    FINN.seed(7)
    s <- simulateForest(mk(draw, sampler), env = env, init_cohort = ic, patches = 2L,
                        patch_size = 0.1, device = "cpu")
    s$long$site[variable == "trees", mean(value, na.rm = TRUE)]
  }
  base <- sim("gamma", "gaussian"); rep_ <- sim("poisson", "refined")
  expect_true(is.finite(base) && is.finite(rep_))
  expect_false(isTRUE(all.equal(base, rep_)))
})

## binomial_from_bernoulli() returned exactly zero for every input until 0.2.0
## (FINNverse/FINN#62): the trial mask was inverted, and the reparameterised draw
## was fed to a second torch_bernoulli() instead of being thresholded. Both moments
## are asserted so neither can come back silently.
test_that("binomial_from_bernoulli matches the binomial in both moments", {
  skip_if_not(torch::torch_is_installed())
  FINN.seed(1)
  for (nn in c(1, 10, 40)) for (pp in c(0.1, 0.3)) {
    n <- torch::torch_full(c(400L, 100L, 1L), nn)
    p <- torch::torch_full(c(400L, 100L, 1L), pp)
    b <- binomial_from_bernoulli(n, p)
    expect_equal(as.numeric(b$mean()), nn * pp, tolerance = 0.02)
    expect_equal(as.numeric(b$var()), nn * pp * (1 - pp), tolerance = 0.05)
  }
})
