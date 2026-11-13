## Per-site path: each response's loss is ONE mean over all units scored in a
## batch (trees, n_at_risk, growth_n, cells), not a sum of per-observation-year
## means. Checked here, on a toy where two sites are observed once (year 10,
## a 10-year interval) and two are observed twice (year 4 and year 10):
##   (i)   one scored observation year: the loss is that year's mean, exactly
##         as before (the weights are n/n = 1);
##   (ii)  two observation years: the loss equals the per-unit pooled mean
##         computed by hand from the model's own predictions, for the per-tree
##         binomial and for the n-weighted cell terms, and differs from the old
##         sum / mean of the years' means;
##   (iii) the "auto" baselines are the same per-unit pooled functional, so the
##         auto-weighted history is the pooled raw loss over the baseline.
skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())
library(data.table)

Nsp <- 2; Nsites <- 4; Tmax <- 10; SEED <- 7
mk_model <- function() finn(N_species = Nsp, recruits_dbh = 5,
  competition_process  = createProcess(~0, func = FINN::competition),
  growth_process       = createProcess(~1 + env1, initEnv = matrix(c(0.5, 0.5, 0.4, -0.4), Nsp, 2),
                                       initSpecies = matrix(c(0.1, 0.2, 0.05, 0.05), Nsp, 2), func = FINN::growth),
  mortality_process    = createProcess(~1 + env1, initEnv = matrix(c(-2.5, -2.0, 0.6, -0.6), Nsp, 2),
                                       initSpecies = matrix(c(0.2, -0.2, 0.1, 0.1, 0, 0), Nsp, 3), func = FINN::mortality),
  regeneration_process = createProcess(~1 + env1, initEnv = matrix(c(1.5, 1.0, 0.5, -0.5), Nsp, 2),
                                       initSpecies = c(0.1, 0.1), func = FINN::regeneration, sample_regeneration = FALSE))
init_trees <- data.table(siteID = rep(1:Nsites, each = 8), patchID = rep(rep(1:2, each = 4), Nsites),
                         species = rep(1:2, 16), dbh = rep(c(10, 20, 30, 40), 8) + rep(0:3, each = 8) * 5,
                         treeName = paste0("t", 1:32), trees = 3)
ic  <- makeInitCohorts(init_trees, Nspecies = Nsp)
env <- data.table(expand.grid(siteID = 1:Nsites, year = 1:Tmax)); env[, env1 := (siteID - 2.5) * 0.5 + sin(year)]
OPT <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam

## observation schedule: sites 1-2 once (year 10, period 10), sites 3-4 twice (4 + 6)
sched <- function(split = TRUE) rbindlist(lapply(1:Nsites, function(s) {
  yrs <- if (split && s >= 3) c(4L, 10L) else 10L
  data.table(siteID = s, year = yrs, period_length = diff(c(0L, yrs)))
}))
mk_obs <- function(split = TRUE) {
  d <- sched(split)[, .(species = 1:Nsp), by = .(siteID, year, period_length)]
  d[, `:=`(dbh = 15 + siteID + species, ba = 8 + species, trees = 40 + 3 * siteID,
           growth = 0.01 * species + 0.002 * siteID, growth_n = siteID + 2 * species,
           mort = 0.03 * species, n_at_risk = 5 * siteID + species, reg = 2 + siteID,
           species_name = paste0("sp", species))]
  d[year == 4L, `:=`(growth_n = 1, n_at_risk = 1)]   # few units in the short interval
  d[order(siteID, year, species)]
}
## trees: many in the decadal year-10 census, few at year 4
mk_trees <- function(split = TRUE) {
  s <- sched(split)
  rbindlist(lapply(seq_len(nrow(s)), function(k) {
    n <- if (s$year[k] == 4L) 2L else if (s$siteID[k] <= 2) 12L else 6L
    data.table(siteID = s$siteID[k], year = s$year[k], patchID = 1L, species = rep_len(1:Nsp, n),
               dbh = seq(8, 40, length.out = n), trees = 1, growth = NA_real_,
               died = rep_len(c(0, 1, 0), n), period_length = s$period_length[k])
  }))
}
run_fit <- function(obs, tree_dt = NULL, family, weights = rep(1, 6), agg = NULL,
                    weight_scale = "raw") {
  FINN.seed(SEED); mm <- mk_model(); mm$period_mode <- "per_site"
  # weight_scale = "raw": this file checks the pooling arithmetic by hand, so the
  # reported loss must be the raw one and not the baseline-scaled objective
  suppressMessages(fit(mm, env = env, data = obs, init_cohort = ic, patches = 2L, patch_size = 0.1,
      env_autoscale = FALSE, loss_family = family, loss_aggregation = agg, tree_data = tree_dt,
      epochs = 1L, lr = 0, loss_weights = weights, weight_scale = weight_scale, batchsize = Nsites, shuffle = FALSE,
      control = finnControl(optimizer = OPT, clip_norm = 1e6), plot_progress = FALSE, device = "cpu"))   # unclipped grads
  mm
}
bce <- function(p, y) { p <- pmin(pmax(p, 1e-6), 1 - 1e-6); -(y * log(p) + (1 - y) * log(1 - p)) }
FAM_TREE <- c(dbh = "none", ba = "none", trees = "none", growth = "none", mortality = "binomial", regeneration = "none")

test_that("(i) one scored observation year: the loss is that year's mean, as before", {
  td <- mk_trees(split = FALSE)
  m  <- run_fit(mk_obs(split = FALSE), td, FAM_TREE, agg = c(mortality = "tree"))
  h  <- as.numeric(m$history[[1]])
  # the reported loss is the step's own loss vector, unweighted
  expect_identical(h, m$loss_raw)
  P  <- predictTrees(m, td, env, patch_size = 0.1)
  expect_equal(h[5], mean(bce(P$mort_pred, td$died)), tolerance = 1e-5)
})

test_that("(ii) per-tree binomial: the loss is the mean over all trees of the batch", {
  td <- mk_trees()
  m  <- run_fit(mk_obs(), td, FAM_TREE, agg = c(mortality = "tree"))
  h5 <- as.numeric(m$history[[1]])[5]
  P  <- predictTrees(m, td, env, patch_size = 0.1)
  l  <- bce(P$mort_pred, td$died)
  pooled <- mean(l)
  by_year <- tapply(l, td$year, mean)
  expect_equal(h5, pooled, tolerance = 1e-5)
  # the toy discriminates: 8 trees at year 4, 36 at year 10
  expect_gt(abs(sum(by_year) - pooled), 1e-2)        # old objective (sum of the years' means)
  expect_gt(abs(mean(by_year) - pooled), 1e-3)       # old report (mean of the years' means)
})

## the model's site-level predictions aggregated over each observation's window,
## as forward() does: growth and mortality the window mean, regeneration the sum
window_preds <- function(m, obs) {
  p <- dcast(m$pred$long$site, siteID + year + species ~ variable)
  for (cc in c("siteID", "year", "species")) p[, (cc) := as.integer(as.character(get(cc)))]
  rows <- unique(obs[, .(siteID, year, period_length)])
  rbindlist(lapply(seq_len(nrow(rows)), function(k) {
    s <- rows$siteID[k]; y <- rows$year[k]; pl <- rows$period_length[k]
    st <- p[siteID == s & year == y, .(species, dbh_p = dbh, ba_p = ba, trees_p = trees)]
    w  <- p[siteID == s & year > y - pl & year <= y,
            .(growth_p = mean(growth), mort_p = mean(mort), reg_p = sum(r_mean_ha)), by = species]
    merge(st, w, by = "species")[, `:=`(siteID = s, year = y)]
  }))
}
FAM_CELL <- c(dbh = "mse", ba = "mse", trees = "mse", growth = "mse", mortality = "binomial", regeneration = "mse")

test_that("(ii) n-weighted cell terms: one pooled mean over n_at_risk / growth_n / cells", {
  obs <- mk_obs()
  m   <- run_fit(obs, family = FAM_CELL)
  h   <- as.numeric(m$history[[1]])
  D   <- merge(obs, window_preds(m, obs), by = c("siteID", "year", "species"))
  hand <- c(D[, mean((dbh_p - dbh)^2)], D[, mean((ba_p - ba)^2)], D[, mean((trees_p - trees)^2)],
            D[, sum(growth_n * (growth_p - growth)^2) / sum(growth_n)],
            D[, sum(n_at_risk * bce(mort_p, mort)) / sum(n_at_risk)],
            D[, mean((reg_p - reg)^2)])
  expect_equal(h[1:6], hand, tolerance = 1e-5)
  # old rule for mortality: the years' n-weighted means, summed
  old_m <- D[, .(v = sum(n_at_risk * bce(mort_p, mort)) / sum(n_at_risk)), by = year][, sum(v)]
  expect_gt(abs(old_m - h[5]), 1e-3)
})

test_that("(iii) the auto baselines use the same per-unit pooling", {
  obs <- mk_obs(); td <- mk_trees()
  fam <- replace(FAM_CELL, c("dbh", "ba", "trees"), "none")
  raw  <- run_fit(obs, td, fam, agg = c(mortality = "tree"))
  auto <- run_fit(obs, td, fam, weights = "auto", agg = c(mortality = "tree"), weight_scale = "auto")
  b <- auto$loss_baseline
  # intercept-only loss over ALL units of all observation years
  y <- td$died; mu <- mean(y)
  expect_equal(b[["mortality"]], mean(bce(rep(mu, length(y)), y)), tolerance = 1e-5)
  mu_g <- obs[, sum(growth_n * growth) / sum(growth_n)]
  expect_equal(b[["growth"]], obs[, sum(growth_n * (growth - mu_g)^2) / sum(growth_n)], tolerance = 1e-5)
  expect_equal(b[["regeneration"]], obs[, mean((reg - mean(reg))^2)], tolerance = 1e-5)
  # the scaled loss is the pooled loss over the pooled baseline: 1 = no better than the intercept
  k <- c(4, 5, 6)
  expect_equal(as.numeric(auto$history[[1]])[k], unname(as.numeric(raw$history[[1]])[k] / b[k]), tolerance = 1e-5)
})

test_that("(ii) the gradient is the gradient of the pooled per-tree mean", {
  # 4 trees at year 4, 36 at year 10: the old rule gave the 4 trees half the objective
  td <- mk_trees(); obs <- mk_obs()
  m  <- run_fit(obs, td, FAM_TREE, agg = c(mortality = "tree"))   # lr = 0: the grads stay on the parameters
  W  <- m$nn_mortality$parameters[[1]]
  g_fit <- as.numeric(W$grad)
  # the same per-tree predictions, rebuilt with autograd through the private
  # operator the loss uses, then pooled by hand
  pr <- m$.__enclos_env__$private
  ev <- pr$extract_env_method(env)
  envt <- list(growth = torch::torch_tensor(ev$growth_env), mort = torch::torch_tensor(ev$mortality_env))
  terms <- lapply(1:2, function(k) {
    i <- c(4L, 10L)[k]
    pl <- merge(data.table(siteID = 1:Nsites), unique(td[year == i, .(siteID, period_length)]), all.x = TRUE)
    p_t <- torch::torch_tensor(ifelse(is.na(pl$period_length), 0, pl$period_length))
    tt <- torch::torch_arange(1, i)
    M  <- tt$unsqueeze(1)$gt((i - p_t)$unsqueeze(2))$to(dtype = torch::torch_float32()) *
          p_t$gt(0.5)$to(dtype = torch::torch_float32())$unsqueeze(2)
    p <- pr$tree_predictions(m$tree_obs, k, envt, i, M, p_t)$mortality
    y <- m$tree_obs$died[, k, , ]; mask <- y$isnan()$bitwise_not()
    list(sum = torch::nnf_binary_cross_entropy(p[mask]$clamp(1e-6, 1 - 1e-6), y[mask], reduction = "sum"),
         n = as.numeric(mask$sum()))
  })
  grad_of <- function(l) { W$grad$zero_(); l$backward(retain_graph = TRUE); as.numeric(W$grad$clone()) }
  g_pooled <- grad_of((terms[[1]]$sum + terms[[2]]$sum) / (terms[[1]]$n + terms[[2]]$n))
  g_old    <- grad_of(terms[[1]]$sum / terms[[1]]$n + terms[[2]]$sum / terms[[2]]$n)
  expect_equal(c(terms[[1]]$n, terms[[2]]$n), c(4, 36))
  expect_equal(g_fit, g_pooled, tolerance = 1e-4)
  expect_gt(max(abs(g_old - g_pooled)), 10 * max(abs(g_fit - g_pooled)))
})
