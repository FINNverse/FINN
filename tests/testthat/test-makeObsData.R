library(testthat)
library(FINN)
library(data.table)

# makeObsData() derives stand metrics + a closed-cohort mortality pair from a
# raw tree list. It needs no torch, so these run everywhere (incl. CRAN).

test_that("annual rates (the default) over a 10-year interval", {
  res <- suppressMessages(
    makeObsData(toy_tree_dt(), plotsize = 1, aggregate_by_site = FALSE)
  )
  obs  <- as.data.table(res$obs_dt)
  oak  <- obs[year == 2010 & species_name == "oak"]
  pine <- obs[year == 2010 & species_name == "pine"]

  # the interval each row closes, in years; NA at the first inventory
  expect_equal(oak$period_length, 10)
  expect_true(all(is.na(obs[year == 2000, period_length])))
  # growth: the yearly rate that compounds to the increment, 1.2^(1/10) - 1
  expect_equal(oak$growth,  1.2^(1/10) - 1,   tolerance = 1e-10)
  expect_equal(pine$growth, 1.125^(1/10) - 1, tolerance = 1e-10)
  # mortality: 1 of 2 died over 10 years -> 1 - 0.5^(1/10) per year
  expect_equal(oak$mort, 1 - 0.5^(1/10), tolerance = 1e-10)
  expect_equal(pine$mort, 0)
  # the counts and recruitment stay quantities of the interval
  expect_equal(c(oak$n_at_risk, oak$n_died, oak$growth_n), c(2L, 1L, 1L))
  expect_equal(pine$reg, 1)
  # the tree list carries the per-tree rate next to the interval increment
  t1 <- res$tree_dt[treeName == "t1" & year == 2010]
  expect_equal(t1$rel_growth, 0.2, tolerance = 1e-10)
  expect_equal(t1$growth, 1.2^(1/10) - 1, tolerance = 1e-10)

  # site aggregation keeps period_length and derives mort from the pooled counts
  agg <- as.data.table(suppressMessages(
    makeObsData(toy_tree_dt(), plotsize = 1, aggregate_by_site = TRUE, minNyears = NULL))$obs_dt)
  expect_equal(agg[year == 2010 & species_name == "oak", mort], 1 - 0.5^(1/10), tolerance = 1e-10)
  expect_equal(agg[year == 2010 & species_name == "oak", period_length], 10)
})

test_that("growth, mortality counts and regeneration are computed correctly (rates = \"interval\")", {
  res <- suppressMessages(
    makeObsData(toy_tree_dt(), plotsize = 1, aggregate_by_site = FALSE, rates = "interval")
  )
  obs <- as.data.table(res$obs_dt)

  # look at the interval that ENDS in 2010 (the first visit has no "before")
  oak  <- obs[year == 2010 & species_name == "oak"]
  pine <- obs[year == 2010 & species_name == "pine"]

  # --- growth: mean relative diameter increment over the living trees ---------
  # oak: only t1 is alive in 2010 -> 24/20 - 1 = 0.2, resting on 1 tree
  expect_equal(oak$growth, 0.2, tolerance = 1e-8)
  expect_equal(oak$growth_n, 1L)
  # pine: t3 -> 45/40 - 1 = 0.125
  expect_equal(pine$growth, 0.125, tolerance = 1e-8)

  # --- closed-cohort mortality: at-risk = alive at interval START -------------
  # oak: t1 + t2 at risk, t2 dies  -> 1/2
  expect_equal(oak$n_at_risk, 2L)
  expect_equal(oak$n_died,    1L)
  expect_equal(oak$mort,      0.5)
  # pine: only t3 at risk (t4 recruited mid-interval, so not in the cohort)
  expect_equal(pine$n_at_risk, 1L)
  expect_equal(pine$n_died,    0L)
  expect_equal(pine$mort,      0)

  # --- regeneration: "new" trees per plot area --------------------------------
  expect_equal(pine$reg, 1)   # one recruit, plotsize = 1
  expect_equal(oak$reg,  0)

  # --- living tree counts -----------------------------------------------------
  expect_equal(oak$trees,  1)  # t2 is dead
  expect_equal(pine$trees, 2)  # t3 + t4

  # period_length is returned in both modes
  expect_equal(oak$period_length, 10)
})

test_that("the two modes agree where the rates cannot differ", {
  a <- as.data.table(suppressMessages(makeObsData(toy_tree_dt(), plotsize = 1, aggregate_by_site = FALSE))$obs_dt)
  i <- as.data.table(suppressMessages(makeObsData(toy_tree_dt(), plotsize = 1, aggregate_by_site = FALSE,
                                                  rates = "interval"))$obs_dt)
  for (v in c("ba", "dbh", "trees", "growth_n", "n_at_risk", "n_died", "reg", "period_length"))
    expect_identical(a[[v]], i[[v]], info = v)
  # annual = interval compounded back over period_length
  ok <- is.finite(a$growth)
  expect_equal((1 + a$growth[ok])^a$period_length[ok] - 1, i$growth[ok], tolerance = 1e-10)
  ok <- is.finite(a$mort)
  expect_equal(1 - (1 - a$mort[ok])^a$period_length[ok], i$mort[ok], tolerance = 1e-10)
})

test_that("n_died never exceeds n_at_risk (closed-cohort invariant)", {
  res <- suppressMessages(
    makeObsData(toy_tree_dt(), plotsize = 1, aggregate_by_site = FALSE)
  )
  obs <- as.data.table(res$obs_dt)
  expect_true(all(obs$n_died <= obs$n_at_risk))
  # mort is NA exactly where nothing was at risk, a rate in [0,1] otherwise
  expect_true(all(is.na(obs$mort) | (obs$mort >= 0 & obs$mort <= 1)))
  expect_true(all(obs[n_at_risk == 0, is.na(mort)]))
})

test_that("Nspecies caps species and lumps the rest into 'other'", {
  # keep only the single most abundant species; the other becomes 'other'.
  res <- suppressMessages(
    makeObsData(toy_tree_dt(), plotsize = 1, aggregate_by_site = FALSE, Nspecies = 1)
  )
  obs <- as.data.table(res$obs_dt)
  # with Nspecies = 1 we expect the top species plus 'other' (2 categories)
  expect_true("other" %in% obs$species_name)
  expect_lte(uniqueN(obs$species_name), 2L)
})

test_that("lumping into 'other' keeps the mortality counts of the lumped species (FINN#73)", {
  # a second pine makes pine the top species, so oak (with its one death) is lumped
  td <- rbind(toy_tree_dt(), data.table(siteName = "s1", patchName = "p1", treeName = "t5",
                                        year = c(2000, 2010), species_name = "pine", dbh = c(25, 27),
                                        status = "alive", living = TRUE))
  full <- as.data.table(suppressMessages(makeObsData(td, plotsize = 1, aggregate_by_site = FALSE))$obs_dt)
  lump <- as.data.table(suppressMessages(makeObsData(td, plotsize = 1, aggregate_by_site = FALSE, Nspecies = 1))$obs_dt)
  expect_equal(lump[year == 2010 & species_name == "other", c(n_at_risk, n_died)],
               full[year == 2010 & species_name == "oak", c(n_at_risk, n_died)])
  expect_equal(lump[year == 2010 & species_name == "other", mort], 1 - 0.5^(1/10), tolerance = 1e-10)
  expect_equal(lump[, sum(n_died, na.rm = TRUE)], full[, sum(n_died, na.rm = TRUE)])
})

test_that("minNyears drops sites with too few inventories", {
  s2 <- data.table(siteName = "s2", patchName = "p1", treeName = "x1",
                   year = 2005, species_name = "oak", dbh = 15,
                   status = "alive", living = TRUE)   # single inventory -> dropped
  res <- suppressWarnings(suppressMessages(
    makeObsData(rbind(toy_tree_dt(), s2), plotsize = 1,
                aggregate_by_site = FALSE, minNyears = 2)
  ))
  sites <- unique(as.data.table(res$obs_dt)$siteName)
  expect_true("s1"  %in% sites)   # two inventories -> kept
  expect_false("s2" %in% sites)   # one inventory  -> dropped
})
