library(testthat)
library(FINN)
library(data.table)

## fia_tree_dt / fia_tree_test: one row per tree alive at the start of a
## remeasurement interval, the trees behind the closed-cohort counts of
## fia_obs_dt / fia_obs_test. No torch needed.

tree_tables <- list(train = list(tree = fia_tree_dt, obs = fia_obs_dt, init = fia_init_trees),
                    holdout = list(tree = fia_tree_test, obs = fia_obs_test, init = fia_init_test))

test_that("the tree tables have the columns fit(tree_data = ) and predictTrees() need", {
  for (s in tree_tables) {
    expect_identical(names(s$tree), c("siteID", "year", "patchID", "species", "dbh", "growth", "died", "period_length"))
    expect_true(all(s$tree$died %in% c(0, 1)))
    expect_true(all(s$tree$dbh > 0))
    expect_true(all(is.na(s$tree$growth[s$tree$died == 1])))     # a dead tree has no growth rate
    expect_true(all(is.finite(s$tree$growth[s$tree$died == 0])))
    expect_true(all(s$tree$period_length >= 10))
  }
})

test_that("the tree tables join to the site tables on siteID, year and species", {
  for (s in tree_tables) {
    key_obs  <- unique(s$obs[, .(siteID, year, species)])
    key_tree <- unique(s$tree[, .(siteID, year, species)])
    expect_equal(nrow(merge(key_tree, key_obs, by = c("siteID", "year", "species"))), nrow(key_tree))
    # the sites of the site table; a site without a tree at risk in any interval
    # (two holdout sites whose first inventory found only dead trees) has no row
    expect_true(all(unique(s$tree$siteID) %in% unique(s$obs$siteID)))
    expect_setequal(unique(s$tree$siteID), s$obs[n_at_risk > 0, unique(siteID)])
    expect_true(all(s$tree$species %in% fia_species_dt$species))
    # the patch ids of the initial tree list (a patch may hold its first living tree
    # only after year 0), and the same window length as the site table
    expect_true(all(s$tree$patchID %in% unique(s$init$patchID)))
    pl <- merge(unique(s$tree[, .(siteID, year, period_length)]),
                unique(s$obs[, .(siteID, year, period_length)]), by = c("siteID", "year"))
    expect_equal(pl$period_length.x, pl$period_length.y)
  }
})

test_that("the deaths per site, year and species are the n_died of the site table (open-cohort tolerance)", {
  for (s in tree_tables) {
    cnt <- merge(s$obs[, .(siteID, year, species, n_at_risk, n_died)],
                 s$tree[, .(n = .N, d = sum(died)), by = .(siteID, year, species)],
                 by = c("siteID", "year", "species"), all.x = TRUE)
    cnt[is.na(n), `:=`(n = 0L, d = 0)]
    # the tree list covers every cohort of the site table ...
    expect_true(all(cnt$n >= cnt$n_at_risk))
    expect_true(all(cnt$d >= cnt$n_died))
    # ... and agrees with it in nearly every cell: the site table lacks the rows
    # of a species whose trees had all died by a census, a few dozen trees
    expect_gt(mean(cnt$n == cnt$n_at_risk & cnt$d == cnt$n_died), 0.99)
    expect_lt(sum(cnt$d - cnt$n_died), 10)
    expect_lt(sum(cnt$n - cnt$n_at_risk), 100)
  }
})
