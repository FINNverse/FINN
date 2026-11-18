# aggregate_species(): the stacked, single-scatter form of aggregate_results()
# used in forward()'s time loop. It must give the same per-site, per-species
# sums as aggregate_results() on fresh zero accumulators, for one and for
# several samples, on several patches, and keep the gradient path.

test_that("aggregate_species equals aggregate_results on zero accumulators", {
  skip_if_no_torch()
  torch::torch_manual_seed(1)
  sites <- 5L; patches <- 3L; cohorts <- 7L; sp <- 4L
  species <- torch::torch_randint(1, sp + 1, c(sites, patches, cohorts), dtype = torch::torch_int64())
  samples <- list(torch::torch_rand(sites, patches, cohorts), torch::torch_rand(sites, patches, cohorts),
                  torch::torch_rand(sites, patches, cohorts))
  ref <- aggregate_results(species, samples, replicate(3, torch::torch_zeros(sites, sp), simplify = FALSE))
  got <- FINN:::aggregate_species(species, samples, sp)
  expect_length(got, 3L)
  for (v in 1:3) {
    expect_equal(got[[v]]$shape, c(sites, sp))
    expect_identical(torch::as_array(got[[v]]), torch::as_array(ref[[v]]))
  }
  # one sample takes the direct path
  got1 <- FINN:::aggregate_species(species, samples[1], sp)
  expect_identical(torch::as_array(got1[[1]]), torch::as_array(ref[[1]]))
  # a species with no cohort sums to 0, and the row sums are the cohort sums
  expect_equal(sum(torch::as_array(got[[2]])), sum(torch::as_array(samples[[2]])), tolerance = 1e-5)
})

test_that("aggregate_species keeps the gradient", {
  skip_if_no_torch()
  sites <- 3L; patches <- 2L; cohorts <- 4L; sp <- 3L
  species <- torch::torch_randint(1, sp + 1, c(sites, patches, cohorts), dtype = torch::torch_int64())
  x <- torch::torch_rand(sites, patches, cohorts, requires_grad = TRUE)
  out <- FINN:::aggregate_species(species, list(x * 2, x), sp)
  (out[[1]]$sum() + out[[2]]$sum())$backward()
  # d/dx of sum(2x) + sum(x) = 3 everywhere
  expect_equal(torch::as_array(x$grad), array(3, dim = c(sites, patches, cohorts)))
})
