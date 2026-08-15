# Integration test for the simulate_managed() driver: it segments the horizon at
# the schedule years, thins the cohort snapshot per patch, re-initialises from the
# thinned (rectangularised) state, and logs the harvest. Toy random-parameter model
# -- this checks the MECHANISM (a realistic managed-vs-unmanaged contrast comes with
# the yield-table fit).
suppressMessages(devtools::load_all(".", quiet = TRUE)); suppressMessages(library(data.table))

env <- data.table(expand.grid(siteID = 1L, year = 1:60)); env[, climate := 0]
FINN.seed(1)
m <- finn(N_species = 1L,
  competition_process  = createProcess(~0, FINN::competition),
  growth_process       = createProcess(~1, FINN::growth,       optimizeSpecies = TRUE),
  regeneration_process = createProcess(~1, FINN::regeneration, optimizeSpecies = TRUE),
  mortality_process    = createProcess(~1, FINN::mortality,    optimizeSpecies = TRUE))

sched <- data.frame(year = c(20, 35, 50))
sched$prescription <- I(list(thin_from_above(0.3), thin_from_above(0.3), thin_from_above(0.3)))
run <- simulate_managed(m, env = env, patches = 25, schedule = sched,
                        patch_size = 0.1, device = "cpu")

stopifnot(nrow(run$harvest) == 3L,                       # a thinning at each entry
          all(run$harvest$year == c(20, 35, 50)),
          all(run$harvest$removed_trees_ha > 0),
          setequal(unique(run$trajectory$year), 1:60))   # full stitched horizon
cat("simulate_managed(): OK -- 3 thinnings logged, trajectory spans years 1-60\n")
print(run$harvest)
