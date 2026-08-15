# Integration test for growth_env: (1) runs in the simulator; (2) with k_env = 0
# it reproduces the built-in growth() exactly (strict opt-in extension); (3) with
# k_env != 0 it changes the trajectory, in the expected direction.
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(data.table))

Nsp <- 2L
env_dt <- data.table(expand.grid(siteID = 1:2, year = 1:80))
env_dt[, climate := ifelse(siteID == 1, -1, 1)]   # two contrasting environments

run <- function(growth_func, k_env = NULL) {
  FINN.seed(1)                                     # identical random init across builds
  gp <- if (is.null(k_env))
    createProcess(~1 + climate, growth_func, optimizeSpecies = TRUE, optimizeEnv = TRUE)
  else
    createProcess(~1 + climate, growth_func, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                  custom_parameters = list(k_env = k_env))
  FINN.seed(1)
  m <- finn(
    N_species            = Nsp,
    competition_process  = createProcess(~0, FINN::competition),
    growth_process       = gp,
    regeneration_process = createProcess(~1, FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = createProcess(~1, FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
  )
  sim <- predict(m, init_cohort = NULL, env = env_dt, disturbance = NULL,
                 patches = 20, device = "cpu")
  d <- sim$long$site
  setDT(d)
  d[variable == "ba", .(ba = sum(value)), by = .(siteID, year)]  # per-site total BA
}

base  <- run(FINN::growth)                    # built-in
env0  <- run(FINN::growth_env, rep(0,  Nsp))  # should equal base
env1  <- run(FINN::growth_env, rep(0.05, Nsp))# should differ

cmp <- function(a, b) {
  m <- merge(a, b, by = c("siteID","year"), suffixes = c("_1","_2"))
  max(abs(m$ba_1 - m$ba_2))
}
cat(sprintf("max |BA| diff  base vs growth_env(k_env=0)   : %.3e   (expect ~0)\n", cmp(base, env0)))
cat(sprintf("max |BA| diff  base vs growth_env(k_env=0.05): %.3e   (expect > 0)\n", cmp(base, env1)))
cat("\nfinal-year BA by site:\n")
fin <- Reduce(function(x,y) merge(x,y,by=c("siteID","year")),
              list(base[year==max(year)], env0[year==max(year)], env1[year==max(year)]))
setnames(fin, c("siteID","year","ba_base","ba_kenv0","ba_kenv05")); print(fin)
