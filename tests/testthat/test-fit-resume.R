library(testthat)
library(FINN)
library(data.table)

# Resuming training on a model restored by torch_load(). The optimizer stored
# in the model does not survive saving: an optim_ignite_* optimizer comes back
# as a dead external pointer (fit() used to crash on it), a pure-R one holds
# copies of the old tensors (fit() used to run but train nothing). fit() must
# rebuild it; within one session it must keep reusing it (Adam state).
# Torch-gated; tiny FIA subset.

resume_setup <- function() {
  ext <- function(f) system.file("extdata", f, package = "FINN")
  obs <- fread(ext("fia_obs_dt.csv")); env <- fread(ext("fia_env_dt.csv")); it <- fread(ext("fia_init_trees.csv"))
  sites <- head(sort(unique(obs$siteID)), 4L)
  Nsp <- max(obs$species)
  list(obs = obs[siteID %in% sites], env = env[siteID %in% sites],
       init = makeInitCohorts(it[siteID %in% sites], Nspecies = Nsp), Nsp = Nsp)
}
resume_model <- function(Nsp) {
  FINN.seed(1)
  finn(N_species = Nsp, recruits_dbh = 12.9,
       competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
       growth_process       = createProcess(~ temp + prec, FINN::growth,       optimizeSpecies = TRUE, optimizeEnv = TRUE),
       mortality_process    = createProcess(~ temp + prec, FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE),
       regeneration_process = createProcess(~ temp + prec, FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE))
}
resume_fit <- function(m, d, opt) {
  suppressMessages(fit(m, env = d$env, data = d$obs, init_cohort = d$init, device = "cpu",
                       epochs = 1L, patches = 2, patch_size = 0.06, lr = 0.02,
                       optimizer = opt, plot_progress = FALSE))
}
param_values <- function(m) unlist(lapply(m$parameters, function(p) as.numeric(p$cpu())))

for (opt_name in c("optim_ignite_adam", "optim_adam")) {
  test_that(sprintf("fit() resumes a torch_load()ed model with %s at the same lr", opt_name), {
    skip_if_no_torch()
    skip_if_not(exists(opt_name, envir = asNamespace("torch")), paste("torch has no", opt_name))
    opt <- getExportedValue("torch", opt_name)
    d <- resume_setup()
    m <- resume_model(d$Nsp)
    resume_fit(m, d, opt)

    f <- tempfile(fileext = ".pt"); on.exit(unlink(f), add = TRUE)
    torch::torch_save(m, f)
    m2 <- torch::torch_load(f)
    before <- param_values(m2)
    expect_no_error(resume_fit(m2, d, opt))                       # ignite: used to crash here
    after <- param_values(m2)
    expect_true(all(is.finite(after)))
    expect_gt(max(abs(after - before)), 0)                        # pure R: used to train nothing
    # and the rebuilt optimizer drives the loaded model's own tensors
    held <- m2$optimizer$param_groups[[1]]$params
    expect_true(identical(held[[1]], m2$parameters[[1]]))
  })
}

test_that("a second fit() in the same session keeps the optimizer (and its Adam state)", {
  skip_if_no_torch()
  d <- resume_setup()
  m <- resume_model(d$Nsp)
  opt <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
  resume_fit(m, d, opt)
  first <- m$optimizer
  resume_fit(m, d, opt)
  expect_true(identical(m$optimizer, first))
})
