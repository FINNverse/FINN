## Release-note figures on the bundled Oregon FIA data: the README model
## (joint fit, 300 epochs) with its 200-plot test set as holdout.
## Run: Rscript demo_fia.R <output folder>   (about 15 minutes, most of it the fit)
suppressMessages({library(FINN); library(data.table); library(torch)})
OUT <- commandArgs(trailingOnly = TRUE)[1]
png_save <- function(p, name, w, h) {
  grDevices::png(file.path(OUT, name), width = w, height = h, units = "in", res = 110, bg = "white")
  print(p); grDevices::dev.off()
}
Nsp  <- max(fia_obs_dt$species)
init <- makeInitCohorts(fia_init_trees, Nspecies = Nsp)
init_test <- makeInitCohorts(fia_init_test, Nspecies = Nsp)
FINN.seed(1)
m <- finn(N_species = Nsp, recruits_dbh = 12.9,
          competition_process  = createProcess(~ 0, competition, optimizeSpecies = TRUE),
          growth_process       = createProcess(~ temp + prec, growth, optimizeSpecies = TRUE),
          mortality_process    = createProcess(~ temp + prec, mortality, optimizeSpecies = TRUE),
          regeneration_process = createProcess(~ temp + prec, regeneration_capped, optimizeSpecies = TRUE))
t0 <- Sys.time()
fit(m, data = fia_obs_dt, env = fia_env_dt, init_cohort = init,
    epochs = 300L, batchsize = 40L, patches = 4, patch_size = 0.06,
    tree_data = fia_tree_dt, plot_progress = FALSE,
    holdout = list(data = fia_obs_test, env = fia_env_test, init_cohort = init_test, tree_data = fia_tree_test))
cat(sprintf("fit: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
torch_save(m, file.path(OUT, "fia_model.pt"))

ev <- m$evaluation
capture.output(print(ev), file = file.path(OUT, "fia_evaluation.txt"))
s <- summary(m)
capture.output(print(s), file = file.path(OUT, "fia_summary.txt"))

pdf(NULL)
png_save(plot(ev), "fia_compare_scatter.png", 11, 7.5)
png_save(plot(ev, type = "heatmap"), "fia_compare_heatmap.png", 13, 5.5)
png_save(plot(m), "fia_plot_model.png", 11, 7)
ale <- ALE(m, plot = FALSE)
png_save(plot(ale), "fia_ale.png", 11, 6)
cold_to_warm <- fia_env_dt[, .(temp = mean(temp)), by = siteID][order(temp)]$siteID[c(10, 100, 190)]
bg <- simulateForest(m, bare_ground = TRUE, years = 500, sites = cold_to_warm)
png_save(plot(bg), "fia_bare_ground.png", 11, 4.5)
sim <- simulateForest(m, env = fia_env_test, init_cohort = init_test)
png_save(plot(sim), "fia_simulation.png", 11, 4)
png_save(comparePlot(sim, fia_obs_test, type = "timeseries", species = 1:3), "fia_timeseries.png", 11, 7)
dev.off()
cat("done\n")
