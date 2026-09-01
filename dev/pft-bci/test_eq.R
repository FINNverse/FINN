# quick de-risk: load a checkpoint model + run a small equilibrium sim
local({ l <- file.path(path.expand("~"), "Rlib-finn020"); if (dir.exists(l)) .libPaths(c(l, .libPaths())) })
suppressWarnings(suppressMessages({ library(data.table); library(FINN) }))
source("dev/pft-bci/finn_membership.R")
D <- "dev/pft-bci"
m <- torch::torch_load(file.path(D, "results/ckpt_learned/epoch_1000model.pt"))
cat("loaded:", paste(class(m)[1]), " N_species:", m$N_species, " K:", m$mm_K, "\n")
env <- fread(file.path(D, "data/env.csv")); coh <- fread(file.path(D, "data/initial_cohorts1985.csv"))
ny <- max(env$year)
env_long <- rbindlist(lapply(0:3, function(k) { e <- copy(env); e[, year := year + k*ny]; e }))
d <- unique(env_long[, .(siteID, year)]); set.seed(1)
d[, intensity := rbinom(.N, 1, 1) * runif(.N, 0.0043, 0.016)]
ch <- FINN::CohortMat(obs_df = coh, sp = m$N_species)
m$eval()
t0 <- Sys.time()
sim <- m$simulate(env = env_long, disturbance = d, init_cohort = ch, patches = 4L, device = "cpu")
w <- as.data.table(sim$wide$site)
cat("sim OK: rows", nrow(w), " vars:", paste(intersect(c("dbh","ba","trees","reg","mort"), names(w)), collapse=","),
    " maxyear", max(w$year), " secs", round(as.numeric(difftime(Sys.time(), t0, units="secs")),1), "\n")
eq <- w[year > max(year) - 30, .(ba = mean(ba, na.rm=TRUE)), by = species]
cat("equilibrium ba per species: n", nrow(eq), " range", round(min(eq$ba),3), "-", round(max(eq$ba),3), "\n")
