# Fit FINN to ONE species' NW-FVA yield table. arg: species [epochs].
# Writes dev/pilots/allsp_<species>.rds. Run 5 in parallel for the all-species set.
args   <- commandArgs(trailingOnly = TRUE)
sp     <- args[1]
epochs <- if (length(args) >= 2) as.integer(args[2]) else 200L
lr     <- if (length(args) >= 3) as.numeric(args[3]) else 0.01
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(data.table))
try(torch::torch_set_num_threads(2L), silent = TRUE)   # avoid oversubscription when parallel

gy <- as.data.table(nwfva_yield_tables())[Art == sp][order(Ekl, Alter)]
classes <- sort(unique(gy$Ekl))
obs <- list(); init <- list(); envm <- list()
for (i in seq_along(classes)) {
  d <- gy[Ekl == classes[i]]; y0 <- d$Alter[1]; d[, yr := Alter - y0]
  obs[[i]]  <- data.table(siteID = i, year = d$yr[-1], species = 1L,
                          dbh = d$Dg[-1], ba = d$G[-1], trees = d$N[-1],
                          growth = NA_real_, mort = NA_real_, reg = NA_real_)
  init[[i]] <- data.table(siteID = i, patchID = 1L, species = 1L, dbh = d$Dg[1], trees = d$N[1])
  envm[[i]] <- data.table(siteID = i, ekl = as.numeric(classes[i]))
}
obs_dt <- rbindlist(obs); Tmax <- max(obs_dt$year)
env_dt <- merge(CJ(siteID = seq_along(classes), year = 1:Tmax), rbindlist(envm), by = "siteID")
init_cohort <- makeInitCohorts(rbindlist(init), Nspecies = 1L)

FINN.seed(1)
m <- finn(N_species = 1L,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(~ ekl, FINN::growth_env, optimizeSpecies = TRUE,
                                       optimizeEnv = TRUE, custom_parameters = list(k_env = 0)),
  mortality_process    = createProcess(~ ekl, FINN::mortality, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~1, FINN::regeneration, initEnv = list(matrix(-12, 1, 1)),
                                       optimizeSpecies = FALSE, optimizeEnv = FALSE))
fit(m, data = obs_dt, env = env_dt, init_cohort = init_cohort,
    patches = 1, patch_size = 1.0, device = "cpu", epochs = epochs, lr = lr,
    batchsize = length(classes),
    weights = c(dbh = 5, ba = 1, trees = 1, growth = 0, mortality = 0, regeneration = 0),
    env_autoscale = TRUE, plot_progress = FALSE)

sim <- predict(m, init_cohort = init_cohort, env = env_dt, patches = 1, patch_size = 1.0, device = "cpu")
s <- as.data.table(sim$long$site)
cmp <- merge(obs_dt[, .(siteID, year, Dg_obs = dbh, N_obs = trees)],
             s[variable == "dbh", .(siteID, year, Dg_sim = value)], by = c("siteID","year"))
cmp <- merge(cmp, s[variable == "trees", .(siteID, year, N_sim = value)], by = c("siteID","year"))
amin <- gy[, .(a = min(Alter)), Ekl]
cmp[, Ekl := classes[siteID]][, age := year + amin$a[match(Ekl, amin$Ekl)]][, Art := sp]
r2 <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
out <- list(cmp = cmp, dgR2 = r2(cmp$Dg_obs, cmp$Dg_sim), nR2 = r2(cmp$N_obs, cmp$N_sim))
saveRDS(out, sprintf("dev/pilots/allsp_%s.rds", sp))
cat(sprintf("%s: Dg R2=%.3f  N R2=%.3f\n", sp, out$dgR2, out$nR2))
