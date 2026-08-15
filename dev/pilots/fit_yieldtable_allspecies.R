# =============================================================================
# Fit FINN to the NW-FVA yield tables for ALL FIVE species (one fit per species,
# each the proven Fichte setup: yield classes = sites, site index = env,
# growth_env + mortality-absorbs-removal, regeneration off). Reports Dg / N R2 per
# species and a faceted obs-vs-fitted figure.
# =============================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages({library(data.table); library(ggplot2)})

gy_all <- as.data.table(nwfva_yield_tables())
species <- c("Eiche", "Buche", "Fichte", "Douglasie", "Kiefer")
r2 <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)

fit_species <- function(sp, epochs = 400L) {
  gy <- gy_all[Art == sp][order(Ekl, Alter)]
  classes <- sort(unique(gy$Ekl))
  obs <- list(); init <- list(); envm <- list()
  for (i in seq_along(classes)) {
    d <- gy[Ekl == classes[i]]; y0 <- d$Alter[1]; d[, yr := Alter - y0]
    obs[[i]]  <- data.table(siteID = i, year = d$yr[-1], species = 1L,
                            dbh = d$Dg[-1], ba = d$G[-1], trees = d$N[-1],
                            growth = NA_real_, mort = NA_real_, reg = NA_real_)
    init[[i]] <- data.table(siteID = i, patchID = 1L, species = 1L,
                            dbh = d$Dg[1], trees = d$N[1])
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
      patches = 1, patch_size = 1.0, device = "cpu", epochs = epochs, lr = 0.02,
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
  list(cmp = cmp, dgR2 = r2(cmp$Dg_obs, cmp$Dg_sim), nR2 = r2(cmp$N_obs, cmp$N_sim))
}

res <- lapply(species, function(sp) { message("fitting ", sp, " ..."); fit_species(sp) })
names(res) <- species
tab <- data.table(Art = species,
                  Dg_R2 = round(sapply(res, `[[`, "dgR2"), 3),
                  N_R2  = round(sapply(res, `[[`, "nR2"), 3))
cat("\n=== ALL-SPECIES FIT QUALITY ===\n"); print(tab)

cmp <- rbindlist(lapply(res, `[[`, "cmp"))
cmp$Ekl <- factor(cmp$Ekl)
p <- ggplot(cmp, aes(age, colour = Ekl)) +
  geom_point(aes(y = Dg_obs), size = .7, alpha = .8) +
  geom_line(aes(y = Dg_sim), linewidth = .6) +
  facet_wrap(~Art, scales = "free") +
  scale_colour_viridis_d(option = "D", end = .9, name = "Ertragsklasse") +
  labs(title = "FINN fitted to the NW-FVA yield tables, all five species",
       subtitle = paste0("points = table Dg(age); lines = fitted FINN. Dg R2: ",
                         paste(sprintf("%s %.3f", species, tab$Dg_R2), collapse = "  ")),
       x = "Alter (yr)", y = expression(D[g]~"(cm)")) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave("dev/pilots/fit_yieldtable_allspecies.png", p, width = 10, height = 6.2, dpi = 130)
fwrite(tab, "dev/pilots/fit_yieldtable_allspecies_r2.csv")
cat("wrote dev/pilots/fit_yieldtable_allspecies.png\n")
