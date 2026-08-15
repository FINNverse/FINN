# =============================================================================
# Pilot 2: does the growth misfit come from ignoring LIGHT/DENSITY (the missing
# self-thinning / graduated thinning), and does an ENVIRONMENT-controlled
# size-decline help?
#
# FINN growth is really  g = shade(light) * exp(pred_env) * exp(-k*dbh).
# Pilot 1 set shade=1 (full light). Here we bring the stand's basal area G back
# in as the competition/light signal (G carries the thinning: it dips when the
# stand is thinned), and let the size-decline k depend on site (Ekl:dbh), with a
# species baseline. Nested models, same species, compared on the reconstructed
# Dg(age):
#   A  log g ~ Ekl + dbh                 (Pilot 1: species k, no density)
#   B  log g ~ Ekl + dbh + G             (+ density/light)
#   C  log g ~ Ekl + dbh + Ekl:dbh       (env-controlled size-decline)
#   D  log g ~ Ekl + dbh + G + Ekl:dbh   (both)
# Reconstruction integrates the fitted law yr-by-yr, taking G(age) from the table
# (linearly interpolated) -- i.e. giving the growth the real managed density.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

gy <- read.csv("inst/extdata/nwfva_gy_tables.csv", stringsAsFactors = FALSE)
gy <- gy[order(gy$Art, gy$Ekl, gy$Alter), ]
gy$g_obs <- gy$dDg_dt / gy$Dg
fitdat <- gy[is.finite(gy$g_obs) & gy$g_obs > 0, ]
species <- unique(gy$Art)

forms <- list(
  A = log(g_obs) ~ Ekl + Dg,
  B = log(g_obs) ~ Ekl + Dg + G,
  C = log(g_obs) ~ Ekl + Dg + Ekl:Dg,
  D = log(g_obs) ~ Ekl + Dg + G + Ekl:Dg
)

ss <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)

# reconstruct Dg(age) for one class under a fitted lm `m`, using table G(age)
recon_traj <- function(m, sub) {
  age <- sub$Alter; Gv <- sub$G; ekl <- sub$Ekl[1]
  dbh <- numeric(length(age)); dbh[1] <- sub$Dg[1]; cur <- sub$Dg[1]
  for (i in 2:length(age)) {
    for (yr in seq_len(age[i] - age[i - 1])) {
      frac <- (age[i - 1] + yr - age[i - 1]) / (age[i] - age[i - 1])
      Gnow <- Gv[i - 1] + frac * (Gv[i] - Gv[i - 1])         # interpolate G within step
      nd <- data.frame(Ekl = ekl, Dg = cur, G = Gnow)
      cur <- cur * (1 + exp(predict(m, nd)))
    }
    dbh[i] <- cur
  }
  data.frame(Art = sub$Art[1], Ekl = ekl, Alter = age, Dg_obs = sub$Dg, Dg_fit = dbh)
}

r2tab <- list(); best_recon <- list()
for (s in species) {
  d  <- fitdat[fitdat$Art == s, ]
  ds <- gy[gy$Art == s & is.finite(gy$Dg), ]
  r2_traj <- numeric(length(forms)); names(r2_traj) <- names(forms)
  rec_by_model <- list()
  for (nm in names(forms)) {
    m <- lm(forms[[nm]], data = d)
    rc <- do.call(rbind, lapply(split(ds, ds$Ekl), function(sub)
      if (nrow(sub) >= 3) recon_traj(m, sub) else NULL))
    rec_by_model[[nm]] <- rc
    r2_traj[nm] <- ss(rc$Dg_obs, rc$Dg_fit)
  }
  r2tab[[s]] <- data.frame(Art = s, t(round(r2_traj, 3)))
  best <- names(which.max(r2_traj))
  bl <- rec_by_model[[best]]; bl$model <- best
  best_recon[[s]] <- bl
}

r2tab <- do.call(rbind, r2tab); rownames(r2tab) <- NULL
cat("=== R2 of reconstructed Dg(age), by model ===\n"); print(r2tab)
cat("\n(A = Pilot 1; B +density G; C env-controlled k; D both)\n")

best_recon <- do.call(rbind, best_recon)
best_recon$Ekl <- factor(best_recon$Ekl)
lab <- unique(best_recon[, c("Art", "model")])
lab_txt <- setNames(sprintf("%s  (best: %s)", lab$Art, lab$model), lab$Art)

p <- ggplot(best_recon, aes(Alter, colour = Ekl)) +
  geom_point(aes(y = Dg_obs), size = 1.1, alpha = 0.8) +
  geom_line(aes(y = Dg_fit), linewidth = 0.7) +
  facet_wrap(~Art, scales = "free", labeller = as_labeller(lab_txt)) +
  scale_colour_viridis_d(option = "D", end = 0.9, name = "Ertragsklasse") +
  labs(title = "Pilot 2: growth with density (G) + environment-controlled size-decline",
       subtitle = "points = yield-table Dg(age); lines = best model, integrated using the table's managed density G(age)",
       x = "Alter (yr)", y = expression(D[g]~"(cm)")) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

ggsave("dev/pilots/growth_fit_density.png", p, width = 9, height = 6, dpi = 130)
write.csv(r2tab, "dev/pilots/growth_fit_density_r2.csv", row.names = FALSE)
cat("\nwrote dev/pilots/growth_fit_density.png\n")
