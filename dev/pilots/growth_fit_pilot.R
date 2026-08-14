# =============================================================================
# Pilot: can FINN's growth response reproduce the NW-FVA yield-table diameter
# growth?  FINN growth (dominant cohort, full light) is
#     g = exp(pred) * exp(-k * dbh)          # g = relative diameter growth / yr
# with pred a linear environmental predictor. We stand site productivity in for
# the environment via the relative yield class Ekl:
#     log g = (b0 + b1 * Ekl) - k * dbh
# -> a LINEAR fit per species. We then INTEGRATE the fitted law year by year and
# compare the reconstructed Dg(age) to the yield table.
#
# This isolates growth from the density/thinning feedback (that needs the
# management operator). It answers: does FINN's growth FORM fit real GY growth?
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

gy <- read.csv("inst/extdata/nwfva_gy_tables.csv", stringsAsFactors = FALSE)
gy <- gy[order(gy$Art, gy$Ekl, gy$Alter), ]
gy$g_obs <- gy$dDg_dt / gy$Dg                       # observed relative growth (1/yr)
fitdat <- gy[is.finite(gy$g_obs) & gy$g_obs > 0, ]

species <- unique(gy$Art)
coefs <- list(); recon <- list()

for (s in species) {
  d <- fitdat[fitdat$Art == s, ]
  # FINN growth form, linear in (b0, b1, k):
  m <- lm(log(g_obs) ~ Ekl + Dg, data = d)
  b0 <- coef(m)[["(Intercept)"]]; b1 <- coef(m)[["Ekl"]]; k <- -coef(m)[["Dg"]]
  gpred <- function(dbh, ekl) exp(b0 + b1 * ekl) * exp(-k * dbh)
  coefs[[s]] <- data.frame(Art = s, b0 = b0, b1 = b1, k = k,
                           R2_logg = summary(m)$r.squared)

  # Reconstruct Dg(age) by integrating the fitted growth from each class's first
  # tabulated point, stepping 1 yr at a time (dbh_{t+1} = dbh_t * (1 + g)).
  for (e in unique(d$Ekl)) {
    tab <- gy[gy$Art == s & gy$Ekl == e & is.finite(gy$Dg), c("Alter", "Dg")]
    if (nrow(tab) < 3) next
    age <- tab$Alter; dbh <- numeric(length(age)); dbh[1] <- tab$Dg[1]
    cur <- tab$Dg[1]
    for (i in 2:length(age)) {
      for (yr in seq_len(age[i] - age[i - 1])) cur <- cur * (1 + gpred(cur, e))
      dbh[i] <- cur
    }
    recon[[paste(s, e)]] <- data.frame(Art = s, Ekl = e, Alter = age,
                                       Dg_obs = tab$Dg, Dg_fit = dbh)
  }
}

coefs <- do.call(rbind, coefs); rownames(coefs) <- NULL
recon <- do.call(rbind, recon)

# fit quality on the reconstructed diameter trajectory
ss <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
r2_traj <- tapply(seq_len(nrow(recon)), recon$Art,
                  function(i) ss(recon$Dg_obs[i], recon$Dg_fit[i]))

cat("=== fitted FINN growth parameters per species ===\n")
print(cbind(coefs, R2_Dg_traj = round(r2_traj[coefs$Art], 3)), digits = 3)

# --- plot: observed (points) vs FINN-growth reconstruction (lines), by class ---
recon$Ekl <- factor(recon$Ekl)
p <- ggplot(recon, aes(Alter, colour = Ekl)) +
  geom_point(aes(y = Dg_obs), size = 1.1, alpha = 0.8) +
  geom_line(aes(y = Dg_fit), linewidth = 0.7) +
  facet_wrap(~Art, scales = "free") +
  scale_colour_viridis_d(option = "D", end = 0.9, name = "Ertragsklasse") +
  labs(title = "Pilot: FINN growth fitted to NW-FVA yield tables",
       subtitle = "points = yield-table Dg(age); lines = FINN growth law integrated from the first point",
       x = "Alter (yr)", y = expression(D[g]~"(cm)")) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

dir.create("dev/pilots", showWarnings = FALSE, recursive = TRUE)
ggsave("dev/pilots/growth_fit_pilot.png", p, width = 9, height = 6, dpi = 130)
write.csv(coefs, "dev/pilots/growth_fit_coefs.csv", row.names = FALSE)
cat("\nwrote dev/pilots/growth_fit_pilot.png\n")
