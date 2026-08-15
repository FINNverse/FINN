# Aggregate the 5 parallel per-species yield-table fits into one facet figure + R2 table.
suppressMessages({library(data.table); library(ggplot2)})
sp <- c("Eiche", "Buche", "Fichte", "Douglasie", "Kiefer")
res <- lapply(sp, function(s) readRDS(sprintf("dev/pilots/allsp_%s.rds", s)))
names(res) <- sp

tab <- data.table(Art = sp,
                  Dg_R2 = round(sapply(res, `[[`, "dgR2"), 3),
                  N_R2  = round(sapply(res, `[[`, "nR2"), 3))
cat("=== ALL-SPECIES FIT QUALITY ===\n"); print(tab)
fwrite(tab, "dev/pilots/fit_yieldtable_allspecies_r2.csv")

good <- tab$Art[tab$Dg_R2 > 0.9]                       # converged species
bad  <- tab$Art[tab$Dg_R2 <= 0.9]
cmp <- rbindlist(lapply(res[good], `[[`, "cmp"))
cmp$Ekl <- factor(cmp$Ekl)
cmp[, Art := factor(Art, levels = good)]
sub <- paste0("points = table Dg(age); lines = fitted FINN.  Converged: ",
              paste(sprintf("%s R2=%.3f", good, tab$Dg_R2[match(good, tab$Art)]), collapse = ", "),
              ".  Diverged (mechanistic fit, degenerate optimum): ",
              paste(bad, collapse = ", "), ".")
p <- ggplot(cmp, aes(age, colour = Ekl)) +
  geom_point(aes(y = Dg_obs), size = .8, alpha = .8) +
  geom_line(aes(y = Dg_sim), linewidth = .65) +
  facet_wrap(~Art, scales = "free") +
  scale_colour_viridis_d(option = "D", end = .9, name = "Ertragsklasse") +
  labs(title = "FINN fitted to the NW-FVA yield tables (through fit())",
       subtitle = sub, x = "Alter (yr)", y = expression(D[g]~"(cm)")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 8.5))
ggsave("dev/pilots/fit_yieldtable_allspecies.png", p, width = 9, height = 4.2, dpi = 130)
cat("wrote dev/pilots/fit_yieldtable_allspecies.png\n")
