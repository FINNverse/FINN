# ============================================================================
# Path B -- analyse the three fitted conditions
# ============================================================================
# Loads pathB_{free,ruger,learned}.rds and answers:
#   PRIMARY   : rare-species prediction gain -- held-out error by rare/common
#               for free vs ruger vs learned.
#   SECONDARY : does the LEARNED membership match Rueger's PFTs? (ARI + confusion)
# Run locally after the cluster fits (results/ synced back).
# ============================================================================

local({ l <- file.path(path.expand("~"), "Rlib-finn020"); if (dir.exists(l)) .libPaths(c(l, .libPaths())) })
suppressMessages({ library(data.table); library(ggplot2) })
D <- "dev/pft-bci"; RES <- file.path(D, "results")
TAG <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else ""   # e.g. "" or "_smoke"
conds <- c("free","ruger","learned")
vars  <- c("ba","dbh","trees","growth","mort","reg")

load1 <- function(c) { f <- file.path(RES, paste0("pathB_", c, TAG, ".rds"))
  if (!file.exists(f)) stop("missing ", f, " -- run 03_pathB_fit.R COND=", c); readRDS(f) }
R <- lapply(conds, load1); names(R) <- conds

# rare/common split from the (shared) training-obs counts
nobs <- R[[1]]$nobs
rare_sp <- nobs[n_train < 20]$species
cat("rare species (<20 train growth obs):", length(rare_sp),
    " | common:", nrow(nobs) - length(rare_sp), "\n\n")

# ---- held-out error per condition ------------------------------------------
held_err <- function(res) {
  pred <- as.data.table(res$pred)                       # siteID, year, species, variable, value
  obs  <- melt(res$obs_te, id.vars = c("siteID","year","species"),
               measure.vars = vars, variable.name = "variable", value.name = "obs")
  obs[, variable := as.character(variable)]
  m <- merge(pred[variable %in% vars], obs, by = c("siteID","year","species","variable"))
  m <- m[is.finite(obs) & is.finite(value)]
  # normalised RMSE per (species, variable), then mean across variables per species
  m[, nrmse := (value - obs)^2]
  ps <- m[, .(rmse = sqrt(mean(nrmse))/ (sd(obs) + 1e-6)), by = .(species, variable)]
  ps[, .(err = mean(rmse, na.rm = TRUE)), by = species]
}
tab <- rbindlist(lapply(conds, function(c) { e <- held_err(R[[c]]); e[, cond := c]; e }))
tab[, grp := ifelse(species %in% rare_sp, "rare", "common")]

summ <- dcast(tab[, .(err = mean(err, na.rm = TRUE)), by = .(cond, grp)],
              grp ~ cond, value.var = "err")
cat("==================== held-out prediction error (lower = better) =========\n")
cat("mean normalised RMSE across variables, by species group:\n\n")
print(summ, digits = 3)
cat("\nRare-species gain (free -> learned):",
    round(100*(1 - summ[grp=="rare", learned]/summ[grp=="rare", free]), 1), "%\n")
cat("Learned vs Rueger on rare species:",
    round(100*(1 - summ[grp=="rare", learned]/summ[grp=="rare", ruger]), 1), "% better\n")

# ---- secondary: learned membership vs Rueger PFTs --------------------------
ari <- function(a,b){tab<-table(a,b);n<-sum(tab);si<-sum(choose(rowSums(tab),2));sj<-sum(choose(colSums(tab),2));sij<-sum(choose(tab,2));e<-si*sj/choose(n,2);(sij-e)/((si+sj)/2-e)}
A <- R$learned$A; ruger <- R$learned$ruger
pft_sp <- which(ruger <= max(ruger) - 1L)                # classified (non-other) species
cl <- max.col(A[pft_sp, ])
cat("\n---- learned membership vs Rueger PFTs (secondary) ----\n")
cat("Adjusted Rand Index:", round(ari(cl, ruger[pft_sp]), 3), "\n")
print(table(Rueger = ruger[pft_sp], learned = cl))

# ---- figure ----------------------------------------------------------------
pl <- melt(summ, id.vars = "grp", variable.name = "cond", value.name = "err")
pl$cond <- factor(pl$cond, levels = c("free","ruger","learned"))
p <- ggplot(pl, aes(cond, err, fill = grp)) +
  geom_col(position = position_dodge()) +
  labs(title = "BCI Path B: held-out prediction error by species group",
       subtitle = "free (per-species) vs Rueger PFTs vs learned membership",
       x = NULL, y = "held-out normalised RMSE", fill = "species") +
  theme_minimal(base_size = 12)
ggsave(file.path(D, paste0("pathB_gain", TAG, ".png")), p, width = 7, height = 4.2, dpi = 130)
cat("\nfigure:", file.path(D, paste0("pathB_gain", TAG, ".png")), "\n")
