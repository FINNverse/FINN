# Dig into the annual overshoot: is it a growth -> overshoot -> mortality-crash
# limit cycle? Plot annual ba, mean dbh, growth, mort together for the learned
# annual model, and print the year-by-year values.
local({ l <- file.path(path.expand("~"), "Rlib-finn020"); if (dir.exists(l)) .libPaths(c(l, .libPaths())) })
suppressWarnings(suppressMessages({ library(data.table); library(FINN); library(ggplot2) }))
source("dev/pft-bci/finn_membership.R")
D <- "dev/pft-bci"
m   <- torch::torch_load(file.path(D, "results/ckpt_learned/epoch_4000model.pt"))
env <- fread(file.path(D, "data/env.csv")); coh <- fread(file.path(D, "data/initial_cohorts1985.csv"))
ch  <- FINN::CohortMat(obs_df = coh, sp = m$N_species); m$eval()
w   <- as.data.table(predict(m, env = env, init_cohort = ch, patches = 10L, device = "cpu")$long$site)
w[variable == "reg", value := value / 0.1]

agg <- w[variable %in% c("ba","trees","growth","mort","reg"),
         .(v = mean(value, na.rm = TRUE)), by = .(year, variable)]     # mean over sites*species
# stand ba/trees are sums over species; growth/mort are per-capita means
stand <- w[variable %in% c("ba","trees"), .(v = sum(value, na.rm=TRUE)), by=.(siteID,year,variable)][
              , .(v = mean(v)), by = .(year, variable)]
rates <- w[variable %in% c("growth","mort"),
           .(v = weighted.mean(value, w = 1, na.rm = TRUE)), by = .(year, variable)]
tab <- dcast(rbind(stand, rates), year ~ variable, value.var = "v")
cat("year-by-year (stand ba/trees; mean growth/mort rate):\n")
print(tab[, lapply(.SD, function(x) round(x, 3))], nrow = 31)

# correlation: does a BA drop follow a mortality spike?
setorder(tab, year); tab[, dBA := ba - shift(ba)]
cat("\ncor(mortality, next-year BA change):",
    round(cor(tab$mort[-nrow(tab)], tab$dBA[-1], use = "complete.obs"), 2), "\n")
cat("cor(growth, BA):", round(cor(tab$growth, tab$ba, use="complete.obs"), 2), "\n")

pl <- melt(tab[, .(year, ba, growth, mort)], id.vars = "year")
p <- ggplot(pl, aes(year, value)) + geom_line() + geom_vline(xintercept = seq(5,30,5), colour="grey80") +
  facet_wrap(~variable, ncol = 1, scales = "free_y") +
  labs(title = "learned annual: BA overshoot vs growth & mortality (grey = census years)") +
  theme_minimal(base_size = 11)
ggsave(file.path(D, "mechanism_check.png"), p, width = 7, height = 6, dpi = 130)
cat("\nfigure: mechanism_check.png\n")
