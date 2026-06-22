# Prototype: ALE comparison (process vs hybrid) for dbh, light, temp, prec, for a
# few species, with MC-dropout CI. States (dbh/light/trees + env) are taken from a
# simulation via record_raws (the FINNetAl approach). 4x4 grid: rows = variable.
suppressMessages({library(FINN); library(torch); library(data.table); library(ggplot2)})
ext <- function(f) system.file("extdata", f, package = "FINN")
env_dt <- fread(ext("fia_env_dt.csv")); sp_dt <- fread(ext("fia_species_dt.csv"))
init   <- fread(ext("fia_init_trees.csv"))
ic <- makeInitCohorts(init, Nspecies = max(fread(ext("fia_obs_dt.csv"))$species))
mp <- torch_load(ext("fia_process_finn.pt")); mh <- torch_load(ext("fia_hybrid_finn.pt"))
env_vars <- as.data.table(mp$env_scaling)$variable

# ---- realized cohort states from a simulation -------------------------------
get_states <- function(model) {
  model$eval(); model$record_raws <- TRUE
  invisible(model$simulate(env = env_dt, init_cohort = ic, patches = 4,
                           patch_size = 0.06, device = "cpu"))
  rg <- model$raw_g; model$record_raws <- FALSE; model$raw_g <- NULL
  env_cols <- c("siteID", env_vars)
  env_site <- unique(env_dt, by = "siteID")[, ..env_cols]   # one env row per site
  st <- rbindlist(lapply(seq_along(rg), function(t) {
    a <- rg[[t]]; d <- dim(a)
    g <- as.data.table(expand.grid(site = 1:d[1], patch = 1:d[2], cohort = 1:d[3]))
    g[, `:=`(dbh = as.vector(a[,,,1]), light = as.vector(a[,,,2]),
             trees = as.vector(a[,,,3]), species = as.vector(a[,,,4]))]
    merge(g, env_site, by.x = "site", by.y = "siteID")
  }))
  st[trees > 0.5 & dbh > 0.5]
}

# ---- predict growth for a table of states (env scaled internally) -----------
make_predf <- function(model) {
  scaling <- as.data.table(model$env_scaling); is_h <- inherits(model$process_growth, "hybrid")
  function(st, sp) {
    n <- nrow(st)
    Z <- do.call(cbind, lapply(env_vars, function(v){ s <- scaling[variable==v]; (st[[v]]-s$center)/s$scale }))
    mk <- function(v, ty) torch_tensor(array(v, c(n,1,1)), dtype = ty)
    pred <- torch_tensor(cbind(1, Z), dtype = torch_float32()); species <- mk(sp, torch_long())
    if (!is_h) { pred <- model$nn_growth(pred); pred <- FINN::index_species(pred, species) }
    as.numeric(as_array(model$growth_func(dbh = mk(st$dbh, torch_float32()),
      trees = mk(st$trees, torch_float32()), light = mk(st$light, torch_float32()),
      species = species, pred = pred, parGrowth = model$par_growth)))
  }
}

# ---- ALE of one feature for one species (all in raw units) ------------------
ale_feature <- function(predf, st_sp, sp, target, breaks) {
  K <- length(breaks) - 1L
  bin <- cut(st_sp[[target]], breaks, include.lowest = TRUE, labels = FALSE)
  local <- numeric(K)
  for (k in seq_len(K)) {
    idx <- which(bin == k); if (!length(idx)) next
    lo <- copy(st_sp[idx]); hi <- copy(st_sp[idx])
    set(lo, j = target, value = breaks[k]); set(hi, j = target, value = breaks[k+1])
    local[k] <- mean(predf(hi, sp) - predf(lo, sp))
  }
  ale <- cumsum(local) - mean(cumsum(local)) + mean(predf(st_sp, sp))
  data.table(x = (breaks[-1] + breaks[-(K+1)]) / 2, growth = pmax(ale, 0))
}

VARS <- c("dbh", "light", "temp", "prec"); SPP <- 1:4
st_p <- get_states(mp); st_h <- get_states(mh)

one_model <- function(model, st, mc, reps) {
  predf <- make_predf(model)
  rbindlist(lapply(VARS, function(v) rbindlist(lapply(SPP, function(sp) {
    st_sp <- st[species == sp]
    if (nrow(st_sp) < 30) return(NULL)
    # trim sparse top 5% of this driver, then bin from species' own range
    st_sp  <- st_sp[get(v) <= quantile(st_sp[[v]], 0.95)]
    breaks <- unique(quantile(st_sp[[v]], probs = seq(0, 1, length.out = 9)))
    if (length(breaks) < 3) return(NULL)
    curves <- rbindlist(lapply(seq_len(reps), function(j) {
      if (mc) model$train() else model$eval()
      d <- if (mc) st_sp[sample(.N, replace = TRUE)] else st_sp
      ale_feature(predf, d, sp, v, breaks)[, rep := j]
    }))
    curves[, .(growth = mean(growth), lo = pmax(mean(growth)-sd(growth),0), hi = mean(growth)+sd(growth)),
           by = x][, `:=`(variable = v, species = sp)]
  }))))
}

t0 <- Sys.time()
hyb <- one_model(mh, st_h, mc = TRUE,  reps = 20)[, model := "Hybrid (growth = NN)"]
proc<- one_model(mp, st_p, mc = FALSE, reps = 1 )[, model := "Process (mechanistic)"]
cat(sprintf("ALE grid computed in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))

both <- merge(rbind(hyb, proc), sp_dt, by = "species")
cols <- c("Process (mechanistic)" = "#1b9e77", "Hybrid (growth = NN)" = "#d95f02")
# one row per driver: facet_wrap over species shares the x-axis across species
# (scales = "free_y"); the variable name becomes the x-axis title below the row.
make_row <- function(v, xlab, top) {
  g <- ggplot(both[variable == v], aes(x, growth, colour = model, fill = model)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~species_name, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
    labs(x = xlab, y = "growth (cm/yr)", colour = NULL, fill = NULL) +
    theme_minimal() + theme(legend.position = "none")
  if (!top) g <- g + theme(strip.text = element_blank())
  g
}
library(patchwork)
p <- (make_row("dbh",  "dbh (cm)", TRUE) /
      make_row("light","available light", FALSE) /
      make_row("temp", expression(temperature~(degree*C)), FALSE) /
      make_row("prec", "precipitation (mm)", FALSE)) +
  plot_layout(guides = "collect") & theme(legend.position = "top")
ggsave("dev/_ale_grid.png", p, width = 9.5, height = 9, dpi = 110)
cat("saved dev/_ale_grid.png\n")
