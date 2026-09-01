# ============================================================================
# Path A -- soft-membership PFTs on real BCI demographic responses (LOCAL)
# ============================================================================
# The PoC mechanism, applied to real BCI data WITHOUT the full FINN gap model.
# Each species' response of growth / mortality / regeneration to the 6 climate
# predictors is modelled as  beta_s = softmax(logits_s) %*% Theta_process,
# with ONE shared membership A = softmax(logits) [N_species x K] across all three
# processes (Rueger's PFT concept: a type has characteristic growth AND mort AND
# recruitment). Fit end-to-end by SGD, then ask: does the learned membership
# recover Rueger et al. 2020's 5 PFTs?
#
# This is the regression-surrogate bridge to Path B (which swaps this linear
# demographic model for FINN's full forward model). No cluster needed.
# ============================================================================

suppressMessages({ library(data.table); library(torch); library(ggplot2) })
set.seed(1); torch_manual_seed(1)
D <- "/Users/yannekkaber/working-directory/FINN/dev/pft-bci"
obs <- fread(file.path(D, "data/obs_species.csv"))
env <- fread(file.path(D, "data/env.csv"))
pft <- fread(file.path(D, "data/species_pft.csv"))          # species -> PFT_2axes
K   <- 5L
predictors <- c("Prec","SR_kW_m2","RH_prc","T_max","T_min","swp")

dat <- merge(obs, env, by = c("siteID","year"))
N   <- max(obs$species)

# ---- per-process observation tensors (presence-gated) ----------------------
# feature matrix X = [1, 6 predictors]; response transformed per process.
make_block <- function(d) {
  X   <- cbind(1, as.matrix(d[, ..predictors]))
  list(X = torch_tensor(X, dtype = torch_float32()),
       y = torch_tensor(d$y, dtype = torch_float32()),
       s = torch_tensor(d$species, dtype = torch_int64()))
}
gd <- make_block(dat[!is.na(growth) & trees > 0][, y := log(growth + 1)])
md <- make_block(dat[!is.na(mort) & trees > 0][, y := qlogis(pmin(pmax(mort, 1e-3), 1-1e-3))])
rd <- make_block(dat[trees > 0][, y := log(pmax(reg,0) + 1)])
cat("obs used -- growth:", gd$y$shape, " mort:", md$y$shape, " reg:", rd$y$shape, "\n")

# ---- model: shared membership A, per-process prototypes --------------------
np <- length(predictors) + 1L
membership_model <- nn_module("mm",
  initialize = function(N, K, np) {
    self$logits <- nn_parameter(torch_randn(N, K) * 0.01)
    self$Wg <- nn_parameter(torch_randn(K, np) * 0.1)
    self$Wm <- nn_parameter(torch_randn(K, np) * 0.1)
    self$Wr <- nn_parameter(torch_randn(K, np) * 0.1)
    self$log_sg <- nn_parameter(torch_zeros(1)); self$log_sm <- nn_parameter(torch_zeros(1))
    self$log_sr <- nn_parameter(torch_zeros(1))
  },
  forward = function() { A <- nnf_softmax(self$logits, dim = 2L)
    list(A = A, bg = A$matmul(self$Wg), bm = A$matmul(self$Wm), br = A$matmul(self$Wr)) }
)
# gaussian NLL of one process block given its per-species coefficient table b [N,np]
block_nll <- function(b, blk, log_s) {
  mu <- (blk$X * b[blk$s, ])$sum(dim = 2L)
  0.5 * (((blk$y - mu) / torch_exp(log_s))^2 + 2*log_s)$mean()
}
fit_model <- function(mod, steps = 4000L, lr = 0.05) {
  opt <- optim_adam(mod$parameters, lr = lr)
  for (i in seq_len(steps)) { opt$zero_grad()
    o <- mod(); L <- block_nll(o$bg, gd, mod$log_sg) +
                     block_nll(o$bm, md, mod$log_sm) +
                     block_nll(o$br, rd, mod$log_sr)
    L$backward(); opt$step() }
  invisible(NULL)
}

mod <- membership_model(N, K, np); fit_model(mod)
o <- mod(); A <- as.matrix(o$A$detach())

# ---- compare learned membership to Rueger's 5 PFTs -------------------------
keep    <- pft$species                                   # non-'other' species
cluster <- max.col(A[keep, ])                            # hard label from argmax
truth   <- pft$PFT_2axes
ari <- function(a, b) {                                  # adjusted Rand index
  tab <- table(a, b); n <- sum(tab)
  si <- sum(choose(rowSums(tab), 2)); sj <- sum(choose(colSums(tab), 2))
  sij <- sum(choose(tab, 2)); e <- si*sj/choose(n,2)
  (sij - e) / ((si+sj)/2 - e)
}
cat("\n=================== Path A: learned PFTs vs Rueger 5 PFTs ===============\n")
cat("Adjusted Rand Index (learned clusters vs Rueger PFT):", round(ari(cluster, truth), 3), "\n")
cm <- table(Rueger = truth, learned = cluster)
# greedily reorder learned columns to best align with PFT rows (for readability)
ord <- integer(0); avail <- seq_len(ncol(cm))
for (r in seq_len(nrow(cm))) { j <- avail[which.max(cm[r, avail])]; ord <- c(ord, j); avail <- setdiff(avail, j) }
ord <- c(ord, avail)
cat("\nConfusion (rows = Rueger PFT, cols = learned cluster, aligned):\n")
print(cm[, ord])
purity <- sum(apply(cm, 1, max)) / sum(cm)
cat("\nMean per-PFT purity (largest learned cluster / PFT size):", round(purity, 3), "\n")
cat("Prototype usage (avg membership mass):", paste(round(colMeans(A), 3), collapse=" "), "\n")

# ---- alikeness heatmap ordered by true PFT ---------------------------------
An <- A[keep, ] / sqrt(rowSums(A[keep, ]^2))
al <- An %*% t(An)
ordS <- order(truth)
same <- outer(truth, truth, `==`)[ordS, ordS]
ut <- upper.tri(same)
cat("Within-PFT vs between-PFT mean alikeness:",
    round(mean((al[ordS,ordS])[ut & same]),3), "vs",
    round(mean((al[ordS,ordS])[ut & !same]),3), "\n")

am <- al[ordS, ordS]
df <- data.table(i = rep(seq_along(keep), each=length(keep)),
                 j = rep(seq_along(keep), times=length(keep)),
                 alike = as.vector(am))
df$pft_i <- truth[ordS][df$i]
p <- ggplot(df, aes(i, j, fill = alike)) + geom_raster() +
  scale_fill_viridis_c(limits = c(0,1)) +
  labs(title = "BCI: learned species alikeness vs Rueger PFTs",
       subtitle = paste0("ordered by Rueger PFT (K=", K, "); blocks = agreement"),
       x = "species (ordered by Rueger PFT)", y = "species", fill = "alikeness") +
  theme_minimal(base_size = 11)
ggsave(file.path(D, "pathA_alikeness.png"), p, width = 6.6, height = 5.8, dpi = 130)

# membership heatmap: species (by PFT) x prototype
dm <- data.table(sp = rep(seq_along(keep), K), k = rep(seq_len(K), each=length(keep)),
                 w = as.vector(A[keep, ]), pft = rep(truth, K))
dm <- dm[order(pft)]; dm$sp_ord <- match(dm$sp, unique(dm$sp))
p2 <- ggplot(dm, aes(factor(k), sp_ord, fill = w)) + geom_raster() +
  scale_fill_viridis_c() +
  labs(title = "BCI membership A (species x prototype)",
       subtitle = "species ordered by Rueger PFT", x = "learned prototype", y = "species") +
  theme_minimal(base_size = 11)
ggsave(file.path(D, "pathA_membership.png"), p2, width = 4.8, height = 6, dpi = 130)
cat("\nFigures: pathA_alikeness.png, pathA_membership.png\n")
