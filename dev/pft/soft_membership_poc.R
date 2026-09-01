# ============================================================================
# Soft-membership PFTs for FINN -- proof of concept (0.3.0 exploration)
# ============================================================================
# Question being tested: does a learned soft membership let well-sampled species
# lend statistical strength to rare ones, purely through shared prototypes?
#
# Design (controlled, so ground truth is known):
#   - K_true = 3 functional types, each a prototype (a, b) of a logistic
#     mortality response  p(x) = sigmoid(a + b * x).
#   - N_species = 12, four per PFT, each species = its PFT prototype + jitter.
#   - Some species are made RARE: few observations, few trees each, and observed
#     only in a NARROW x band -> their slope b is nearly unidentifiable alone.
#   - Fit the SAME data three ways and compare parameter recovery vs truth:
#       (1) per-species   : independent (a_s, b_s)          [no pooling]
#       (2) soft-membership: theta_s = softmax(logits_s) %*% Proto  [learned PFTs]
#       (3) global         : one (a, b) for all             [full pooling]
#
# Reusable piece = `soft_membership` nn_module: this is exactly the layer that
# would drop into FINN's parameter machinery (theta_s used in place of the
# per-species par_<type> table).
# ============================================================================

suppressPackageStartupMessages({ library(torch); library(ggplot2) })

set.seed(1); torch_manual_seed(1)
outdir <- "dev/pft"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- ground truth ----------------------------------------------------------
K_true <- 3L
proto_true <- rbind(c(a = -1.5, b =  2.0),   # PFT 1: low baseline, strong response
                    c(a =  0.0, b = -1.0),   # PFT 2
                    c(a =  1.0, b =  0.5))    # PFT 3
N_sp  <- 12L
pft   <- rep(1:3, each = 4L)                  # true assignment
jit   <- 0.10
theta_true <- proto_true[pft, ] + cbind(rnorm(N_sp, 0, jit), rnorm(N_sp, 0, jit))

# rare species: one per PFT (indices 2, 6, 10) + species 12; the rest common
rare <- rep(FALSE, N_sp); rare[c(2L, 6L, 10L, 12L)] <- TRUE

# ---- generate observations -------------------------------------------------
sigmoid <- function(z) 1 / (1 + exp(-z))
make_obs <- function(s) {
  if (rare[s]) { n_obs <- 4L; n_tr <- 8L;  x <- runif(n_obs, -0.5, 0.5) }  # scarce + narrow band
  else         { n_obs <- 150L; n_tr <- 30L; x <- runif(n_obs, -2, 2) }    # plenty + full range
  p <- sigmoid(theta_true[s, 1] + theta_true[s, 2] * x)
  data.frame(species = s, x = x, n = n_tr, k = rbinom(n_obs, n_tr, p))
}
dat <- do.call(rbind, lapply(seq_len(N_sp), make_obs))

x_t   <- torch_tensor(dat$x,  dtype = torch_float32())
n_t   <- torch_tensor(dat$n,  dtype = torch_float32())
k_t   <- torch_tensor(dat$k,  dtype = torch_float32())
sp_ix <- torch_tensor(dat$species, dtype = torch_int64())   # 1-based

# binomial NLL given per-observation logit
nll <- function(logit) {
  p <- torch_sigmoid(logit)$clamp(1e-6, 1 - 1e-6)
  -(k_t * torch_log(p) + (n_t - k_t) * torch_log(1 - p))$mean()
}

# ---- reusable module: learned soft membership ------------------------------
soft_membership <- nn_module(
  "soft_membership",
  initialize = function(N_species, K, n_par, init_logits = NULL) {
    self$logits <- nn_parameter(if (is.null(init_logits))
      torch_randn(N_species, K) * 0.01 else init_logits)   # [N_species, K]
    self$proto  <- nn_parameter(torch_randn(K, n_par) * 0.1) # [K, n_par]
  },
  forward = function() {
    A     <- nnf_softmax(self$logits, dim = 2L)             # rows sum to 1
    theta <- A$matmul(self$proto)                           # [N_species, n_par]
    list(theta = theta, A = A)
  }
)

train <- function(step_fn, params, steps = 3000L, lr = 0.05) {
  opt <- optim_adam(params, lr = lr)
  for (i in seq_len(steps)) { opt$zero_grad(); loss <- step_fn(); loss$backward(); opt$step() }
  invisible(NULL)
}
gather_theta <- function(theta) theta[sp_ix, ]   # [N_obs, n_par] per-row species params

# (1) per-species -------------------------------------------------------------
theta_ps <- nn_parameter(torch_zeros(N_sp, 2L))
train(function() { th <- gather_theta(theta_ps)
  nll(th[, 1] + th[, 2] * x_t) }, list(theta_ps))
est_ps <- as.matrix(theta_ps$detach())

# (2) soft membership (K = 5, one MORE than truth: must discover it needs 3) ---
K <- 5L
mod <- soft_membership(N_sp, K, 2L)
train(function() { out <- mod()
  th <- gather_theta(out$theta)
  nll(th[, 1] + th[, 2] * x_t) }, mod$parameters)
out    <- mod(); A <- as.matrix(out$A$detach()); est_sm <- as.matrix(out$theta$detach())

# (3) global ------------------------------------------------------------------
theta_g <- nn_parameter(torch_zeros(1L, 2L))
train(function() { nll(theta_g[1, 1] + theta_g[1, 2] * x_t) }, list(theta_g))
est_g <- matrix(as.numeric(theta_g$detach()), N_sp, 2L, byrow = TRUE)

# ---- recovery error (Euclidean distance to true theta) ---------------------
err <- function(est) sqrt(rowSums((est - theta_true)^2))
tab <- data.frame(species = 1:N_sp, pft = pft, rare = rare,
                  per_species = err(est_ps), soft_membership = err(est_sm),
                  global = err(est_g))

summ <- aggregate(cbind(per_species, soft_membership, global) ~ rare, tab, mean)
cat("\n================ parameter-recovery error (lower = better) ============\n")
cat("mean ||theta_hat - theta_true||, split by whether the species is rare:\n\n")
print(within(summ, { rare <- ifelse(rare, "RARE (4 species)", "common (8 species)") }),
      row.names = FALSE, digits = 3)

# split by parameter: intercept a is identifiable even when rare; slope b is the
# one a narrow-band rare species cannot pin down alone -> where pooling must help.
aerr <- function(est) abs(est[, 1] - theta_true[, 1])
berr <- function(est) abs(est[, 2] - theta_true[, 2])
par_tab <- data.frame(rare = rare,
  a_per = aerr(est_ps), a_soft = aerr(est_sm),
  b_per = berr(est_ps), b_soft = berr(est_sm))
psum <- aggregate(cbind(a_per, a_soft, b_per, b_soft) ~ rare, par_tab, mean)
cat("\nBy parameter (mean |error|): intercept a vs slope b\n")
print(within(psum, { rare <- ifelse(rare, "RARE", "common") }),
      row.names = FALSE, digits = 3)

cat("\nPer rare species (this is where borrowing strength shows up):\n")
print(tab[tab$rare, c("species","pft","per_species","soft_membership","global")],
      row.names = FALSE, digits = 3)

# effective # prototypes actually used (avg membership mass per proto)
use <- colMeans(A)
cat("\nPrototype usage (avg membership mass; K was set to", K,
    "but truth is 3):\n"); print(round(use, 3))

# ---- alikeness: cosine similarity of learned membership rows ----------------
An   <- A / sqrt(rowSums(A^2))
alike <- An %*% t(An)                         # [N_sp, N_sp] in [0,1]
cat("\nWithin-PFT vs between-PFT mean alikeness (from learned A):\n")
ut <- upper.tri(alike)
same <- outer(pft, pft, `==`)
cat("  within-PFT :", round(mean(alike[ut & same]), 3),
    "\n  between-PFT:", round(mean(alike[ut & !same]), 3), "\n")

# who did each rare species end up most aligned with? (identifiability of the
# membership itself -- a rare species CAN latch onto the wrong prototype)
cat("\nRare-species membership check (best-matching common species by alikeness):\n")
for (s in which(rare)) {
  cand <- setdiff(which(!rare), s)
  best <- cand[which.max(alike[s, cand])]
  cat(sprintf("  sp%2d (true PFT %d) most alike -> sp%2d (PFT %d), alike=%.2f  %s\n",
              s, pft[s], best, pft[best], alike[s, best],
              ifelse(pft[best] == pft[s], "correct", "*** WRONG PFT ***")))
}

# ---- figures ---------------------------------------------------------------
# recovery
long <- reshape(tab[, c("species","rare","per_species","soft_membership","global")],
                varying = c("per_species","soft_membership","global"),
                v.names = "error", timevar = "model",
                times = c("per-species","soft-membership","global"),
                direction = "long")
long$grp <- ifelse(long$rare, "rare", "common")
agg <- aggregate(error ~ model + grp, long, mean)
agg$model <- factor(agg$model, levels = c("per-species","soft-membership","global"))
p1 <- ggplot(agg, aes(model, error, fill = grp)) +
  geom_col(position = position_dodge()) +
  labs(title = "Parameter recovery: soft membership rescues rare species",
       subtitle = "mean ||theta_hat - theta_true||  (lower is better)",
       x = NULL, y = "recovery error", fill = "species") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "recovery_error.png"), p1, width = 7, height = 4.2, dpi = 130)

# alikeness heatmap (ordered by true PFT)
ord <- order(pft)
am  <- alike[ord, ord]
df  <- expand.grid(i = seq_len(N_sp), j = seq_len(N_sp))
df$alike <- as.vector(am)
df$si <- factor(paste0("sp", ord[df$i], "(PFT", pft[ord][df$i], ")"),
                levels = paste0("sp", ord, "(PFT", pft[ord], ")"))
df$sj <- factor(paste0("sp", ord[df$j], "(PFT", pft[ord][df$j], ")"),
                levels = rev(paste0("sp", ord, "(PFT", pft[ord], ")")))
p2 <- ggplot(df, aes(si, sj, fill = alike)) + geom_tile() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(title = "Learned species alikeness (cosine of membership rows)",
       subtitle = "ordered by TRUE PFT -> block structure = emergent PFTs",
       x = NULL, y = NULL, fill = "alikeness") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path(outdir, "alikeness.png"), p2, width = 6.8, height = 6, dpi = 130)

cat("\nFigures written to", file.path(outdir, c("recovery_error.png","alikeness.png")), "\n")
