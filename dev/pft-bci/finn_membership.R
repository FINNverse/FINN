# ============================================================================
# finn_membership -- soft-membership PFT reparameterization of FINN (0.3.0 PoC)
# ============================================================================
# A dev-only subclass of FINN's internal finn_class. It leaves the entire FINN
# forward model, fit(), predict() and data pipeline untouched, and only changes
# HOW the per-species mechanistic parameters are produced:
#
#     par_<type>_s  =  forward( A[s,] %*% Proto_<type> , lower, upper )
#     A = softmax(mm_logits)              [N_species x K]   (ONE shared membership)
#     Proto_<type>                        [K x n_par]        per-process prototypes
#
# So a species' growth / mortality / regeneration / competition parameters are a
# convex mix of K shared prototypes ("PFTs"). K = N_species recovers free
# per-species; small K forces species to borrow strength through shared
# prototypes. Env-effect parameters are left per-species (free) -- this isolates
# the effect of pooling the *demographic* parameters, which is what defines a PFT.
#
# Because par_<type> is an ACTIVE BINDING in finn_class and par_<type>_unconstrained
# is referenced nowhere else, overriding the four bindings fully controls the
# parameterization; the leaf par_<type>_unconstrained tensors are simply never
# created, so fit()'s optimizer (self$parameters) trains mm_logits + Proto_* +
# the untouched env/loss/theta parameters.
# ============================================================================

library(torch)
library(FINN)

.FINN_finn_class <- getFromNamespace("finn_class", "FINN")
.FINN_forward    <- getFromNamespace("forward",    "FINN")   # sigmoid(x)*(hi-lo)+lo
.FINN_backward   <- getFromNamespace("backward",   "FINN")   # qlogis((v-lo)/(hi-lo))

# getter shared by the four overridden bindings
.mm_get_par <- function(self, type) {
  upper <- self[[paste0("par_", type, "_upper")]]
  lower <- self[[paste0("par_", type, "_lower")]]
  proto <- self[[paste0("mm_Proto_", type)]]
  if (is.null(proto)) {                        # membership not installed yet
    return(torch::torch_tensor(self[[paste0(".mm_init_", type)]], dtype = torch::torch_float32()))
  }
  A   <- torch::nnf_softmax(self$mm_logits, dim = 2L)
  unc <- A$matmul(proto)                       # [N_species, n_par] in logit space
  .FINN_forward(unc, upper, lower)
}
# setter: capture the init that finn_class's setup passes in, for prototype init
.mm_set_par <- function(self, type, value) {
  self[[paste0(".mm_init_", type)]] <- as.matrix(value)
}

finn_membership <- torch::nn_module(
  classname = "finn_membership",
  inherit   = .FINN_finn_class,

  # mm_K         : number of prototypes (ignored if mm_membership is given).
  # mm_membership: NULL -> membership A is LEARNED (soft). An integer vector of
  #                length N_species -> A is FROZEN to that hard assignment
  #                (one-hot), so only the prototypes are learned. This emulates a
  #                fixed PFT scheme (e.g. Rueger) at species resolution, as a
  #                control for the learned grouping.
  initialize = function(..., mm_K = 5L, mm_membership = NULL) {
    self$mm_membership <- mm_membership
    self$mm_K <- if (is.null(mm_membership)) as.integer(mm_K) else length(unique(mm_membership))
    super$initialize(...)          # runs the full finn_class setup (captures inits)
    private$install_membership()   # then build shared membership + prototypes
  },

  active = list(
    par_growth       = function(value) if (missing(value)) .mm_get_par(self, "growth")       else .mm_set_par(self, "growth", value),
    par_mortality    = function(value) if (missing(value)) .mm_get_par(self, "mortality")    else .mm_set_par(self, "mortality", value),
    par_regeneration = function(value) if (missing(value)) .mm_get_par(self, "regeneration") else .mm_set_par(self, "regeneration", value),
    par_competition  = function(value) if (missing(value)) .mm_get_par(self, "competition")  else .mm_set_par(self, "competition", value),
    mm_A             = function() torch::nnf_softmax(self$mm_logits, dim = 2L)
  ),

  private = list(
    install_membership = function() {
      K <- self$mm_K; N <- self$N_species
      if (is.null(self$mm_membership)) {
        # learned soft membership: near-uniform init
        self$register_parameter("mm_logits", torch::nn_parameter(torch::torch_randn(N, K) * 0.01))
      } else {
        # frozen hard membership: one-hot logits as a buffer (not optimized)
        grp  <- as.integer(as.factor(self$mm_membership))     # 1..K
        oneh <- matrix(-10, N, K); oneh[cbind(seq_len(N), grp)] <- 10
        self$register_buffer("mm_logits", torch::torch_tensor(oneh, dtype = torch::torch_float32()))
      }
      for (type in c("growth", "mortality", "regeneration", "competition")) {
        init  <- self[[paste0(".mm_init_", type)]]                 # [N, n_par], bounded
        upper <- as.numeric(self[[paste0("par_", type, "_upper")]])
        lower <- as.numeric(self[[paste0("par_", type, "_lower")]])
        umean <- as.numeric(.FINN_backward(matrix(colMeans(init), 1), upper, lower))  # [n_par] logit space
        proto <- matrix(rep(umean, each = K), nrow = K) +          # [K, n_par]: all rows ~ pooled mean
                 matrix(stats::rnorm(K * length(umean), 0, 0.05), nrow = K)
        self$register_parameter(paste0("mm_Proto_", type),
                                torch::nn_parameter(torch::torch_tensor(proto, dtype = torch::torch_float32())))
      }
    }
  )
)
