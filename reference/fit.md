# Fit FINN

Fit (calibrate) a FINN model end-to-end by gradient-descent
optimization. All mechanistic process parameters, environmental-response
coefficients and, if present, the weights of any hybrid (DNN-based)
processes are estimated jointly in a single optimization, rather than
calibrating each process in isolation.

## Usage

``` r
fit(
  model,
  data = NULL,
  env,
  disturbance = NULL,
  patches = 100L,
  patch_size = 0.1,
  init_cohort = NULL,
  epochs = 20L,
  lr = 0.01,
  lr_scheduler = "none",
  lr_scheduler_params = list(),
  loss = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality =
    "binomial", regeneration = "nbinom"),
  weights = "auto",
  optimizer = optim_ignite_adam,
  batchsize = NULL,
  device = c("cpu", "gpu"),
  update_step = 1L,
  start_time = 1L,
  plot_progress = TRUE,
  folder = NULL,
  checkpoints = 100L,
  shuffle = TRUE,
  record_gradients = FALSE,
  env_autoscale = TRUE,
  clip_norm = 2,
  ...
)
```

## Arguments

- model:

  (`finn_class`)  
  Object of class `finn_class` created by
  [finn](https://finnverse.github.io/FINN/reference/finn.md).

- data:

  (`data.table|data.frame`)  
  Data about demographic rates and stand variables must be passed as
  `data.table` or `data.frame`.

  Optional columns that change how the responses are scored:

  - `period_length` — how many simulated timesteps each observation's
    inventory interval spans. The rate responses (`growth`, `mortality`)
    are then compared against the model's **mean** over that interval,
    and `regeneration` against its **sum** (recruits accumulate). Absent
    or `NA` means each observation is compared against a single
    timestep. NOTE: this must currently be the SAME for every site — a
    per-site interval is not supported yet, so inventories with varying
    remeasurement gaps need filtering to a constant interval first.

  - `n_at_risk` — trees behind each observed `mortality` proportion;
    used to weight the binomial term (see `weights`).

  - `growth_n` — trees behind each observed `growth` mean; used to
    weight the squared-error term, since the variance of a mean is
    \\\sigma^2/n\\.

- env:

  (`data.table|data.frame`)  
  Data with environmental covariates must be passed as `data.table` or
  `data.frame`.

- disturbance:

  (`data.table|data.frame`)  
  Data with disturbance rates must be passed as `data.table` or
  `data.frame`.

- patches:

  (`integer(1)`)  
  Number of patches.

- patch_size:

  (`numeric(1)`)  
  Patch size.

- init_cohort:

  (`CohortMat`)  
  Initial cohort matrix of class `CohortMat`, created by
  [CohortMat](https://finnverse.github.io/FINN/reference/CohortMat.md)

- epochs:

  (`integer(1)`)  
  Number of iteration steps.

- lr:

  (`numeric(1)`)  
  Learning rate of the optimizer.

- lr_scheduler:

  (`character(1)|function`)  
  Learning-rate schedule applied on top of the constant `lr` above.
  `"none"` (default) keeps `lr` constant, reproducing the previous
  behavior exactly. Built-ins: `"step"` (decay by `gamma` every
  `step_size` epochs), `"exponential"` (multiply by `gamma` every
  epoch), `"cosine"` (cosine-anneal to `eta_min` over `T_max` epochs),
  `"plateau"` (reduce by `factor` after `patience` epochs without
  improvement in the total loss). Advanced: pass a `function(optimizer)`
  returning a `torch` lr scheduler object (must implement `$step()`).

- lr_scheduler_params:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  Tuning overrides for the scheduler chosen via `lr_scheduler`, e.g.
  `list(step_size = 50, gamma = 0.5)` for `"step"`, `list(T_max = 200)`
  for `"cosine"`, `list(factor = 0.5, patience = 20)` for `"plateau"`.
  Unset entries fall back to defaults scaled off `epochs`/`lr`.

- loss:

  (`character(6)`)  
  Named vector of the different losses. Names should be `dbh`, `ba`,
  `trees`, `growth`, `mortality`, and `regeneration`. Supported losses
  are `mse`, `poisson`, `nbinom`, `gaussian`, and `binomial`. `binomial`
  is a Bernoulli/binomial negative log-likelihood intended for
  `mortality`; it expects both the prediction and the observation to be
  proportions in \[0, 1\].

- weights:

  (`"auto"` or `numeric(6)`)  
  Weights of the six losses, in the order `dbh`, `ba`, `trees`,
  `growth`, `mortality`, `regeneration`.

  **Weights must account for the raw scale of each loss, not just its
  importance.** The six terms are summed, and their raw magnitudes
  differ by orders of magnitude: on the bundled FIA data `dbh` is an MSE
  in cm^2 (~430 at the start of training) while `growth` is an MSE on a
  ratio (~0.03) — a factor of ~1e4. Weights chosen by importance alone
  are therefore dominated by whichever response happens to have the
  largest units, and the rest receive too little gradient to learn from.
  With equal weights, `dbh` takes ~87% of the FIA objective and `growth`
  ~1%; `growth` then never improves.

  `"auto"` (the default) fixes this without needing to know the units.
  Each loss is divided by its **intercept-only baseline** — the loss you
  would get by predicting the single best constant (the null model).
  Every term then measures the same thing on the same scale: the
  fraction of its own null deviance. For a squared-error term the
  baseline is exactly \\var(y)\\, so this reduces to scaling by
  \\1/sd(y)^2\\; but it generalises to the Poisson, negative-binomial
  and binomial terms, where a standard deviation is not the right scale.
  It is the same idea as `env_autoscale = TRUE`, applied to the
  responses rather than the predictors.

  A useful side effect: the per-response values in `m$history` become
  directly interpretable as **"how much better than the mean"** — `1` is
  no better than the intercept, below `1` is better, above `1` is worse.
  The baselines are stored in `m$loss_baseline`, the resulting weights
  in `m$loss_weights`, and both are reported once when fitting starts.

  Passing a `numeric(6)` uses those weights as-is and disables the
  scaling. On the FIA example, balancing the objective is worth ~+0.24
  Spearman on held-out `growth` — far more than tuning `lr` or `epochs`.
  See the "Fitting FINN to forest inventory data" vignette.

- optimizer:

  (`torch_optimizer_generator`)  
  Optimizer from the `torch` package.

- batchsize:

  (`integer(1)`)  
  Batch size, model will be trained in random batch sizes of the data to
  preserve memory and improve convergence.

- device:

  (`character(1)`)  
  Should the model be fitted on the CPU or the GPU (Graphic card).
  Support is only for NVIDIA GPUs available.

- update_step:

  (`integer(1)`)  
  Number of steps for which the gradient should be calculated. Automatic
  differentation becomes slow for larger update steps and the risk of
  vanishing gradients increases.

- start_time:

  (`integer(1)`)  
  Starting from which year should the model be fitted. Can be used to
  use on burn-in.

- plot_progress:

  (`logical(1)`)  
  Plot fitting progress (losses) or not.

- folder:

  (`character(1)`)  
  Path to folder for saving checkpoint models. If `NULL`, no models will
  be saved during the training.

- checkpoints:

  (`integer(1)`)  
  Interval size in epochs for saving checkpoint models.

- shuffle:

  (`logical(1)`)  
  Shuffle data or not.

- record_gradients:

  (`logical(1)`)  
  Record the gradients of all parameters or not. Can get large for many
  epochs.

- env_autoscale:

  (`logical(1)`)  
  If `TRUE`, FINN z-standardizes the environmental predictors in `env`
  internally: the per-variable mean and standard deviation are learned
  from the training `env` and stored on the model, then re-applied
  automatically at every
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)/[`simulate()`](https://rdrr.io/r/stats/simulate.html)
  call. This lets you pass raw (untransformed) `env` for both
  calibration and prediction; FINN guarantees an identical
  transformation at both stages. Recommended (and the default) for
  numerical stability when predictors are on different scales. Set
  `FALSE` to use `env` exactly as supplied (e.g. if you have already
  standardized it yourself). The learned constants are available as
  `model$env_scaling`.

- clip_norm:

  (`numeric(1)|list()`)  
  Gradient-norm budget passed to
  [`torch::nn_utils_clip_grad_norm_()`](https://torch.mlverse.org/docs/reference/nn_utils_clip_grad_norm_.html),
  applied separately to each of three parameter groups (mechanistic
  per-species rates, env-effect networks, loss-distribution nuisance
  parameters) rather than once globally across all parameters. A single
  number (default `2.0`) applies the same budget to every group; a named
  list/vector keyed by `"mechanistic"`/`"nn"`/`"loss"` overrides
  individual groups, e.g. `clip_norm = list(loss = 5, nn = 1)`.

- ...:

  Additional arguments passed to `optimizer`.

## Value

The fitted `model`, invisibly. `fit()` trains the model in place, so the
returned object is the one passed in, now carrying the training results
(e.g. `$history`, `$loss_weights`, `$loss_baseline`).

## Details

Calibration relies on the fact that FINN is implemented in `torch` for R
and is therefore fully differentiable: at each simulated time step the
predicted stand variables and demographic rates are compared to the
observed data through a joint loss function, and gradients of this loss
with respect to every model parameter are obtained via automatic
differentiation (backpropagation) and used to update the parameters with
a `torch` optimizer (Adam,
[torch::optim_ignite_adam](https://torch.mlverse.org/docs/reference/optim_ignite_adam.html),
by default).

The joint loss is the sum of per-variable losses (akin to negative
log-likelihoods), one for each of `dbh`, `ba` (basal area), `trees`
(number of trees), `growth`, `mortality` and `regeneration`. Following
Pichler & Käber (2026), reasonable choices are mean squared error
(`"mse"`, equivalent to a Gaussian likelihood) for `dbh` and `ba`,
Poisson likelihood for `trees`, negative binomial (`"nbinom"`) for
`regeneration`, and `"mse"` for `growth` as a continuous rate.
`mortality` is an observed *proportion* of trees that died and the model
predicts it through a sigmoid, so it defaults to `"binomial"` — a
Bernoulli/binomial likelihood (binary cross-entropy, which admits
fractional targets). This respects the \[0, 1\] support and the
mean-variance link of a proportion, both of which `"mse"` ignores (the
model also supports `"gaussian"` and `"poisson"` as alternatives via the
`loss` argument; see Appendix B of the paper for details). Each loss can
be weighted individually via `weights`, and missing values in the
observed data are masked out of the corresponding loss term. The model
is trained for `epochs` iterations over the (optionally batched and
shuffled) data using `optimizer` with learning rate `lr`.

Backpropagating gradients through a long simulated time series is prone
to vanishing gradients and is computationally expensive. FINN therefore
uses truncated backpropagation through time: the computational graph is
detached every `update_step` simulated years, gradients are accumulated
and the parameters are updated, before the simulation continues. Smaller
`update_step` values are faster and avoid vanishing gradients but
provide a more short-sighted learning signal; larger values let the loss
integrate over a longer trajectory at increased computational cost.
`start_time` allows discarding an initial burn-in period of the
simulation from the loss, and `checkpoints`/`folder` allow periodically
saving the model state during training.

`growth` is compared as a **relative** rate (`dbh/dbh_before - 1`),
which is the model's native parameter (`dbh_new = dbh * (1 + g)`).
Training against the absolute diameter increment (`dbh * g`) instead was
tested over three seeds and was worse on every one — even when scored on
the absolute scale — and about four times more seed-variable, because it
couples the growth parameter to whichever trees happen to be present.
