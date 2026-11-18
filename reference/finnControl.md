# Settings of the optimizer used by fit()

Collects the settings of the optimization behind
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) in one
object, passed as `fit(control = finnControl(...))`. The defaults are
robust in most cases: a fit needs no `control` argument at all. The
learning rate `lr`, the number of `epochs` and the `batchsize` are
changed much more often, so they stay arguments of
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) itself.

## Usage

``` r
finnControl(
  optimizer = torch::optim_ignite_adam,
  lr_scheduler = "none",
  lr_scheduler_params = list(),
  clip_norm = 2,
  ...
)

# S3 method for class 'finnControl'
print(x, ...)
```

## Arguments

- optimizer:

  (`function`)  
  A `torch` optimizer generator. Default
  [torch::optim_ignite_adam](https://torch.mlverse.org/docs/reference/optim_ignite_adam.html).

- lr_scheduler:

  (`character(1)` or `function`)  
  Learning-rate schedule: `"none"` (default), `"step"`, `"exponential"`,
  `"cosine"`, `"plateau"`, or a `function(optimizer)` returning a
  `torch` scheduler. See Details.

- lr_scheduler_params:

  (`list`)  
  Settings of the chosen schedule. See Details.

- clip_norm:

  (`numeric(1)` or named `list`)  
  Gradient-norm budget per parameter group. Default `2`.

- ...:

  Further arguments for `optimizer`, e.g. `weight_decay = 1e-4`.

- x:

  (`finnControl`)  
  The object to print.

## Value

A list of class `finnControl`, to be passed to `fit(control = )`.

## Details

[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) estimates
the parameters by gradient descent: in every epoch it simulates the
sites, compares the simulation with the observations (see
[finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md)),
and moves every parameter a small step in the direction that makes the
comparison better. `finnControl()` sets how those steps are taken.

- `optimizer` is the rule that turns a gradient into a step. The
  default, Adam
  ([torch::optim_ignite_adam](https://torch.mlverse.org/docs/reference/optim_ignite_adam.html)),
  adapts the step size per parameter and works well without tuning. Any
  `torch` optimizer generator can be passed; extra arguments for it,
  such as `weight_decay`, go into `...`.

- `lr_scheduler` lowers the learning rate during the fit. The default
  `"none"` keeps `lr` constant. `"step"` multiplies the learning rate by
  `gamma` every `step_size` epochs, `"exponential"` multiplies it by
  `gamma` every epoch, `"cosine"` lets it fall smoothly to `eta_min`
  over `T_max` epochs, and `"plateau"` lowers it by `factor` after
  `patience` epochs without an improvement of the total loss. A schedule
  helps when the losses keep oscillating late in a fit. Advanced: a
  `function(optimizer)` that returns a `torch` learning-rate scheduler
  (an object with a `$step()` method).

- `lr_scheduler_params` tunes the chosen schedule, e.g.
  `list(step_size = 50, gamma = 0.5)` for `"step"`, `list(T_max = 200)`
  for `"cosine"`, `list(factor = 0.5, patience = 20)` for `"plateau"`.
  Entries that are not given take defaults scaled to `epochs` and `lr`.

- `clip_norm` limits how large a single step can be. Gradients are
  clipped to this norm, separately for three groups of parameters (the
  species parameters of the mechanistic processes, the
  environmental-response networks, and the parameters of the
  likelihoods), so that a group with large gradients cannot swamp the
  others. A single number applies the same budget to every group; a
  named list, e.g. `list(mechanistic = 2, nn = 1, loss = 5)`, sets
  groups individually.

Before 0.2.0 these were arguments of
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) itself. The
old form is still accepted with a warning, e.g.
`fit(lr_scheduler = "step")` becomes
`fit(control = finnControl(lr_scheduler = "step"))`.

## See also

[`fit()`](https://finnverse.github.io/FINN/reference/fit.md),
[finn_loss](https://finnverse.github.io/FINN/reference/finn_loss.md)

## Examples

``` r
# the defaults, which most fits can keep
finnControl()
#> FINN optimizer settings (fit(control = ))
#>   optimizer:            Adam (torch::optim_ignite_adam) 
#>   lr_scheduler:         none 
#>   clip_norm:            2 

# halve the learning rate every 50 epochs and allow larger steps
ctrl <- finnControl(lr_scheduler = "step",
                    lr_scheduler_params = list(step_size = 50, gamma = 0.5),
                    clip_norm = 5)

library(data.table)
sites <- 1:3
Nsp   <- max(fia_obs_dt$species)
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
m <- finn(N_species = Nsp, recruits_dbh = 12.9,
          growth_process       = createProcess(~ temp + prec, growth),
          mortality_process    = createProcess(~ temp + prec, mortality),
          regeneration_process = createProcess(~ temp + prec, regeneration_capped))
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06,
    control = ctrl, plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.768)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%
```
