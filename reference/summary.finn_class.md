# Summarise a fitted FINN model

Collects what is needed to judge a fitted model, in this order:

1.  **Skill**: per response, R2 (`1 - SSE/SST`) and the correlation
    (Pearson r) of observed and simulated first, then RMSE, the bias,
    the ratio of the mean prediction to the mean observation (is the
    magnitude right?; for mortality weighted by the trees at risk, where
    the observations carry `n_at_risk` and `n_died`) and n; for
    mortality also the tree-level AUC and the ratio of predicted to
    observed deaths. On the training data and, when
    [`evaluate()`](https://finnverse.github.io/FINN/reference/evaluate.md)
    or the argument `holdout` of
    [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) stored
    one, on the holdout. The observations are compared with the
    simulation aggregated the way
    [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) scores
    it (see
    [`comparePlot()`](https://finnverse.github.io/FINN/reference/comparePlot.md)).
    Without a stored evaluation,
    [`summary()`](https://rdrr.io/r/base/summary.html) runs
    [`evaluate()`](https://finnverse.github.io/FINN/reference/evaluate.md)
    on the training data first (one simulation) and stores it on the
    model.

2.  **Convergence**: by how many percent each loss term changed over the
    last `epochs_back` epochs (the mean of the last `epochs_back` epochs
    against the mean of the `epochs_back` before), as a table with one
    sentence on the largest change; and how much the parameters moved
    over the last `epochs_back` epochs, as a table with one row per
    process (growth, mortality, regeneration, competition) and part (its
    species parameters, where they are fitted (`optimizeSpecies = TRUE`
    in
    [`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md)),
    its environmental response, i.e. the network or the linear climate
    effects, and for regeneration the recruit dispersion and the cap),
    naming the parameter that moved most. The change is the size of the
    step the values took as a percentage of their size,
    `100 * |values now - values then| / |values then|`, over all species
    of a parameter (on the scale it acts on, e.g. "mortality: growth")
    or all weights of a network layer; 1 % means the values moved by
    about a hundredth of their size.
    [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) keeps
    the parameters of the last 51 epochs for this, so `epochs_back` can
    be at most 50; for a model fitted before this version the parameter
    snapshots of the fit (`model$param_history`, every `checkpoints`
    epochs) are used instead; they carry no epoch, so the epochs between
    them are estimated from the length of the loss history and their
    number (exact for one call of
    [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) whose
    `epochs` is a multiple of `checkpoints`). Below the explanation, one
    sentence names the species parameter that moved most and the largest
    change overall. A fit has settled when both are small; network
    weights can keep moving while the loss stays flat.

3.  **Bare ground**, when
    [`evaluate()`](https://finnverse.github.io/FINN/reference/evaluate.md)
    stored a bare-ground run: per site, the total basal area and the
    leading species with its share, as the mean over the last tenth of
    the run.

4.  **Environmental importance and average conditional effects**, per
    process and species, from the model's conditional effects (computed
    once and cached, see
    [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md)). Where
    they cannot be computed for a model, a warning says why and the
    other sections are still returned.

[`print()`](https://rdrr.io/r/base/print.html) shows it;
`print(s, epochs_back = 10)` redoes the convergence check over another
window without recomputing anything.
[`plot.finn_class()`](https://finnverse.github.io/FINN/reference/plot.finn_class.md)
draws the fitted parameters themselves.

## Usage

``` r
# S3 method for class 'finn_class'
summary(
  object,
  env = NULL,
  init_cohort = NULL,
  importance = c("ale", "permutation"),
  env_autoscale = TRUE,
  sim_seed = 42L,
  nperm = 20L,
  scale = TRUE,
  epochs_back = 50L,
  ...,
  species_names = NULL,
  parameter_names = NULL
)

# S3 method for class 'summary.finn_class'
print(x, epochs_back = x$epochs_back, species = FALSE, ...)
```

## Arguments

- object:

  (`finn_class`)  
  fitted model.

- env, init_cohort:

  (`NULL` or data)  
  interpretation data for the importance and the conditional effects;
  `NULL` uses the cached training data (see
  [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md)).

- importance:

  (`character(1)`)  
  `"ale"` (default) = analytical ALE-variance importance (cheap, from
  the cache); `"permutation"` = re-simulating permutation importance
  (see
  [`feature_importance()`](https://finnverse.github.io/FINN/reference/feature_importance.md)),
  cached on the model under the same input key.

- env_autoscale:

  (`logical(1)`)  
  see [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md).

- sim_seed:

  (`integer(1)`)  
  seed for the simulation (see
  [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md)).

- nperm:

  (`integer(1)`)  
  replicates for `importance = "permutation"`.

- scale:

  (`logical(1)`)  
  for `importance = "ale"`, divide `Var(ALE)` by the process x species
  rate variance so the importances are dimensionless and comparable
  across processes and species (Sobol-style). Default `TRUE`.

- epochs_back:

  (`integer(1)`)  
  the window of the convergence check, in epochs. Default `50`.

- ...:

  [`summary()`](https://rdrr.io/r/base/summary.html): passed to
  [`feature_importance()`](https://finnverse.github.io/FINN/reference/feature_importance.md)
  with `importance = "permutation"` (and from there to
  [`predict()`](https://finnverse.github.io/FINN/reference/predict.finn_class.md)).
  [`print()`](https://rdrr.io/r/base/print.html): unused.

- species_names, parameter_names:

  Names of the species and of the parameters of a process function of
  one's own, as for
  [`plot.finn_class()`](https://finnverse.github.io/FINN/reference/plot.finn_class.md);
  `NULL` (default): the species names the model keeps from
  [`fit()`](https://finnverse.github.io/FINN/reference/fit.md) and the
  documented parameter names of FINN's own process functions.

- x:

  A `summary.finn_class` object.

- species:

  (`logical(1)`)  
  also print the skill per species. Default `FALSE`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html): an object of class
`summary.finn_class`, a list with `skill` (per set, overall),
`skill_species` (per set and species), `convergence` (the loss history,
the recent parameters, the parameter snapshots with their estimated
spacing, and the bounds and names of the species parameters),
`bare_ground`, `importance`, `average_conditional_effects`, `method`
(the `importance` option), `importance_label`, `species_names`,
`environment_error` (`NULL`, or why item 4 could not be computed) and
`epochs_back`. [`print()`](https://rdrr.io/r/base/print.html): `x`,
invisibly.

## See also

[`evaluate()`](https://finnverse.github.io/FINN/reference/evaluate.md),
[`comparePlot()`](https://finnverse.github.io/FINN/reference/comparePlot.md),
[`plot.finn_class()`](https://finnverse.github.io/FINN/reference/plot.finn_class.md),
[`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md),
[`averageConditionalEffects()`](https://finnverse.github.io/FINN/reference/averageConditionalEffects.md),
[`feature_importance()`](https://finnverse.github.io/FINN/reference/feature_importance.md)

## Examples

``` r
library(data.table)
# a small model fitted on two Oregon FIA sites (see ?fit)
sites <- 1:2
Nsp   <- max(fia_obs_dt$species)
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
FINN.seed(1)
m <- finn(N_species = Nsp, recruits_dbh = 12.9,
          growth_process       = createProcess(~ temp + prec, growth),
          mortality_process    = createProcess(~ temp + prec, mortality),
          regeneration_process = createProcess(~ temp + prec, regeneration_capped))
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 3L, patches = 4, patch_size = 0.06,
    plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00301, ba=8.77, trees=0.852, growth=1.58e+04, mortality=32.4, regeneration=0.677)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%

s <- summary(m)
#> Evaluating on the training data (one simulation); evaluate() can add a holdout and a bare-ground run.
s
#> FINN model summary
#> ==================
#> 
#> ## Skill (observed vs simulated, aggregated as fit() scores)
#> 
#> Training data:
#>  response     R2 correlation    RMSE    bias pred/obs  n   AUC deaths pred/obs
#>       dbh  0.512      0.7764 12.7264 4.31644     1.06 14    NA              NA
#>        ba  0.588      0.8468  0.2168 0.06657     1.36 44    NA              NA
#>     trees  0.812      0.9181  0.4827 0.00000     1.00 44    NA              NA
#>    growth -1.299      0.5841  0.0124 0.00648     1.88 14    NA              NA
#>      mort -1.180      0.0198  0.0178 0.01305     3.66 14 0.517             3.8
#>       reg -0.111      0.0191  3.7885 0.38021     1.33 44    NA              NA
#> 
#> ## Convergence
#>   loss: fewer than 50 + 1 epochs recorded
#> 
#> Parameters: how much they moved over the last 2 epochs (from the parameters of every epoch).
#>   The change is the size of the step the values took, as a percentage of the size of the values:
#>   100 x |values now - values then| / |values then|, over all species of a parameter (or all weights
#>   of a network layer). 1 % means the values moved by about a hundredth of their size.
#>       process                   part               moved most change
#>        growth environmental response       nn_growth.0.weight 3.91 %
#>     mortality environmental response    nn_mortality.0.weight 0.64 %
#>  regeneration environmental response nn_regeneration.0.weight 6.57 %
#>  regeneration     recruit dispersion   par_theta_recruits_raw 3.62 %
#>  regeneration        recruitment cap                 reg_logK 3.95 %
#> 
#> ## Environment
#> 
#> ### GROWTH
#> 
#> Analytical ALE importance (rate-normalised) (species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      temp                0.3035         0.00295         0.5032        4.7e-05
#>      prec                0.0572         0.03089         0.0923        1.0e-01
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>         0.0831             0.0013                 0                 0.15438
#>         0.0289             0.1267                 0                 0.00172
#>  Alnus rubra Abies amabilis   other
#>            0          0.421 0.05146
#>            0          0.283 0.00655
#> 
#> Average conditional effects (mean; species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      temp              -0.01516        -0.00461         0.0501      -0.000741
#>      prec              -0.00658         0.01491        -0.0215      -0.034198
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>       -0.01361           -0.00439            0.0455                 0.07978
#>        0.00802           -0.04334           -0.0223                -0.00842
#>  Alnus rubra Abies amabilis    other
#>       0.0288         0.0482 -0.01724
#>       0.0115         0.0395 -0.00615
#> 
#> ### MORTALITY
#> 
#> Analytical ALE importance (rate-normalised) (species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      prec                  1.03           0.174          0.188          0.342
#>      temp                  1.03           0.174          0.188          0.342
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>           0.23              0.172                 0                   0.286
#>           0.23              0.172                 0                   0.286
#>  Alnus rubra Abies amabilis other
#>            0        0.00117 0.834
#>            0        0.00117 0.834
#> 
#> Average conditional effects (mean; species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      prec             -0.000529        -0.00051       -0.00015       0.000275
#>      temp              0.000529         0.00051        0.00015      -0.000275
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>      -0.000475          -0.000131           0.00029                 0.00044
#>       0.000475           0.000131          -0.00029                -0.00044
#>  Alnus rubra Abies amabilis     other
#>     0.000478      -2.47e-05 -0.000551
#>    -0.000478       2.47e-05  0.000551
#> 
#> ### REGENERATION
#> 
#> Analytical ALE importance (rate-normalised) (species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      prec               0.17032          0.0127       8.56e-02         0.1392
#>      temp               0.00311          0.0912       2.18e-05         0.0457
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>          0.147           7.03e-03             0.317                  0.1860
#>          0.071           5.08e-05             0.178                  0.0196
#>  Alnus rubra Abies amabilis   other
#>       0.0288         0.0839 0.04778
#>       0.2910         0.2422 0.00564
#> 
#> Average conditional effects (mean; species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      prec               -0.0823         -0.0166       0.050466         0.0567
#>      temp                0.0111          0.0446       0.000806         0.0325
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>        -0.0641          -0.009463           -0.1103                  0.0795
#>         0.0445           0.000805           -0.0828                 -0.0259
#>  Alnus rubra Abies amabilis   other
#>        0.028        -0.0439 -0.0420
#>       -0.089        -0.0745  0.0144
# the convergence check over the last two epochs instead
print(s, epochs_back = 2)
#> FINN model summary
#> ==================
#> 
#> ## Skill (observed vs simulated, aggregated as fit() scores)
#> 
#> Training data:
#>  response     R2 correlation    RMSE    bias pred/obs  n   AUC deaths pred/obs
#>       dbh  0.512      0.7764 12.7264 4.31644     1.06 14    NA              NA
#>        ba  0.588      0.8468  0.2168 0.06657     1.36 44    NA              NA
#>     trees  0.812      0.9181  0.4827 0.00000     1.00 44    NA              NA
#>    growth -1.299      0.5841  0.0124 0.00648     1.88 14    NA              NA
#>      mort -1.180      0.0198  0.0178 0.01305     3.66 14 0.517             3.8
#>       reg -0.111      0.0191  3.7885 0.38021     1.33 44    NA              NA
#> 
#> ## Convergence
#> 
#> Loss: the largest change is dbh, 122.21 % over the last 2 epochs (mean of the last 2 epochs against the 2 before).
#>          term loss_now   change
#>           dbh    1.927 122.21 %
#>            ba    0.162 -26.21 %
#>         trees    0.464   8.19 %
#>        growth    4.752 -19.85 %
#>     mortality    1.228  -1.10 %
#>  regeneration    1.275   0.37 %
#> 
#> Parameters: how much they moved over the last 2 epochs (from the parameters of every epoch).
#>   The change is the size of the step the values took, as a percentage of the size of the values:
#>   100 x |values now - values then| / |values then|, over all species of a parameter (or all weights
#>   of a network layer). 1 % means the values moved by about a hundredth of their size.
#>       process                   part               moved most change
#>        growth environmental response       nn_growth.0.weight 3.91 %
#>     mortality environmental response    nn_mortality.0.weight 0.64 %
#>  regeneration environmental response nn_regeneration.0.weight 6.57 %
#>  regeneration     recruit dispersion   par_theta_recruits_raw 3.62 %
#>  regeneration        recruitment cap                 reg_logK 3.95 %
#> 
#> ## Environment
#> 
#> ### GROWTH
#> 
#> Analytical ALE importance (rate-normalised) (species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      temp                0.3035         0.00295         0.5032        4.7e-05
#>      prec                0.0572         0.03089         0.0923        1.0e-01
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>         0.0831             0.0013                 0                 0.15438
#>         0.0289             0.1267                 0                 0.00172
#>  Alnus rubra Abies amabilis   other
#>            0          0.421 0.05146
#>            0          0.283 0.00655
#> 
#> Average conditional effects (mean; species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      temp              -0.01516        -0.00461         0.0501      -0.000741
#>      prec              -0.00658         0.01491        -0.0215      -0.034198
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>       -0.01361           -0.00439            0.0455                 0.07978
#>        0.00802           -0.04334           -0.0223                -0.00842
#>  Alnus rubra Abies amabilis    other
#>       0.0288         0.0482 -0.01724
#>       0.0115         0.0395 -0.00615
#> 
#> ### MORTALITY
#> 
#> Analytical ALE importance (rate-normalised) (species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      prec                  1.03           0.174          0.188          0.342
#>      temp                  1.03           0.174          0.188          0.342
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>           0.23              0.172                 0                   0.286
#>           0.23              0.172                 0                   0.286
#>  Alnus rubra Abies amabilis other
#>            0        0.00117 0.834
#>            0        0.00117 0.834
#> 
#> Average conditional effects (mean; species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      prec             -0.000529        -0.00051       -0.00015       0.000275
#>      temp              0.000529         0.00051        0.00015      -0.000275
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>      -0.000475          -0.000131           0.00029                 0.00044
#>       0.000475           0.000131          -0.00029                -0.00044
#>  Alnus rubra Abies amabilis     other
#>     0.000478      -2.47e-05 -0.000551
#>    -0.000478       2.47e-05  0.000551
#> 
#> ### REGENERATION
#> 
#> Analytical ALE importance (rate-normalised) (species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      prec               0.17032          0.0127       8.56e-02         0.1392
#>      temp               0.00311          0.0912       2.18e-05         0.0457
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>          0.147           7.03e-03             0.317                  0.1860
#>          0.071           5.08e-05             0.178                  0.0196
#>  Alnus rubra Abies amabilis   other
#>       0.0288         0.0839 0.04778
#>       0.2910         0.2422 0.00564
#> 
#> Average conditional effects (mean; species in columns):
#>  variable Pseudotsuga menziesii Pinus ponderosa Pinus contorta Abies concolor
#>      prec               -0.0823         -0.0166       0.050466         0.0567
#>      temp                0.0111          0.0446       0.000806         0.0325
#>  Abies grandis Tsuga heterophylla Tsuga mertensiana Lithocarpus densiflorus
#>        -0.0641          -0.009463           -0.1103                  0.0795
#>         0.0445           0.000805           -0.0828                 -0.0259
#>  Alnus rubra Abies amabilis   other
#>        0.028        -0.0439 -0.0420
#>       -0.089        -0.0745  0.0144
# the tables are also returned
s$skill$train
#>    variable     n         r2        rmse          r        bias    ratio
#>      <fctr> <int>      <num>       <num>      <num>       <num>    <num>
#> 1:      dbh    14  0.5123935 12.72637644 0.77639667 4.316443225 1.061483
#> 2:       ba    44  0.5876016  0.21683171 0.84684849 0.066568246 1.362183
#> 3:    trees    44  0.8116517  0.48265365 0.91805722 0.000000000 1.000000
#> 4:   growth    14 -1.2989358  0.01235288 0.58410426 0.006476158 1.878507
#> 5:     mort    14 -1.1797813  0.01779307 0.01979373 0.013052526 3.659553
#> 6:      reg    44 -0.1114828  3.78851236 0.01911046 0.380210490 1.334585
#>          auc deaths_ratio
#>        <num>        <num>
#> 1:        NA           NA
#> 2:        NA           NA
#> 3:        NA           NA
#> 4:        NA           NA
#> 5: 0.5170455     3.797192
#> 6:        NA           NA
```
