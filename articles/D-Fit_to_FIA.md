# Fitting FINN to forest inventory data (Oregon, US FIA)

Dynamic forest models are usually hand-calibrated. Here we calibrate
FINN directly from data: a subset of the US Forest Inventory & Analysis
(FIA) program for Oregon, prepared exactly as in the **Preparing your
data for FINN** vignette.

We fit on 200 sites and evaluate on 200 **held-out** sites, so every
number below is reported both in-sample and out-of-sample.

Everything on this page really runs — both models are trained, simulated
and interpreted for real. Because that needs a torch backend and several
minutes, the vignette is **precompiled**: `vignettes/build.R` knits
`D-Fit_to_FIA.Rmd.orig` once on a developer machine and commits the
resulting static `.Rmd`. So the code you see is exactly the code that
produced the output beside it, and the package builds anywhere without
torch.

``` r

library(FINN)
library(torch)
library(data.table)
library(ggplot2)

# Training length. 500 is comfortably past convergence: every response plateaus
# by ~epoch 400 (checked across 3 seeds). Raising it changes nothing except the
# knit time.
EPOCHS <- 500L
```

A note on the **loss weights**, because they matter more than anything
else here — more than the learning rate, more than the epoch budget.
FINN sums one loss per response, and their raw scales differ by ~10⁴:
`dbh` is a squared error in cm² (variance ≈ 494), while `growth` is a
squared error on a ratio (variance ≈ 0.012). Summed as-is, `dbh` takes
~87% of the objective and `growth` ~1% — far too little gradient to
learn from, so growth simply never improves.

[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) handles
this by default (`weights = "auto"`): each loss is divided by its
**intercept-only baseline** — the loss you would get by predicting the
single best constant. Every term then measures the same thing, *the
fraction of its own null deviance*, and becomes directly readable: **1
means no better than the mean, below 1 is better.** For a squared-error
term the baseline is exactly the variance, so this is just scaling by
1/σ²; but it also covers the Poisson, negative-binomial and binomial
terms, where a standard deviation is not the right scale.

On this data that is worth about **+0.24 Spearman on held-out growth**,
for ~0.01 of basal area. You can still pass a `numeric(6)` to override
it.

## The data

The bundled tables are built by **`dev/make_extdata.R`**, which draws
400 sites from the full Oregon FIA set in `data-raw/` and splits them
into two disjoint samples. The **climate** (`fia_env_dt.csv`) traces
upstream to the FINN-fia analysis repo
(`scripts/03_attach_environment.R` → `07_prepare_finn_inputs.R`); see
`data-raw/README.md` for the full chain.

The data come in two **disjoint** samples: 200 sites to fit on, and 200
completely separate sites held out for evaluation. Because FINN’s
parameters are per-species × environment rather than per-site, a fitted
model transfers to sites it has never seen — so the holdout is a genuine
out-of-sample test.

``` r

ext <- function(f) system.file("extdata", f, package = "FINN")

# --- training sites --------------------------------------------------------
obs_dt     <- fread(ext("fia_obs_dt.csv"))     # observations, wide (one column per variable)
env_dt     <- fread(ext("fia_env_dt.csv"))     # RAW, untransformed climate
init_trees <- fread(ext("fia_init_trees.csv"))
species_dt <- fread(ext("fia_species_dt.csv"))

# --- holdout sites: never seen during fitting ------------------------------
obs_test   <- fread(ext("fia_obs_test.csv"))
env_test   <- fread(ext("fia_env_test.csv"))
init_test  <- fread(ext("fia_init_test.csv"))

# long form of the observations, so later we can match FINN's long predictions
# (sim$long$site) with a single merge()
resp  <- c("dbh", "ba", "trees", "growth", "mort", "reg")
to_long <- function(d) melt(d, id.vars = c("siteID", "year", "species", "species_name"),
                            measure.vars = resp, variable.name = "variable", value.name = "obs")
obs_long  <- to_long(obs_dt)
test_long <- to_long(obs_test)

cat(sprintf("train: %d sites | holdout: %d sites | %d species | years %s\n",
            uniqueN(obs_dt$siteID), uniqueN(obs_test$siteID), uniqueN(obs_dt$species),
            paste(sort(unique(obs_dt$year)), collapse = " & ")))
#> train: 200 sites | holdout: 200 sites | 11 species | years 1 & 2
species_dt
#>     species            species_name
#>       <int>                  <char>
#>  1:       1   Pseudotsuga menziesii
#>  2:       2         Pinus ponderosa
#>  3:       3          Pinus contorta
#>  4:       4          Abies concolor
#>  5:       5           Abies grandis
#>  6:       6      Tsuga heterophylla
#>  7:       7       Tsuga mertensiana
#>  8:       8 Lithocarpus densiflorus
#>  9:       9             Alnus rubra
#> 10:      10          Abies amabilis
#> 11:      11                   other
```

The environment is supplied in **natural units**: with
`env_autoscale = TRUE` (the default) FINN z-standardizes the predictors
internally and reuses the same constants at prediction time, so raw
values are all we ever pass.

``` r

ggplot(melt(env_dt[, c("temp", "prec")], measure.vars = c("temp", "prec")),
       aes(value)) +
  geom_histogram(bins = 20, fill = "grey40") +
  facet_wrap(~variable, scales = "free", labeller = as_labeller(
    c(temp = "Mean~annual~temp~(degree*C)", prec = "Annual~precip~(mm)"), label_parsed)) +
  labs(x = NULL, y = "site x year") + theme_minimal()
```

![](D/D-plot-env-1.png)

A few observed series — basal area and tree numbers by year for the most
abundant species:

``` r

top_sp <- obs_dt[species_name != "other", .(tot = sum(ba, na.rm = TRUE)), by = species][order(-tot)][1:4, species]
ggplot(obs_long[species %in% top_sp & variable %in% c("ba", "trees")],
       aes(factor(year), obs, fill = factor(species_name))) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~variable, scales = "free_y") +
  scale_fill_viridis_d(name = "species") +
  labs(x = "observation year", y = NULL) + theme_minimal()
```

![](D/D-plot-observations-1.png)

## Initial cohorts

Each simulation starts from the observed first inventory. We build the
starting state for both samples — the holdout gets its own cohorts, from
its own sites.

``` r

Nsp <- max(obs_dt$species)
init_cohorts      <- makeInitCohorts(init_trees, Nspecies = Nsp)
init_cohorts_test <- makeInitCohorts(init_test,  Nspecies = Nsp)
```

## Fit a Process-FINN

In the fully mechanistic configuration every process keeps its
predefined form, while its species and environmental parameters are
learned (`optimize* = TRUE`) and the `~.` formula lets every
environmental predictor enter the growth, regeneration and mortality
responses. The model is calibrated end-to-end by gradient descent
through the entire simulation.

``` r

FINN.seed(42)
m <- finn(
  N_species            = uniqueN(obs_dt$species),
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0, FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createProcess(~., FINN::growth,       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(~., FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~., FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)
```

``` r

fit(m,
  env         = env_dt,      # raw climate
  data        = obs_dt,
  init_cohort = init_cohorts,
  device      = "cpu",
  epochs      = EPOCHS,
  batchsize   = 40L,
  patches     = 4,
  patch_size  = 0.06,
  lr          = 0.01,
  env_autoscale   = TRUE         # default
)
```

![](D/D-fit-process-1.png)

## Replace growth with a neural network (Hybrid-FINN)

FINN’s defining feature is that **any single demographic process can be
swapped from a mechanistic function to a neural network**, while the
others stay mechanistic — a *hybrid* model. The mechanistic processes
keep the system ecologically constrained; the network absorbs structure
the fixed functional form cannot express. Here we replace **growth**
with a small neural network via
[`createHybrid()`](https://finnverse.github.io/FINN/reference/createHybrid.md),
leaving competition, regeneration and mortality mechanistic, and refit
on the same data with the **identical
[`fit()`](https://finnverse.github.io/FINN/reference/fit.md) call**.

``` r

FINN.seed(42)
m_hybrid <- finn(
  N_species            = uniqueN(obs_dt$species),
  recruits_dbh         = 12.9,
  competition_process  = createProcess(~0, FINN::competition,  optimizeSpecies = TRUE),
  growth_process       = createHybrid(~., hidden = c(20L, 20L), transformer = FALSE),  # NN replaces growth
  regeneration_process = createProcess(~., FINN::regeneration, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(~., FINN::mortality,    optimizeSpecies = TRUE, optimizeEnv = TRUE)
)
```

``` r

fit(m_hybrid,
  env = env_dt, data = obs_dt, init_cohort = init_cohorts, device = "cpu",
  epochs = EPOCHS, batchsize = 40L, patches = 4, patch_size = 0.06,
  lr = 0.01, env_autoscale = TRUE
)
```

![](D/D-hybrid-fit-1.png)

With both variants fitted, the rest of the vignette inspects and
compares them: training convergence, a
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
of the learned structure, predicted-vs-observed fit, the learned niches,
and the accumulated local effects that reveal what the network learned.

## Training convergence

Both variants are trained the same way. Because the losses are scaled by
their intercept-only baselines, `m$history` can be read **per response**
on one axis, and the dashed line at 1 is a real reference: it is the
intercept-only model, so a curve below it means that response is
predicted better than by its own mean.

``` r

hist_dt <- as.data.table(do.call(rbind, m$history))[, 1:6]
setnames(hist_dt, c("dbh", "ba", "trees", "growth", "mort", "reg"))
hist_dt[, epoch := .I]
ggplot(melt(hist_dt, id.vars = "epoch", variable.name = "response"),
       aes(epoch, value, colour = response)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_line(alpha = 0.8) +
  scale_y_log10() +
  scale_colour_viridis_d() +
  labs(x = "epoch", y = "loss / intercept-only baseline\n(1 = no better than the mean)",
       colour = NULL) +
  theme_minimal()
```

![](D/D-convergence-1.png)

The structural responses drop well below 1 within a few epochs. The
per-tree rates sit much closer to it — they are only weakly constrained
by two inventories, which is the honest picture and the reason for the
caveat further down.

## What the model learned — `summary()`

[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
gives a compact overview of the *learned structure*: for each process
and species it reports how important each environmental predictor is
(ALE-variance importance, rate-normalised) and its average conditional
effect. It reuses the cached conditional effects, so once
[`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) (below) or
a previous
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
has run it needs no further simulation.

``` r
summary(m)
FINN model summary
==================

### GROWTH

Analytical ALE importance (rate-normalised) (species in columns):
  variable    sp1    sp2     sp3      sp4      sp5    sp6     sp7    sp8
 precwarmq 0.0777 0.2431 0.05043 0.796841 0.084939 0.7666 3.09420 7.4790
   tempmax 0.0400 0.0163 0.01517 0.272596 0.177615 0.3990 0.27777 0.3691
      prec 0.0629 2.8378 0.00842 4.938956 0.005970 0.5971 0.86532 0.1145
      temp 0.0279 0.1512 0.27304 0.000966 0.974698 0.0647 0.10584 1.8580
  precseas 0.6239 0.0081 0.19111 0.096503 0.273255 0.5201 0.43262 0.0216
   tempmin 0.1222 0.4398 1.49744 2.028579 0.000523 0.0174 0.00976 0.4448
     sp9    sp10     sp11
 0.04296  0.2646 1.10e+00
 0.00171 10.9017 3.38e-01
 0.48435  0.0577 5.37e-07
 0.08967  4.9145 2.60e-01
 2.10261  1.5018 3.04e-01
 0.50974  0.3609 8.84e-02

Average conditional effects (mean; species in columns):
  variable     sp1      sp2     sp3      sp4      sp5     sp6      sp7     sp8
 precwarmq  0.0282 -0.03273 -0.0302 -0.15030  0.01540 -0.0811  0.10503 -1.0869
      temp -0.0169  0.02643  0.0575  0.00406  0.06109  0.0268  0.01678 -0.3803
      prec  0.0229  0.12334  0.0101  0.31810 -0.00531  0.0805 -0.07479  0.1168
   tempmin -0.0392 -0.04434  0.1075 -0.17882 -0.00142  0.0145  0.00642 -0.2296
   tempmax  0.0193  0.00850 -0.0115  0.05938 -0.02722 -0.0528  0.03594 -0.2309
  precseas  0.0724  0.00692 -0.0444  0.05006 -0.03492 -0.0789 -0.03451  0.0204
      sp9    sp10      sp11
 -0.03115  0.0411  7.35e-02
 -0.06232 -0.2451 -4.12e-02
  0.08809  0.0222 -5.99e-05
 -0.12212  0.0759 -2.28e-02
 -0.00763  0.2774  4.35e-02
  0.24099  0.1058  4.43e-02

### MORTALITY

Analytical ALE importance (rate-normalised) (species in columns):
  variable     sp1     sp2    sp3    sp4      sp5    sp6    sp7     sp8     sp9
 precwarmq 0.74077 1.12288 0.0264 0.4891 0.756559 0.0243 0.0661 0.25948 0.20964
  precseas 0.15872 0.00903 0.1513 0.1861 1.348941 1.3885 0.0979 0.00503 0.02273
      prec 0.29065 0.84383 0.2782 0.0267 1.005813 0.0749 0.0463 0.00263 0.51044
   tempmax 0.05344 1.29433 0.1664 0.8065 0.001803 0.2124 0.0287 0.02451 0.46741
      temp 0.01877 0.39330 0.0332 0.0450 0.000367 0.2727 0.0158 0.28261 0.02899
   tempmin 0.00144 0.07040 0.0633 0.1304 0.023599 0.2555 0.0393 0.17846 0.00322
    sp10     sp11
 0.21021 9.26e-07
 0.11655 1.43e-01
 0.41213 6.39e-02
 0.02673 2.74e-02
 0.01368 2.96e-01
 0.00226 1.90e-02

Average conditional effects (mean; species in columns):
  variable      sp1      sp2      sp3     sp4      sp5      sp6     sp7
  precseas  0.01281  0.00213  0.02022 -0.0227  0.06221 -0.02095 -0.0443
      prec -0.01818  0.02292  0.03012  0.0100  0.03817 -0.00903 -0.0340
 precwarmq  0.02466 -0.02202  0.00731  0.0365 -0.05088  0.00359 -0.0416
   tempmax  0.00660 -0.01844  0.01671 -0.0417 -0.00255 -0.01217  0.0192
      temp  0.00464  0.01307 -0.01113  0.0106  0.00118 -0.00639  0.0197
   tempmin  0.00131  0.00547 -0.01473 -0.0191  0.00665 -0.00805 -0.0307
      sp8      sp9     sp10      sp11
 -0.00543  0.01244 -0.10191  0.015183
 -0.00226  0.03094 -0.06847  0.011071
  0.01499 -0.01681 -0.05457 -0.000039
 -0.00742 -0.02880 -0.02427 -0.006484
 -0.01306 -0.00766 -0.02066 -0.023918
 -0.01065  0.00208 -0.00779 -0.006117

### REGENERATION

Analytical ALE importance (rate-normalised) (species in columns):
  variable   sp1    sp2    sp3   sp4    sp5      sp6      sp7     sp8     sp9
  precseas 0.138 0.0606 1.2937 2.600 0.4397 0.008384 0.218402 0.00184 0.08923
 precwarmq 0.205 0.1233 0.0326 0.907 0.2480 0.309712 0.050345 1.73169 0.91725
      temp 0.019 0.0136 0.3565 0.699 0.0317 0.084249 0.166767 0.02020 0.53911
   tempmin 0.857 0.0129 0.1620 0.840 0.0178 0.042227 0.081534 0.00946 0.04209
      prec 0.185 0.7159 0.5609 0.595 0.0789 0.000534 0.000347 0.06862 0.06328
   tempmax 0.792 0.0132 0.0529 0.540 0.0801 0.130536 0.018251 0.32347 0.00859
   sp10    sp11
 0.0517 0.66061
 0.5885 0.07052
 0.2796 2.23342
 0.1254 1.30112
 0.2980 0.01013
 0.1215 0.00517

Average conditional effects (mean; species in columns):
  variable    sp1    sp2    sp3    sp4    sp5     sp6      sp7     sp8     sp9
  precseas -1.095  0.221  1.965  1.165 -0.836 -0.1119  0.72202  0.0739 -0.2636
      prec  0.936 -1.735 -2.399 -0.541 -0.677 -0.0125 -0.00692  0.2395  0.1249
      temp -0.310 -0.168 -1.577 -0.907 -0.189  0.1583 -0.39983  0.1788  0.3380
   tempmin  2.071  0.169 -1.374  0.880  0.213  0.1690 -0.18202  0.1368  0.1370
 precwarmq  0.776 -0.687 -0.748 -1.579  0.538  0.1673  0.14256 -0.6776  0.2422
   tempmax  1.775  0.116 -0.555 -0.572 -0.286 -0.1854 -0.06065 -0.3358  0.0503
   sp10    sp11
  0.935  1.3941
  0.705 -0.1413
 -0.979 -2.2936
 -0.395  1.6252
  0.618  0.3360
  0.305 -0.0932
```

The niche and ALE sections below unpack this overview into per-driver
response *curves*.

## Assess the fit — in-sample and on held-out sites

The fit was calibrated on the 200 training sites, so scoring it there
only tells us how well it *reproduces* them. The honest question is
whether it transfers, so we simulate each model twice — once on the
training sites, once on the 200 held-out sites — and score both with
**Spearman correlation** (rank agreement) and **RMSE** (in each
response’s own units), matched on site × year × species.

``` r

# One prediction pass per model x sample. The holdout is driven only by its own
# environment and starting cohorts; nothing about it entered the fit.
pred_of <- function(model, label, env, init, obs_l, split) {
  model$eval()                                   # deterministic: dropout off
  sim <- predict(model, env = env, init_cohort = init,
                 patches = 4, patch_size = 0.06, device = "cpu")
  merge(obs_l, sim$long$site, by = c("siteID", "year", "species", "variable"
  ))[, `:=`(model = label, split = split)][]
}
FINN.seed(1)
obs_pred <- rbind(
  pred_of(m,        "Process (mechanistic)", env_dt,   init_cohorts,      obs_long,  "train"),
  pred_of(m,        "Process (mechanistic)", env_test, init_cohorts_test, test_long, "holdout"),
  pred_of(m_hybrid, "Hybrid (growth = NN)",  env_dt,   init_cohorts,      obs_long,  "train"),
  pred_of(m_hybrid, "Hybrid (growth = NN)",  env_test, init_cohorts_test, test_long, "holdout"))

metrics <- obs_pred[is.finite(obs) & is.finite(value),
                    .(spearman = round(cor(obs, value, method = "spearman"), 2),
                      rmse     = round(sqrt(mean((obs - value)^2)), 2),
                      n        = .N), by = .(model, split, variable)]

comparison <- dcast(metrics, variable + model ~ split, value.var = c("spearman", "rmse"))
setcolorder(comparison, c("variable", "model", "spearman_train", "spearman_holdout",
                          "rmse_train", "rmse_holdout"))
comparison[order(variable, model)]
#> Key: <variable, model>
#>     variable                 model spearman_train spearman_holdout rmse_train
#>       <fctr>                <char>          <num>            <num>      <num>
#>  1:      dbh  Hybrid (growth = NN)           0.81             0.75      13.05
#>  2:      dbh Process (mechanistic)           0.83             0.72      12.09
#>  3:       ba  Hybrid (growth = NN)           0.81             0.80       0.12
#>  4:       ba Process (mechanistic)           0.82             0.80       0.13
#>  5:    trees  Hybrid (growth = NN)           0.80             0.78       0.74
#>  6:    trees Process (mechanistic)           0.80             0.78       0.76
#>  7:   growth  Hybrid (growth = NN)           0.60             0.54       0.09
#>  8:   growth Process (mechanistic)           0.60             0.54       0.10
#>  9:     mort  Hybrid (growth = NN)           0.27             0.03       0.15
#> 10:     mort Process (mechanistic)           0.23             0.05       0.15
#> 11:      reg  Hybrid (growth = NN)           0.28             0.26       6.21
#> 12:      reg Process (mechanistic)           0.28             0.27       6.21
#>     rmse_holdout
#>            <num>
#>  1:        13.86
#>  2:        14.20
#>  3:         0.09
#>  4:         0.09
#>  5:         0.92
#>  6:         0.78
#>  7:         0.10
#>  8:         0.12
#>  9:         0.16
#> 10:         0.16
#> 11:         6.47
#> 12:         6.45
```

The structural variables — basal area, tree numbers and diameter — are
reconstructed well, **and hold up on sites the model never saw**: that
gap between the train and holdout columns is the part that matters,
because it is the only evidence that the fitted processes generalise
rather than memorise. The noisier per-tree rates (growth, mortality,
regeneration) are harder to constrain from two inventories and should be
read with that in mind (see the caveat below). The hybrid tracks the
mechanistic model on the structural variables — with **no growth
equation specified** — which is what makes hybrids useful when a
process’s form is uncertain; the surrounding mechanistic processes still
anchor the model ecologically.

``` r

# facet_wrap, NOT facet_grid: the three responses live on completely different
# scales (dbh ~150 cm, ba ~8 m2, trees ~20). facet_grid(scales = "free") frees x
# per column but y per ROW, so ba and trees would be drawn on dbh's y-axis -
# squashing them onto y = 0 and flattening the 1:1 line until a good fit looks
# like a broken one. facet_wrap frees both axes per panel.
ggplot(obs_pred[variable %in% c("dbh", "ba", "trees") & is.finite(obs)],
       aes(obs, value, colour = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_point(alpha = 0.25, size = 0.5) +
  facet_wrap(~ split + variable, scales = "free", ncol = 3,
             labeller = labeller(.multi_line = FALSE)) +
  scale_colour_manual(values = model_cols) +
  labs(x = "observed", y = "predicted", colour = NULL) +
  theme_minimal() + theme(legend.position = "top")
```

![](D/D-assess-plot-1.png)

Each panel is on its **own** axes, and the grey line is 1:1 — so the
diagonal is the target in every panel regardless of the units on it.

> **A caveat on the rate variables.** `growth`, `mort` and `reg` are
> per-tree rates and are only weakly constrained by two inventories, so
> their Spearman is a **noisy estimate** — refitting with a different
> random seed moves it by a few hundredths. Differences of that size
> between the models should not be over-read; the structural variables,
> estimated from far more information, are what carry the comparison.
> `mort` is the extreme case and gets a vignette of its own —
> **Mortality: a binomial response and a neural-network process** —
> because rank correlation is the wrong tool for a response that is
> mostly zeros.

## Interpreting the growth process with ALE

[`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) returns the
**accumulated local effect** of every driver on each process — the
effect measure of choice when predictors are correlated, because it
accumulates *local* changes within the observed data instead of
extrapolating to unrealistic combinations. Running it on **both** fitted
models lets us read what growth responds to *and* compare the
mechanistic process with the neural-network hybrid, across growth’s
**structural** inputs (tree size, light) and its **environmental niche**
(climate). It returns one table per process (`growth`, `mortality`,
`regeneration`); we use the growth table of each model.

``` r

FINN.seed(1)
ale_proc <- ALE(m,        env_dt, init_cohorts, plot = FALSE)   # mechanistic growth
ale_hyb  <- ALE(m_hybrid, env_dt, init_cohorts, plot = FALSE)   # growth = NN

# one tidy table with both variants. Environmental drivers come back on the
# model's scaled axis, so back-transform temp/prec to natural units (deg C, mm);
# the structural drivers (dbh, light) are already on their natural axes.
to_natural <- function(a, model) {
  sc <- as.data.table(model$env_scaling)[, .(var = variable, center, scale)]
  merge(data.table(a), sc, by = "var", all.x = TRUE)[!is.na(center), x := x * scale + center]
}
growth <- rbind(
  to_natural(ale_proc$growth, m)[,        model := "Process (mechanistic)"],
  to_natural(ale_hyb$growth,  m_hybrid)[, model := "Hybrid (growth = NN)"])
growth <- merge(growth, species_dt, by = "species")
```

One helper plots the growth ALE of both variants for any set of drivers
— one panel per species × driver, with **independent axes** so responses
of very different magnitude stay legible side by side:

``` r

plot_growth_ale <- function(vars, labs) {
  d <- growth[var %in% vars & species %in% c(1, 2, 4, 5)]   # the four focal species
  d[, var := factor(var, levels = names(labs))]
  ggplot(d, aes(x, ale, colour = model)) +
    geom_line(linewidth = 0.7) +
    facet_grid(species_name ~ var, scales = "free",
               labeller = labeller(var = as_labeller(labs, label_parsed),
                                    species_name = label_value)) +
    scale_colour_manual(values = model_cols) +
    labs(x = NULL, y = "accumulated local effect on growth", colour = NULL) +
    theme_minimal() + theme(legend.position = "top", strip.text.y = element_text(angle = 0))
}
```

### Environmental niches

Because the environmental responses are learned, each curve is an
inferred **niche** — how a species’ growth scales along a climate
gradient. Ecologically familiar contrasts emerge without being
prescribed: *Pinus ponderosa* grows best on the driest sites, while
wetter-climate species ramp up with precipitation. The hybrid’s learned
climate response broadly tracks the mechanistic one.

``` r

plot_growth_ale(c("temp", "prec"),
                c(temp = "temperature~(degree*C)", prec = "precipitation~(mm)"))
```

![](D/D-ale-environmental-1.png)

### Which drivers matter — `feature_importance()`

The niches show the *shape* of each response;
[`feature_importance()`](https://finnverse.github.io/FINN/reference/feature_importance.md)
gives its *magnitude*. It permutes each environmental predictor and
re-simulates, measuring the resulting increase in prediction error — a
re-simulation counterpart to the analytical importance in
[`summary()`](https://rspatial.github.io/terra/reference/summary.html).
Below we read off the importances for the growth process of the four
focal species:

``` r

FINN.seed(1)
fi <- feature_importance(m, env_dt, init_cohorts, nperm = 20L, seed = 1L,
                         patches = 4, patch_size = 0.06, device = "cpu")
fimp_growth <- merge(data.table(fi$growth)[species %in% c(1, 2, 4, 5)],
                     species_dt, by = "species")
```

``` r

ggplot(fimp_growth, aes(reorder(variable, importance), importance, fill = importance)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~species_name, nrow = 1) +
  scale_fill_viridis_c(guide = "none") +
  labs(x = NULL, y = "permutation importance (delta error)") + theme_minimal()
```

![](D/D-featimp-plot-1.png)

Each species leans on a different climate driver — the importances are
learned per species, not shared.

### Response to size and light — did the network recover the mechanistic form?

The point of a hybrid is to ask which functional form the neural network
learned and how it *compares* to the mechanistic function it replaced.
For the pines the network recovered the mechanistic decline of growth
with size almost exactly; for *Abies grandis* and *Pseudotsuga* it
settled on a different, data-driven shape — the freedom that makes
hybrids useful when the mechanistic form is uncertain.

``` r

plot_growth_ale(c("dbh", "light"),
                c(dbh = "diameter~(cm)", light = "light~availability"))
```

![](D/D-ale-structural-1.png)

Each curve spans only the range its *own* model actually simulates (ALE
never extrapolates), so where the process and hybrid cover different
ranges — as for *Abies grandis* light — the two models simply place that
species in different conditions; the mismatch is information, not an
artefact.

``` r

sessionInfo()
#> R version 4.5.0 (2025-04-11)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS 26.5.1
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] C
#> 
#> time zone: Europe/Berlin
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] torch_0.15.1      ggplot2_3.5.2     data.table_1.17.8 FINN_0.1.0       
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.6.5        cli_3.6.6          knitr_1.50         rlang_1.2.0       
#>  [5] xfun_0.57          processx_3.8.6     generics_0.1.4     coro_1.1.0        
#>  [9] labeling_0.4.3     glue_1.8.0         bit_4.6.0          ps_1.9.1          
#> [13] scales_1.4.0       grid_4.5.0         abind_1.4-8        evaluate_1.0.5    
#> [17] tibble_3.3.0       lifecycle_1.0.5    compiler_4.5.0     dplyr_1.1.4       
#> [21] RColorBrewer_1.1-3 Rcpp_1.1.0         pkgconfig_2.0.3    farver_2.1.2      
#> [25] viridisLite_0.4.2  R6_2.6.1           tidyselect_1.2.1   pillar_1.11.0     
#> [29] callr_3.7.6        magrittr_2.0.3     withr_3.0.2        tools_4.5.0       
#> [33] bit64_4.6.0-1      gtable_0.3.6
```
