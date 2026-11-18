# Growth

FINN's mechanistic growth process: the relative diameter growth of a
cohort in one simulation step,
`growth = shade * exp(pred) * exp(-parGrowth[, 2] * dbh)`, so that the
diameter after the step is `dbh * (1 + growth)`. `shade` is the light
response, a sigmoid in the available light centred on the species' light
requirement `parGrowth[, 1]` and rescaled to run from 0 in full shade to
1 in full light; `exp(pred)` is the species' environmental response; and
the last term lets the relative growth decline with size at the species'
rate `parGrowth[, 2]` (per cm).

## Usage

``` r
growth(
  dbh,
  species,
  parGrowth,
  pred,
  light,
  light_steepness = 10,
  debug = FALSE,
  trees = NULL
)
```

## Arguments

- dbh:

  torch.Tensor Diameter at breast height of each cohort, in cm.

- species:

  torch.Tensor (int64) Species index of each cohort.

- parGrowth:

  torch.Tensor Growth parameters per species, two columns: the light
  requirement of the shade response (between 0 and 1) and the size decay
  (per cm).

- pred:

  torch.Tensor Environmental response per species, on the log scale.

- light:

  torch.Tensor Available light for each cohort (0 to 1), from
  [`competition()`](https://finnverse.github.io/FINN/reference/competition.md).

- light_steepness:

  numeric Steepness of the light-response sigmoid. Defaults to 10.

- debug:

  logical If TRUE, return the intermediate components as a list
  (`shade`, `light`, `environment`, `growth`). Defaults to FALSE.

- trees:

  torch.Tensor (Optional) Number of trees per cohort. Not part of the
  equation; the model passes it for its records. Defaults to NULL.

## Value

torch.Tensor The relative diameter growth of each cohort per simulation
step, `[sites, patches, cohorts]`, non-negative; or, if `debug = TRUE`,
a list of the intermediate components.

## See also

[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md),
[`mortality()`](https://finnverse.github.io/FINN/reference/mortality.md),
which takes this growth as an input,
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md),
[`predictTrees()`](https://finnverse.github.io/FINN/reference/predictTrees.md)
to evaluate the fitted function on observed trees

## Examples

``` r
# growth as the growth process of a model: the equation above with its two
# species parameters, and a linear response to one environmental predictor
m <- finn(N_species = 2, growth_process = createProcess(~ 1 + temp, growth))

# a 30-year simulation from bare ground on one site. A model that was not
# fitted has no scaling constants for `env` yet and uses it as given, so the
# predictor is passed standardised (0 = an average site).
env <- data.frame(siteID = 1L, year = 1:30, temp = 0)
sim <- simulateForest(m, env = env, patches = 5, patch_size = 0.1)
tail(sim$wide$site)
#>    siteID species  year      dbh        ba trees    growth       mort   reg
#>     <int>   <int> <int>    <num>     <num> <num>     <num>      <num> <num>
#> 1:      1       2    25 43.95226 0.4515558   1.8 0.3176936 0.02061903   0.2
#> 2:      1       2    26 35.49084 0.4853425   2.4 0.3634461 0.02086896   0.6
#> 3:      1       2    27 37.95831 0.5201546   2.4 0.5138360 0.02171769   0.2
#> 4:      1       2    28 37.86688 0.5574864   2.6 0.4682905 0.02145828   0.2
#> 5:      1       2    29 41.08892 0.5983314   2.6 0.4132676 0.02119704   0.0
#> 6:      1       2    30 44.88969 0.6459334   2.6 0.3383577 0.02078469   0.0
#>    r_mean_ha
#>        <num>
#> 1: 1.0368954
#> 2: 1.0248361
#> 3: 1.0110019
#> 4: 0.9956975
#> 5: 0.9801638
#> 6: 0.9627833
```
