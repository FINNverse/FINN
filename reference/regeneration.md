# Regeneration without a recruitment cap

FINN's uncapped mechanistic regeneration process: the mean number of
recruits a species contributes to one patch in one simulation step, in
stems per hectare, `mean = regP * exp(pred) + reg_floor`. `regP` is the
patch's light response, a sigmoid in the available light centred on the
species' light requirement `parReg` and shifted to be 0 without light,
so a shade-tolerant species (small `parReg`) recruits under a closed
canopy and a light-demanding one does not. `exp(pred)` is the species'
environmental response, and `reg_floor`
([`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md),
default `1e-3`) is an additive floor that keeps a species from being
lost for good.

The mean is unbounded in `pred`, so a steep climate coefficient can turn
a favourable site into a recruitment burst of hundreds of stems in one
year.
[`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md)
adds a Beverton-Holt cap and is the default regeneration process since
0.2.0; this function is the uncapped process of Pichler & Käber (2026)
and stays available as the opt-in.

## Usage

``` r
regeneration(species, parReg, pred, light, debug = FALSE)
```

## Arguments

- species:

  torch.Tensor species information.

- parReg:

  torch.Tensor Regeneration parameters. 0 \<= parReg \<= 1 This
  parameter denotes the fraction of light needed for a species to
  regenerate. In general low values for high regeneration and high
  values for low regeneration.

- pred:

  torch.Tensor Prediction values.

- light:

  torch.Tensor Available light variable for calculation.

- debug:

  logical If TRUE, return the intermediate components as a list.
  Defaults to FALSE.

## Value

torch.Tensor Mean recruitment per species, stems ha^-1 step^-1, shaped
`[sites, patches, species]` – a per-HECTARE rate, not a per-patch count:
`forward()` multiplies it by `patch_size_ha` itself before the recruit
draw. It includes the additive `reg_floor` and is unbounded above.

## See also

[`regeneration_capped()`](https://finnverse.github.io/FINN/reference/regeneration_capped.md)
(the default regeneration process since 0.2.0),
[`regeneration_adult()`](https://finnverse.github.io/FINN/reference/regeneration_adult.md),
[recruitment](https://finnverse.github.io/FINN/reference/recruitment.md)
for the recruitment options and how to combine them

## Examples

``` r
# the uncapped regeneration of the paper instead of the default
# regeneration_capped(); recruits enter at 1 cm (recruits_dbh)
m <- finn(N_species = 2, recruits_dbh = 1,
          regeneration_process = createProcess(~ 1 + temp, regeneration))

# a 30-year simulation from bare ground; an unfitted model uses `env` as
# given (no scaling constants yet), so the predictor is passed standardised
env <- data.frame(siteID = 1L, year = 1:30, temp = 0)
sim <- simulateForest(m, env = env, patches = 5, patch_size = 0.1)
tail(sim$wide$site)
#>    siteID species  year      dbh        ba trees    growth       mort   reg
#>     <int>   <int> <int>    <num>     <num> <num>     <num>      <num> <num>
#> 1:      1       2    25 45.65352 0.3339909   1.4 0.2475858 0.01785392   0.0
#> 2:      1       2    26 34.61153 0.3624062   2.0 0.2108145 0.01796514   0.6
#> 3:      1       2    27 37.13254 0.3924276   2.0 0.3475665 0.01763394   0.0
#> 4:      1       2    28 39.73211 0.4234293   2.0 0.3203096 0.01771779   0.0
#> 5:      1       2    29 42.50273 0.4554290   2.0 0.2906452 0.01779613   0.0
#> 6:      1       2    30 45.51809 0.4888201   2.0 0.2563899 0.01789678   0.0
#>    r_mean_ha
#>        <num>
#> 1: 0.7366613
#> 2: 0.7304133
#> 3: 0.7230549
#> 4: 0.7198278
#> 5: 0.7110088
#> 6: 0.7006934
```
