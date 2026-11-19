# Regeneration with an effect of conspecific adult density

A variant of
[regeneration_adult](https://finnverse.github.io/FINN/reference/regeneration_adult.md)
in which the density of conspecific adults is an EFFECT on recruitment
rather than a gate: a species without adults on the site still recruits
a share `f0` of its potential, so bare ground can be colonised from an
empty stand.

## Usage

``` r
regeneration_adult_effect(
  species,
  parReg,
  pred,
  light,
  dbh = NULL,
  trees = NULL,
  debug = FALSE
)
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

- dbh:

  torch.Tensor cohort diameters, `[sites, patches, cohorts]`.

- trees:

  torch.Tensor stems per cohort, `[sites, patches, cohorts]`.

- debug:

  logical If TRUE, return the intermediate components as a list.
  Defaults to FALSE.

## Value

torch.Tensor mean recruitment, stems ha^-1 step^-1,
`[sites, patches, species]`. With `debug = TRUE`, a list that also holds
`adults` (A, stems ha^-1, `[sites, 1, species]`), `f0` and `seed`
(f(A)).

## Details

The conspecific adult density `A_k` (stems ha^-1) is computed as in
[regeneration_adult](https://finnverse.github.io/FINN/reference/regeneration_adult.md):
the stems of species `k` weighted by a soft adult threshold
`sigmoid((d - x_k) / tau)` (`tau` = 2 dbh units), pooled over the
patches of the site. The adult term is \$\$f(A_k) = f\_{0,k} + (1 -
f\_{0,k}) \frac{A_k}{A_k + A\_{1/2,k}}\$\$ with \\f_0 =
\mathrm{sigmoid}(\mathrm{reg\\adult\\logitf0}) \in (0, 1)\\, and the
mean recruitment is \$\$m = r(\mathrm{light}) \\ e^{\mathrm{env}} f(A) +
\mathrm{reg\\floor},\$\$ capped by the Beverton-Holt term of
[regeneration_capped](https://finnverse.github.io/FINN/reference/regeneration_capped.md),
`K m / (K + m)`.

`f0` is the share of the potential recruitment a species gets on a site
without conspecific adults – seed rain from outside the site, advance
regeneration. It rises to 1 as adults become abundant. `f0 -> 1` means
no adult effect (the process becomes
[regeneration_capped](https://finnverse.github.io/FINN/reference/regeneration_capped.md));
`f0 -> 0` recovers the gate of
[regeneration_adult](https://finnverse.github.io/FINN/reference/regeneration_adult.md),
where a species without adults recruits only at `reg_floor`.

The soft threshold counts small stems partly: a 1-cm stem is about 1.5 %
of an adult at `x = 9.4` cm, so a low `A_half` gives a non-zero adult
term even in a stand of seedlings.

The cap `reg_logK` is created by
[`finn()`](https://finnverse.github.io/FINN/reference/finn.md) as for
[regeneration_capped](https://finnverse.github.io/FINN/reference/regeneration_capped.md)
(resolution from
[`createProcess(cap = )`](https://finnverse.github.io/FINN/reference/createProcess.md),
prior from `prior`/`prior_lambda`), so it need not be declared. Declare
the three adult parameters with
[`createProcess(custom_parameters = )`](https://finnverse.github.io/FINN/reference/createProcess.md);
a length-1 init gives one shared value and a length-`N_species` init one
value per species:

- `reg_adult_logx`: the log adult threshold, `x = exp(reg_adult_logx)`
  in the units of `dbh`;

- `reg_adult_logA`: the log half-saturation adult density,
  `A_half = exp(reg_adult_logA)` stems ha^-1;

- `reg_adult_logitf0`: the logit of the share without adults,
  `f0 = plogis(reg_adult_logitf0)`.

As for
[regeneration_adult](https://finnverse.github.io/FINN/reference/regeneration_adult.md),
`forward()` passes `dbh` and `trees` because the signature has them.
Called without `dbh`/`trees` (as the ALE code in `xAI.R` does), the
function returns recruitment with the seed source present (`f = 1`).

## See also

[regeneration_adult](https://finnverse.github.io/FINN/reference/regeneration_adult.md),
[regeneration_capped](https://finnverse.github.io/FINN/reference/regeneration_capped.md),
[`createProcess()`](https://finnverse.github.io/FINN/reference/createProcess.md),
[recruitment](https://finnverse.github.io/FINN/reference/recruitment.md)
for the recruitment options and how to combine them

## Examples

``` r
library(data.table)
# three sites of the bundled Oregon FIA sample (?fia_obs_dt)
sites <- 1:3
Nsp   <- max(fia_obs_dt$species)
init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)

FINN.seed(1)
m <- finn(
  N_species = Nsp, recruits_dbh = 12.9,
  regeneration_process = createProcess(
    ~ temp + prec, FINN::regeneration_adult_effect,
    custom_parameters = list(reg_adult_logx    = rep(log(20), Nsp),   # adult from 20 cm
                             reg_adult_logA    = rep(log(20), Nsp),   # half effect at 20 adults/ha
                             reg_adult_logitf0 = qlogis(0.2)),        # 20 % without adults, shared
    optimizeSpecies = TRUE, optimizeEnv = TRUE))
fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
    init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06, plot_progress = FALSE)
#> ℹ weights: intercept-only baseline x your vector (dbh=0.00243, ba=12, trees=0.968, growth=1.84e+04, mortality=38.5, regeneration=0.768)
#> ℹ share of the objective at the intercept-only baseline (scaled): dbh=17%, ba=17%, trees=17%, growth=17%, mortality=17%, regeneration=17%
# the fitted share of the potential recruitment without conspecific adults
plogis(as.numeric(m$reg_adult_logitf0))
#> [1] 0.1968171
```
