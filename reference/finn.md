# Forest Informed Neural Network

Creates a Forest Informed Neural Network (FINN), a differentiable,
cohort-based dynamic forest (gap) model in the tradition of
JABOWA/ForClim-style models, in which any of the four demographic
processes (growth, mortality, regeneration, competition for light) can
either be specified mechanistically or replaced by a deep neural network
(DNN). Mechanistic and DNN-based processes are calibrated jointly,
end-to-end, via gradient descent (see
[fit](https://finnverse.github.io/FINN/reference/fit.md)).

## Usage

``` r
finn(
  N_species,
  mortality_process = NULL,
  growth_process = NULL,
  regeneration_process = NULL,
  competition_process = NULL,
  recruits_dbh = 1
)
```

## Arguments

- N_species:

  (`integer(1)`)  
  Number of species.

- mortality_process:

  (`function`)  
  Mortality process, created by
  [createProcess](https://finnverse.github.io/FINN/reference/createProcess.md)
  (mechanistic) or
  [createHybrid](https://finnverse.github.io/FINN/reference/createHybrid.md)
  (DNN-based). If `NULL` (default), the default mechanistic mortality
  process is used.

- growth_process:

  (`function`)  
  Growth process, created by
  [createProcess](https://finnverse.github.io/FINN/reference/createProcess.md)
  (mechanistic) or
  [createHybrid](https://finnverse.github.io/FINN/reference/createHybrid.md)
  (DNN-based). If `NULL` (default), the default mechanistic growth
  process is used.

- regeneration_process:

  (`function`)  
  Regeneration process, created by
  [createProcess](https://finnverse.github.io/FINN/reference/createProcess.md)
  (mechanistic) or
  [createHybrid](https://finnverse.github.io/FINN/reference/createHybrid.md)
  (DNN-based). If `NULL` (default), the default mechanistic regeneration
  process is used.

- competition_process:

  (`function`)  
  Competition (light availability) process, created by
  [createProcess](https://finnverse.github.io/FINN/reference/createProcess.md).
  If `NULL` (default), the default mechanistic competition process is
  used.

- recruits_dbh:

  (`numeric(1)`)  
  Starting dbh for recruits. Has value 1.0 as default.

## Details

FINN represents the forest as cohorts of trees, grouped by site, patch
and cohort, each characterized by diameter at breast height (dbh),
number of trees, and species identity. Starting from an initial state,
FINN simulates the forest forward in discrete annual time steps by
sequentially applying the four demographic processes: competition (light
availability, based on basal area and species-specific shading), growth
(diameter increment as a function of light, size and environment),
mortality (binomial death of trees as a function of growth, light, size
and environment) and regeneration (recruitment of new cohorts as a
function of light and environment, drawn from a negative binomial
distribution). Each process additionally depends on a species- and
process-specific environmental-response function that maps site-level
environmental predictors to a scalar effect on the process.

Every process and its environmental-response function can be configured
in one of two ways: (1) mechanistically, with an explicit functional
form and interpretable parameters (e.g. light-response thresholds,
allometric coefficients), created via
[createProcess](https://finnverse.github.io/FINN/reference/createProcess.md);
or (2) as a hybrid process, in which the environmental-response function
or the entire process equation is replaced by a DNN, created via
[createHybrid](https://finnverse.github.io/FINN/reference/createHybrid.md).
The remaining mechanistic processes constrain the DNN to ecologically
plausible behaviour, while the DNN absorbs misalignments and structural
simplifications that would otherwise bias the mechanistic processes;
both are estimated jointly rather than calibrated in isolation and
plugged in afterwards. If a process argument is left at its default
`NULL`, the corresponding default mechanistic process (as described in
Pichler & Käber, 2026, Methods in Ecology and Evolution) is used.

`finn()` only assembles the model architecture (analogous to
instantiating a `torch` `nn_module`); none of the process or
environmental-response parameters are estimated yet. Use
[fit](https://finnverse.github.io/FINN/reference/fit.md) to calibrate
the returned object against observed data, and
[predict.finn_class](https://finnverse.github.io/FINN/reference/predict.finn_class.md)/[simulateForest](https://finnverse.github.io/FINN/reference/simulateForest.md)
to simulate forest dynamics from a (fitted) model.

## References

Pichler, M., & Käber, Y. (2026). Inferring processes within dynamic
forest models using hybrid modelling. *Methods in Ecology and
Evolution*.
[doi:10.1111/2041-210x.70347](https://doi.org/10.1111/2041-210x.70347)
