#' Calculate the basal area of a tree given the diameter at breast height (dbh)
#'
#' This function calculates the basal area of a tree given the diameter at breast height (dbh): the
#' area of a circle of that diameter, `pi * (dbh / 200)^2`, in square metres for a `dbh` in centimetres.
#' It is the tensor counterpart of [FINN::dbh2ba()].
#'
#' @param dbh torch.Tensor (or numeric vector) The diameter at breast height of the tree, in cm.
#'
#' @return An object of the type of `dbh` (a torch.Tensor for a tensor) with the basal area of one
#'   stem, in m^2.
#'
#' @seealso [FINN::BA_stand()] for the basal area per hectare of a cohort, [FINN::dbh2ba()]
#'
#' @examplesIf torch::torch_is_installed()
#' # one stem of 10, 30 and 50 cm: 0.008, 0.071 and 0.196 m^2
#' dbh = torch::torch_tensor(c(10, 30, 50))
#' basal_area = BA_stem(dbh)
#' print(basal_area)
#'
#' # the same numbers without torch
#' BA_stem(c(10, 30, 50))
#'
#' @import torch
#' @export
BA_stem = function(dbh) {
  return(pi*(dbh/100./2.)^2.0)
}


#' Calculate the Basal Area of a Stand
#'
#' This function calculates the basal area of a stand based on the diameter at breast height (dbh), the number of trees, and the patch size in hectares.
#'
#' @param dbh A torch tensor or numeric vector representing the diameter at breast height of the trees in centimeters.
#' @param trees A torch tensor or numeric vector representing the number of trees.
#' @param patch_size_ha A numeric value representing the size of the patch in hectares.
#'
#' @details
#' The basal area of a stand is the cross-sectional area of all trees in a stand per unit area. This function calculates the basal area per ha using the formula:
#' \deqn{BA = \left( \frac{\pi \left( \frac{\text{dbh}}{100} \right)^2}{4} \right) \times \text{trees} \div \text{patch\_size\_ha}}
#'
#' The formula takes into account the diameter at breast height (dbh) in centimeters, the number of trees, and the size of the patch in hectares to calculate the basal area in square meters per hectare.
#'
#' This plot illustrates the basal area for different combinations of dbh and number of trees.
#'
#' \figure{BA_stand_plot2.png}{dbh, trees, basal area}
#'
#' Sensitivity of basal area for different combinations of dbh, number of trees to patch size.
#'
#' \figure{BA_stand_plot1.png}{Patch size, dbh, trees, basal area}
#'
#' @return An object of the type and shape of `dbh` (a torch.Tensor for tensors): the basal area of
#'   each cohort in square meters per hectare.
#' @seealso [FINN::BA_stem()] for one stem, [FINN::competition()], which turns the basal area of the
#'   taller cohorts into the light a cohort receives
#' @examplesIf torch::torch_is_installed()
#' # three cohorts on a patch of 0.1 ha: two 50-cm trees, ten 20-cm trees and
#' # forty 5-cm trees
#' cohorts <- data.frame(siteID = 1, patchID = 1, cohortID = 1:3, species = 1,
#'                       dbh = c(50, 20, 5), trees = c(2, 10, 40))
#' cohort <- CohortMat(obs_df = cohorts)
#' ba <- BA_stand(cohort$dbh, cohort$trees, patch_size_ha = 0.1)
#' ba          # m^2/ha per cohort, [sites, patches, cohorts]
#' ba$sum()    # the stand: 7.9 m^2/ha
#'
#' # the same trees on a patch half the size are twice the basal area per hectare
#' BA_stand(cohort$dbh, cohort$trees, patch_size_ha = 0.05)$sum()
#'
#' # numeric vectors work as well
#' BA_stand(dbh = c(50, 20, 5), trees = c(2, 10, 40), patch_size_ha = 0.1)
#' @export
BA_stand <- function(dbh, trees, patch_size_ha) {
  return((pi * (dbh / 100 / 2)^2 * trees) / patch_size_ha)
}



#' Calculate the height of a tree based on its diameter at breast height and an allometry parameter
#'
#' @param dbh A numeric value representing the diameter at breast height of the tree in cm.
#' @param parHeight A numeric value representing the species height allometry.
#'
#' @details
#'
#' This function calculates the height of a tree based on the diameter at breast height (dbh) and a parameter parHeight.
#'
#' The height is calculated using the formula:
#' \deqn{height = \left( \exp \left( \frac{(\text{dbh} \times \text{parHeight})}{(\text{dbh} + 100)} \right) - 1 \right) \times 100 + 0.001}
#' where dbh is the diameter at breast height of the tree in cm and parHeight is an allometric species specific parameter.
#'
#' As the diameter grows, the height approaches \eqn{(\exp(\text{parHeight}) - 1) \times 100}{(exp(parHeight) - 1) * 100}
#' m: 35 m for parHeight = 0.3, 65 m for 0.5, 101 m for 0.7. When the height
#' parameter is fitted (the first column of the competition parameters), FINN keeps
#' it between 0.3 and \eqn{\log(2.3) \approx 0.83}{log(2.3) = 0.83}, so no species
#' can grow taller than 130 m, the physical limit of tree height (Koch et al. 2004).
#' Values below 0.3 suit small trees and shrubs; set `createProcess(lower = )` for
#' them.
#'
#' @references Koch, G. W., Sillett, S. C., Jennings, G. M. & Davis, S. D. (2004).
#'   The limits to tree height. Nature 428, 851-854.
#'
#' \figure{height_plot1.png}{Parameter range}
#'
#' @return A numeric vector (or torch.Tensor, for tensor input) with the height of the tree in m.
#' @seealso [FINN::competition()], which ranks the cohorts of a patch by this height
#' @examples
#' # a 30-cm tree of a species with parHeight = 0.5 is about 12 m tall
#' height(30, 0.5)
#' # the same diameter for a tall-growing and a short-growing species
#' height(c(30), c(0.5,0.3))
#' # two diameters of one species
#' height(c(30,20), c(0.5))
#'
#' # the height curve of three species; the asymptote is (exp(parHeight) - 1) * 100 m
#' dbh <- seq(0, 150, by = 10)
#' round(sapply(c(low = 0.3, mid = 0.5, tall = 0.7), function(p) height(dbh, p)), 1)
#'
#' @export
height = function(dbh, parHeight) {
  height = (exp((((dbh * parHeight) / (dbh+100))))-1)*100 + 0.001
  return(height)
}


#' Compute the fraction of available light (light) for each cohort based on the given parameters
#'
#' This function calculates the fraction of available light for each cohort of trees based on their diameter at breast height (dbh), species, number of trees, and global parameters.
#'
#' The light of a cohort is one minus the weighted basal area of the cohorts of its patch that are
#' taller than it (by more than 0.1 m), `light = max(0, 1 - sum(BA_taller * parComp[, 2] * 0.1))`,
#' where `BA` is the basal area in m^2/ha ([FINN::BA_stand()]) and the height comes from
#' [FINN::height()] with the species' height parameter `parComp[, 1]`. With a competition
#' strength of 0.2 (about the start value) a cohort under 50 m^2/ha of taller trees receives no light. In a model the function
#' is the competition process, `createProcess(func = competition)`; it takes no environment.
#'
#' @param dbh torch.Tensor Diameter at breast height for each cohort, in cm, `[sites, patches, cohorts]`.
#' @param species torch.Tensor (int64) species index for each cohort, `[sites, patches, cohorts]`.
#' @param trees torch.Tensor Number of trees in each cohort, `[sites, patches, cohorts]`.
#' @param parComp torch.Tensor Competition parameters, one row per species and two columns: the
#'   height parameter of [FINN::height()] and the competition strength.
#' @param h (Optional) Defaults to NULL: the light each cohort receives at its own height. Any
#'   non-`NULL` value returns the light at the forest floor instead, one value per patch
#'   (`[sites, patches, 1]`), which is what the regeneration process is given; the value itself is
#'   not used.
#' @param patch_size_ha numeric Patch size in hectares.
#' @param ba torch.Tensor (Optional) Pre-computed basal area. Defaults to NULL.
#' @param cohortHeights torch.Tensor (Optional) Pre-computed cohort heights. Defaults to NULL.
#' @param n_quantiles integer Number of height quantiles used when `continuous = FALSE`. Defaults to 10.
#' @param continuous logical Use the continuous competition formulation. Defaults to FALSE.
#' @param light_mask logical Take the height quantiles over the live cohorts only
#'   (`trees > 0.5`), leaving every cohort's real height as the bin query. Defaults
#'   to TRUE since 0.2.0: unmasked, the padded cohort width is set by the widest
#'   site in the batch, so the same stand's light depends on which other sites
#'   happen to share its batch. `FALSE` reproduces earlier results. See
#'   [FINN::createProcess()] for the measurement.
#'
#' @return torch.Tensor Fraction of available light (between 0 and 1) for each cohort,
#'   `[sites, patches, cohorts]`; with `h` given, the light at the forest floor, `[sites, patches, 1]`.
#'
#' @seealso [FINN::height()], [FINN::BA_stand()], [FINN::createProcess()] (`n_quantiles`,
#'   `continuous`, `light_mask`), [FINN::finn()]
#'
#' @examplesIf torch::torch_is_installed()
#' # the light of three cohorts on one 0.1 ha patch: two 50-cm trees, ten 20-cm
#' # trees and forty 5-cm trees of one species (height parameter 0.5,
#' # competition strength 0.2)
#' shape   <- c(1, 1, 3)
#' dbh     <- torch::torch_tensor(array(c(50, 20, 5), shape))
#' trees   <- torch::torch_tensor(array(c(2, 10, 40), shape))
#' species <- torch::torch_tensor(array(1L, shape), dtype = torch::torch_long())
#' parComp <- torch::torch_tensor(matrix(c(0.5, 0.2), nrow = 1))
#' # the tallest cohort is in full light, each smaller one is shaded by the
#' # basal area above it
#' competition(dbh, species, trees, parComp, patch_size_ha = 0.1, continuous = TRUE)
#' # the light that reaches the forest floor, where the recruits start
#' competition(dbh, species, trees, parComp, h = 0, patch_size_ha = 0.1)
#'
#' # the light model as the competition process of a model. It takes no
#' # environment; `continuous = TRUE` compares every pair of cohorts instead of
#' # binning them into height quantiles
#' m <- finn(N_species = 2,
#'           competition_process = createProcess(func = competition, continuous = TRUE))
#'
#' # a 30-year simulation from bare ground; an unfitted model uses `env` as
#' # given (no scaling constants yet), so the predictor is passed standardised
#' env <- data.frame(siteID = 1L, year = 1:30, temp = 0)
#' sim <- simulateForest(m, env = env, patches = 5, patch_size = 0.1)
#' tail(sim$wide$site)
#'
#' @import torch
#' @export
competition = function(dbh, species, trees, parComp, h = NULL, patch_size_ha, ba = NULL, cohortHeights = NULL, n_quantiles = 10, continuous = FALSE, light_mask = TRUE){
  parHeight = parComp[,1]
  parCompStr = parComp[,2]
  if(is.null(ba)) ba = BA_stand(dbh = dbh, trees = trees, patch_size_ha = patch_size_ha)*parCompStr[species]*0.1
  if(is.null(cohortHeights)) cohortHeights = height(dbh, parHeight[species])$unsqueeze(4)

  if(is.null(h)) {
    if(continuous) {
      h = cohortHeights
      ba_height = (ba$unsqueeze_(4)$multiply(((cohortHeights - h$permute(c(1,2, 4, 3)) - 0.1)/1e-2)$sigmoid_() ))$sum(-2)
    } else {
      device = dbh$device
      quants = torch::torch_linspace(0, 1, steps = n_quantiles+1, dtype = dbh$dtype, device = device)
      # Quantiles and their indices. `light_mask` takes the bin edges over the LIVE cohorts
      # only: the cohort dimension is padded to the widest site x patch of the batch, and an
      # empty slot has height 0.001, so with many of them every live cohort lands in the top
      # bin, where nothing shades anything. Only the quantile INPUT is masked; the
      # searchsorted query keeps every cohort's real height, so a dead-but-tracked cohort is
      # still placed in a bin and still receives light.
      hq = if(isTRUE(light_mask))
             torch::torch_where(trees$gt(0.5)$unsqueeze(4), cohortHeights,
                                torch::torch_full_like(cohortHeights, NaN))
           else cohortHeights
      cohortHeights_quant = if(isTRUE(light_mask))
        # a patch whose cohorts are all dead has an all-NaN slice and no quantile;
        # its basal area is zero, so edges of 0 put everything in the top bin and
        # give the empty patch full light, which is the right answer.
        torch::torch_nanquantile(hq, q = quants, dim = 3)$squeeze(4)$permute(c(2, 3, 1))$nan_to_num(nan = 0)
        else torch::torch_quantile(hq, q = quants, dim = 3)$squeeze(4)$permute(c(2, 3, 1))
      indices = torch::torch_searchsorted(cohortHeights_quant$narrow(dim = 3, start = 2L, length = n_quantiles - 1L)$contiguous(), cohortHeights$squeeze(4), right = TRUE) + 1L
      # Calculation of comp
      cohortHeights_quant = cohortHeights_quant$unsqueeze(4)
      h_quant = cohortHeights_quant
      threshold_indices = ( (cohortHeights_quant - h_quant$permute(c(1,2, 4, 3)) - 0.1)  /1e-2)$sigmoid_()
      # ba must be grouped by indices
      ba_accum = torch_zeros(ba$shape[1], ba$shape[2], quants$shape, device = device)$scatter_add(3, index=indices, src = ba)
      ba_height_accum = (ba_accum$unsqueeze(4)$multiply( threshold_indices )) $sum(-2)
      ba_height = torch_gather(ba_height_accum,3, index = indices)
    }


  }else{
    ba_height = (ba$unsqueeze_(4)$multiply_(((cohortHeights - 0.1)/1e-2)$sigmoid_() ))$sum(-2)
  }
  light = 1-ba_height
  light = torch_clamp(light, min = 0)
  return(light)
}

#' Mortality
#'
#' @description
#' FINN's mechanistic mortality process: the probability that a tree of a cohort
#' dies in one simulation step. It is a logistic function of four terms,
#' `p = sigmoid(pred + light * parMort[, 1] + parMort[, 2] * dbh/100 + growth/growth_scale * parMort[, 3])`,
#' with one coefficient per species for each of shade, size and growth. `pred` is
#' the species' environmental response on the logit scale, so a species with no
#' environmental predictors dies at a constant rate. Growth enters as the rate the
#' growth process produced in the same step, which is what makes mortality depend
#' on a tree's recent performance rather than on its size alone.
#'
#' The probability is clamped to `[mortality_floor, mortality_ceiling]` (default
#' `0.001` / `0.999` per simulation step) and the number of deaths is then drawn
#' with `mortality_draw`; both are settings of [FINN::createProcess()].
#'
#' @param dbh torch.Tensor Diameter of each cohort, in cm.
#' @param species torch.Tensor Species index of each cohort.
#' @param trees torch.Tensor Number of trees in each cohort.
#' @param parMort torch.Tensor Mortality parameters per species, three columns:
#'   the shade, size and growth coefficients of the logit.
#' @param pred torch.Tensor Environmental response per species, on the logit scale.
#' @param light torch.Tensor Available light for each cohort.
#' @param base_steepness numeric Not used by this function; kept in the signature for custom
#'   mortality functions written against it.
#' @param debug logical If TRUE, return a list with the light (`light`), the mortality
#'   probability before any further step (`environment`) and the returned probability
#'   (`predM`, identical to `environment`). Defaults to FALSE.
#' @param growth torch.Tensor (Optional) Growth entering the mortality response; defaults to the model's current growth.
#' @param growth_scale numeric Growth enters the logit as `growth / growth_scale`.
#'   The default 0.01 expresses the annual relative growth in %/yr (about one
#'   standard deviation of observed growth), so that the growth coefficient
#'   `parMort[, 3]` acts per 1 %/yr and its bounds of +-10 allow a meaningful
#'   effect. With the raw rate (`growth_scale = 1`, the behaviour before 0.2.0)
#'   the term spans only a few tenths on the logit between slow and fast trees
#'   even at the bound. Set it for a fitted model through
#'   `createProcess(func = FINN::mortality, growth_scale = 1)`, not through a
#'   wrapper function: the kernel reads `self$record_raws` from its own
#'   environment and a wrapper leaves that unbound.
#'
#' @return A `torch` tensor of per-cohort mortality probabilities per simulation step,
#'   `[sites, patches, cohorts]`.
#'
#' @seealso [FINN::createProcess()] (`mortality_draw`, `mortality_floor`, `growth_scale`),
#'   [FINN::growth()], [FINN::finn()], [FINN::predictTrees()] to evaluate the fitted function on
#'   observed trees
#'
#' @section The starting value:
#' The environmental response is `sigmoid(pred)`. Before any fitting it starts
#' at about 2 % mortality per time step: the intercept of the response is
#' `qlogis(0.02)` and its other weights are 0 (with `createProcess(hidden = )`
#' the bias of the output layer carries the intercept), and the species terms
#' for light, size and growth move it by at most about one logit unit. Until
#' 0.2.0 the weights started at torch's random values, i.e. at `sigmoid(0) =
#' 0.5` per step, which on yearly steps emptied an unfitted stand within a
#' decade. `createProcess(initEnv = )` sets another start; a formula without
#' an intercept (`~ 0 + temp`) leaves torch's start in place.
#'
#' @examplesIf torch::torch_is_installed()
#' # mortality as the mortality process of a model, responding to one predictor;
#' # growth enters in %/yr (growth_scale = 0.01, the default)
#' m <- finn(N_species = 2,
#'           mortality_process = createProcess(~ 1 + temp, mortality, growth_scale = 0.01))
#'
#' # a 30-year simulation from bare ground; an unfitted model uses `env` as
#' # given (no scaling constants yet), so the predictor is passed standardised
#' env <- data.frame(siteID = 1L, year = 1:30, temp = 0)
#' sim <- simulateForest(m, env = env, patches = 5, patch_size = 0.1)
#' tail(sim$wide$site)
#'
#' # the start of the environmental response: qlogis(0.02) and 0 per species
#' as.matrix(m$nn_mortality$parameters[[1]])
#'
#' @export
mortality = function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL, growth_scale = 0.01) {
  if(is.null(growth)) growth = self$g
  if(self$record_raws) {
    self$raw_m = c(self$raw_m,  list(as_array( torch::torch_cat(list(dbh$unsqueeze(4),
                                                                     growth$unsqueeze(4),
                                                                     light$unsqueeze(4),
                                                                     trees$unsqueeze(4),
                                                                     species$unsqueeze(4)$float()), dim = 4)) ))
  }

  environment = torch::torch_sigmoid(pred+growth/growth_scale*parMort[,3][species]+light*parMort[,1][species] +  parMort[,2][species]*(dbh / ( 100)))
  predM = environment
  if(debug == TRUE) out = list(light = light, environment = environment, predM = predM)
  else out = predM
  return(out)
}



mortality_hybrid = function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL) {
  if(is.null(growth)) growth = self$g
  m = self$nn_mortality(dbh = dbh, growth = growth, trees = trees, light = light, species = species, env = pred)$sigmoid()
  if(self$record_raws) {
    self$raw_m = c(self$raw_m,  list(as_array( torch::torch_cat(list(dbh$unsqueeze(4),
                                                                     growth$unsqueeze(4),
                                                                     light$unsqueeze(4),
                                                                     trees$unsqueeze(4),
                                                                     species$unsqueeze(4)$float()), dim = 4)) ))
  }
  return(m)
}



#' Growth
#'
#' @description
#' FINN's mechanistic growth process: the relative diameter growth of a cohort in
#' one simulation step, `growth = shade * exp(pred) * exp(-parGrowth[, 2] * dbh)`,
#' so that the diameter after the step is `dbh * (1 + growth)`. `shade` is the
#' light response, a sigmoid in the available light centred on the species' light
#' requirement `parGrowth[, 1]` and rescaled to run from 0 in full shade to 1 in
#' full light; `exp(pred)` is the species' environmental response; and the last
#' term lets the relative growth decline with size at the species' rate
#' `parGrowth[, 2]` (per cm).
#'
#' @param dbh torch.Tensor Diameter at breast height of each cohort, in cm.
#' @param species torch.Tensor (int64) Species index of each cohort.
#' @param parGrowth torch.Tensor Growth parameters per species, two columns: the
#'   light requirement of the shade response (between 0 and 1) and the size decay
#'   (per cm).
#' @param pred torch.Tensor Environmental response per species, on the log scale.
#' @param light torch.Tensor Available light for each cohort (0 to 1), from
#'   [FINN::competition()].
#' @param light_steepness numeric Steepness of the light-response sigmoid. Defaults to 10.
#' @param debug logical If TRUE, return the intermediate components as a list
#'   (`shade`, `light`, `environment`, `growth`). Defaults to FALSE.
#' @param trees torch.Tensor (Optional) Number of trees per cohort. Not part of
#'   the equation; the model passes it for its records. Defaults to NULL.
#'
#' @return torch.Tensor The relative diameter growth of each cohort per
#'   simulation step, `[sites, patches, cohorts]`, non-negative; or, if
#'   `debug = TRUE`, a list of the intermediate components.
#'
#' @seealso [FINN::createProcess()], [FINN::mortality()], which takes this
#'   growth as an input, [FINN::finn()], [FINN::predictTrees()] to evaluate the
#'   fitted function on observed trees
#'
#' @examplesIf torch::torch_is_installed()
#' # growth as the growth process of a model: the equation above with its two
#' # species parameters, and a linear response to one environmental predictor
#' m <- finn(N_species = 2, growth_process = createProcess(~ 1 + temp, growth))
#'
#' # a 30-year simulation from bare ground on one site. A model that was not
#' # fitted has no scaling constants for `env` yet and uses it as given, so the
#' # predictor is passed standardised (0 = an average site).
#' env <- data.frame(siteID = 1L, year = 1:30, temp = 0)
#' sim <- simulateForest(m, env = env, patches = 5, patch_size = 0.1)
#' tail(sim$wide$site)
#'
#' @import torch
#'
#' @export
growth = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = FALSE, trees = NULL){

  if(self$record_raws) {
    self$raw_g = c(self$raw_g,  list(as_array( torch::torch_cat(list(dbh$unsqueeze(4),
                                                                     light$unsqueeze(4),
                                                                     trees$unsqueeze(4),
                                                                     species$unsqueeze(4)$float()), dim = 4)) ))
  }

  shade = ((1 / (1 + torch::torch_exp(-light_steepness * (light - parGrowth[,1][species]))) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]))) /
         (1 / (1 + torch::torch_exp(-light_steepness * (1 - parGrowth[,1][species]))) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]))))

  environment = torch::torch_exp(pred) # inverse link function
  # growth = shade * environment * (torch::torch_exp(-(dbh / (parGrowth[,2][species] * 100))))
  growth = shade * environment * (torch::torch_exp(-parGrowth[,2][species] * dbh))
  if(debug == TRUE) out = list(shade = shade, light = light, environment = environment,growth = growth) else out = growth
  return(out)
}


growth_hybrid= function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = FALSE, trees = NULL) {

  if(self$record_raws) {
    self$raw_g = c(self$raw_g,  list(as_array( torch::torch_cat(list(dbh$unsqueeze(4),
                                                                     light$unsqueeze(4),
                                                                     trees$unsqueeze(4),
                                                                     species$unsqueeze(4)$float()), dim = 4)) ))
  }

  g = (self$nn_growth(dbh = dbh, trees = trees, light = light, species = species, env = pred)- exp(1))$exp()
  return(g)
}


#' Regeneration without a recruitment cap
#'
#' @description
#' FINN's uncapped mechanistic regeneration process: the mean number of recruits
#' a species contributes to one patch in one simulation step, in stems per
#' hectare, `mean = regP * exp(pred) + reg_floor`. `regP` is the patch's light
#' response, a sigmoid in the available light centred on the species' light
#' requirement `parReg` and shifted to be 0 without light, so a shade-tolerant species
#' (small `parReg`) recruits under a closed canopy and a light-demanding one does
#' not. `exp(pred)` is the species' environmental response, and `reg_floor`
#' ([FINN::createProcess()], default `1e-3`) is an additive floor that keeps a
#' species from being lost for good.
#'
#' The mean is unbounded in `pred`, so a steep climate coefficient can turn a
#' favourable site into a recruitment burst of hundreds of stems in one year.
#' [FINN::regeneration_capped()] adds a Beverton-Holt cap and is the default
#' regeneration process since 0.2.0; this function is the uncapped process of
#' Pichler & Käber (2026) and stays available as the opt-in.
#'
#' @param species torch.Tensor species information.
#' @param parReg torch.Tensor Regeneration parameters. 0 <= parReg <= 1
#' This parameter denotes the fraction of light needed for a species to regenerate.
#' In general low values for high regeneration and high values for low regeneration.
#' @param pred torch.Tensor Prediction values.
#' @param light torch.Tensor Available light variable for calculation.
#' @param debug logical If TRUE, return the intermediate components as a list. Defaults to FALSE.
#'
#' @return torch.Tensor Mean recruitment per species, stems ha^-1 step^-1,
#'   shaped `[sites, patches, species]` -- a per-HECTARE rate, not a per-patch
#'   count: `forward()` multiplies it by `patch_size_ha` itself before the
#'   recruit draw. It includes the additive `reg_floor` and is unbounded above.
#'
#' @seealso [FINN::regeneration_capped()] (the default regeneration process since
#'   0.2.0), [FINN::regeneration_adult()], [recruitment] for the recruitment
#'   options and how to combine them
#'
#' @examplesIf torch::torch_is_installed()
#' # the uncapped regeneration of the paper instead of the default
#' # regeneration_capped(); recruits enter at 1 cm (recruits_dbh)
#' m <- finn(N_species = 2, recruits_dbh = 1,
#'           regeneration_process = createProcess(~ 1 + temp, regeneration))
#'
#' # a 30-year simulation from bare ground; an unfitted model uses `env` as
#' # given (no scaling constants yet), so the predictor is passed standardised
#' env <- data.frame(siteID = 1L, year = 1:30, temp = 0)
#' sim <- simulateForest(m, env = env, patches = 5, patch_size = 0.1)
#' tail(sim$wide$site)
#'
#' @import torch
#' @importFrom torch torch_sigmoid
#' @export
regeneration = function(species, parReg, pred, light, debug = FALSE) {
  if(self$record_raws) {
    self$raw_r = c(self$raw_r,  list(as_array(light$unsqueeze(4))))
  }


  if("matrix" %in% class(pred)) pred = torch::torch_tensor(pred)
  environment = torch::torch_exp(pred) # Environmental inverse link function
  regP = (1 / (1 + torch_exp(-10 * (light - parReg))) - 1 / (1 + torch_exp(10 * parReg))) / (1 - 1 / (1 + torch_exp(10 * (1 - parReg))))
  ## additive floor: createProcess(reg_floor = ); models saved before the argument existed
  ## carry no field and keep the historical 0.2
  reg_floor = if(is.null(self$reg_floor)) 0.2 else self$reg_floor
  mean = (regP*(environment[,NULL])$`repeat`(c(1, species$shape[2], 1)) + reg_floor)
  #regP = torch_sigmoid((light + (1-parReg) - 1)/1e-3) # TODO masking? better https://pytorch.org/docs/stable/generated/torch.masked_select.html
  if(debug == TRUE) out = list(regP = regP, mean = mean) else out = mean
  return(out)
}

regeneration_hybrid = function(species, parReg, pred, light, debug = FALSE) {

  if(self$record_raws) {
    self$raw_r = c(self$raw_r,  list(as_array(light$unsqueeze(4))))
  }

  # clamp the pre-exp logit to avoid float overflow (-> Inf -> NaN gradients)
  # if the net's raw output drifts up during training; exp(10) ~ 22000/ha is
  # already far beyond any plausible recruitment rate
  r = self$nn_regeneration(light = light, species = species, env = pred)$clamp(max = 10)$exp()
}

#' Regeneration with a Beverton-Holt recruitment cap
#'
#' The default regeneration process since 0.2.0, and a drop-in alternative to
#' the uncapped [FINN::regeneration]. It caps the otherwise
#' unbounded mean recruitment with a Beverton-Holt density dependence,
#' `recruits = K * m / (K + m)`, where `m` is the mechanistic mean recruitment
#' (identical to [FINN::regeneration]) and `K` is a fitted carrying capacity
#' (stems ha^-1 step^-1). As `m` grows, recruitment saturates towards `K`.
#' `regeneration_saturation()` is the same function under its older name, kept
#' as a deprecated alias for one release.
#'
#' The cap is the parameter `reg_logK` (`K = exp(reg_logK)`), which **FINN
#' creates itself** -- there is nothing to declare. `createProcess(cap = )`
#' decides its resolution: `"all"` (the default) gives one cap shared by every
#' species, `"species"` one cap per species. Read the fitted value back with
#' `exp(as.numeric(m$reg_logK))`.
#'
#' @section The prior on the cap:
#' A cap that only the data constrain runs to whatever the fitted window
#' supports, which on a national inventory is far above anything observed. So
#' the cap carries a prior, added to the objective as
#' `prior_lambda * (log K - log prior)^2` per cap -- the same form as the L2
#' pull of a staged calibration. It is added **after** the `loss_weights = "auto"`
#' baseline scaling, so `prior_lambda` means the same thing on every dataset.
#'
#' `createProcess(prior = "max_obs_recruits")`, the default, resolves the prior
#' from the data at the start of [FINN::fit()]: the largest observed recruitment
#' rate per timestep, `max(reg / period_length)`. The division matters on any
#' inventory with a remeasurement interval, because `reg` is the SUM over that
#' interval while the cap bounds the mean of one timestep -- taking the column
#' as supplied would place the prior a factor `period_length` too high, and the
#' cap would never bind. A number sets the prior directly, on that same
#' per-timestep scale. `prior_lambda = 0` switches the prior off and leaves the
#' cap to the data alone. Scoring no recruitment at all is an error rather than
#' a silent default, because nothing else in the objective constrains a cap.
#'
#' The resolved value is also the cap's starting value, and it is resolved only
#' on the FIRST [FINN::fit()] of a model: a second call continues from the
#' trained cap instead of resetting it to the prior.
#'
#' With `cap = "species"` the prior is **hierarchical**: each species' log cap is
#' pulled towards the mean of all current log caps rather than towards a fixed
#' value, so the prior penalises how far the species spread apart and leaves the
#' overall level to the data. Every species still starts at the observed maximum.
#'
#' @param species torch.Tensor species information.
#' @param parReg torch.Tensor Regeneration parameters. 0 <= parReg <= 1
#' This parameter denotes the fraction of light needed for a species to regenerate.
#' In general low values for high regeneration and high values for low regeneration.
#' @param pred torch.Tensor Prediction values.
#' @param light torch.Tensor Available light variable for calculation.
#' @param debug logical If TRUE, return the intermediate components as a list. Defaults to FALSE.
#'
#' @return torch.Tensor Mean recruitment, stems ha^-1 step^-1,
#'   `[sites, patches, species]` -- a per-HECTARE rate, not a per-patch count.
#'   `forward()` multiplies it by `patch_size_ha` itself before the recruit draw,
#'   and `reg` is scored on the per-hectare value.
#'
#' @seealso [FINN::regeneration], [FINN::regeneration_adult()],
#'   [FINN::regeneration_adult_effect()],
#'   [FINN::createProcess()], [recruitment] for the recruitment options and how
#'   to combine them
#'
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' Nsp <- max(fia_obs_dt$species)
#' # the default: one cap for all species, prior from the data
#' m <- finn(N_species = Nsp, recruits_dbh = 12.9)
#'
#' # one cap per species, a fixed prior (stems per ha and year), a stronger pull
#' m2 <- finn(N_species = Nsp, recruits_dbh = 12.9,
#'            regeneration_process = createProcess(~ temp + prec, regeneration_capped,
#'                                                 cap = "species", prior = 500, prior_lambda = 1))
#'
#' # the cap is created and resolved by fit(); read it back on the natural scale
#' sites <- 1:3
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#' fit(m2, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 1L, patch_size = 0.06, plot_progress = FALSE)
#' exp(as.numeric(m2$reg_logK))   # one carrying capacity per species
#'
#' @import torch
#' @importFrom torch torch_sigmoid
#' @export
regeneration_capped = function(species, parReg, pred, light, debug = FALSE) {
  if(is.null(self$reg_logK))
    stop("`regeneration_capped` needs the cap parameter `reg_logK`, which finn() ",
         "creates from createProcess(cap = ). This model has none, so it was built ",
         "before 0.2.0 or with a custom func; rebuild it with finn().",
         call. = FALSE)

  if(self$record_raws) {
    self$raw_r = c(self$raw_r,  list(as_array(light$unsqueeze(4))))
  }

  if("matrix" %in% class(pred)) pred = torch::torch_tensor(pred)
  # The cap below saturates, so nothing in the loss holds the network's output down:
  # under Adam it drifts upward until exp() or the product K * mean overflows float32
  # (pred ~ 83-89), which gives NaN gradients and then NaN recruits. At pred = 70 the
  # capped mean is already K to the last float32 digit for any regP > 0 (regP is a
  # difference of two float32 sigmoids, so its smallest non-zero value is ~1e-12) and
  # any K below 1e11 per ha, so the clamp changes no finite result.
  environment = torch::torch_exp(pred$clamp(max = 70)) # Environmental inverse link function
  regP = (1 / (1 + torch_exp(-10 * (light - parReg))) - 1 / (1 + torch_exp(10 * parReg))) / (1 - 1 / (1 + torch_exp(10 * (1 - parReg))))
  reg_floor = if(is.null(self$reg_floor)) 0.2 else self$reg_floor
  mean = (regP*(environment[,NULL])$`repeat`(c(1, species$shape[2], 1)) + reg_floor)

  # Beverton-Holt cap: K = exp(reg_logK), reshaped to broadcast along the species
  # (last) dimension of `mean`. cap = "all" -> length 1 -> (1,1,1); cap =
  # "species" -> length N_species -> (1,1,N_species).
  K = torch::torch_exp(self$reg_logK)$reshape(c(1L, 1L, -1L))
  mean = K * mean / (K + mean)

  if(debug == TRUE) out = list(regP = regP, mean = mean) else out = mean
  return(out)
}

#' Regeneration limited by conspecific adult density
#'
#' A mechanistic extension of [FINN::regeneration_capped] in which
#' recruitment needs a seed source: the mean recruitment of a species scales
#' with the density of its own ADULTS on the site, the stems whose diameter
#' exceeds an adult threshold `x`.
#'
#' For species `k` on a site,
#' \deqn{A_k = \frac{1}{P a} \sum_{\mathrm{cohorts\ of\ } k} n \, \sigma\!\left(\frac{d - x_k}{\tau}\right)}{A_k = sum(n * sigmoid((d - x_k)/tau)) / (P * a)}
#' is the conspecific adult density (stems ha^-1), pooled over the `P` patches
#' of the site (area `a` each) because seeds disperse beyond one patch. The soft
#' threshold (`tau` = 2 dbh units) keeps `x_k` differentiable; a hard `d > x`
#' would give it no gradient. The seed-limitation term
#' \deqn{f(A_k) = \frac{A_k}{A_k + A_{1/2,k}}}{f(A_k) = A_k / (A_k + A_half_k)}
#' is 0 without adults and approaches 1 once adults are abundant, so the mean
#' recruitment
#' \deqn{m = r(\mathrm{light}) \, e^{\mathrm{env}} f(A) + \mathrm{reg\_floor}}{m = regP(light) * exp(env) * f(A) + reg_floor}
#' falls to the floor (read it as immigration from outside the site) where the
#' species has no adult, and rises to the light- and environment-driven
#' potential of [FINN::regeneration] where adults are many. The result is capped
#' with the same Beverton-Holt term as [FINN::regeneration_capped],
#' `K m / (K + m)`.
#'
#' The cap `reg_logK` is created by [FINN::finn()] as for
#' [FINN::regeneration_capped], with its resolution from `createProcess(cap = )`
#' and its prior from `prior`/`prior_lambda`, so it need not be declared. A
#' `reg_logK` declared in `custom_parameters` (as in the example) is kept and
#' sets the resolution by its length, but its start value is replaced by the
#' resolved prior at the first [FINN::fit()]. Declare the two adult parameters
#' with `createProcess(custom_parameters = )`; a length-1 init gives one shared
#' value and a length-`N_species` init one value per species:
#' * `reg_logK` (optional, see above): the log Beverton-Holt cap, stems ha^-1
#'   step^-1;
#' * `reg_adult_logx`: the log adult threshold, `x = exp(reg_adult_logx)` in the
#'   units of `dbh`. Once `x` falls below the smallest diameter present, every
#'   stem counts as an adult and `x` stops receiving gradient;
#' * `reg_adult_logA`: the log half-saturation adult density,
#'   `A_half = exp(reg_adult_logA)` stems ha^-1.
#'
#' `forward()` passes `dbh` and `trees` to a regeneration function only when its
#' signature has them, so the other regeneration functions are unaffected.
#' Called without `dbh`/`trees` (as the ALE code in `xAI.R` does), the function
#' returns recruitment with the seed source present (`f = 1`).
#'
#' @inheritParams regeneration_capped
#' @param dbh torch.Tensor cohort diameters, `[sites, patches, cohorts]`.
#' @param trees torch.Tensor stems per cohort, `[sites, patches, cohorts]`.
#'
#' @return torch.Tensor mean recruitment, stems ha^-1 step^-1,
#'   `[sites, patches, species]`. With `debug = TRUE`, a list that also holds
#'   `adults` (A, stems ha^-1, `[sites, 1, species]`) and `seed` (f(A)).
#'
#' @seealso [FINN::regeneration_capped], [FINN::regeneration_adult_effect()]
#'   (adult density as an effect, not a gate), [FINN::createProcess()],
#'   [recruitment] for the recruitment options and how to combine them
#'
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # three sites of the bundled Oregon FIA sample (?fia_obs_dt)
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#'
#' FINN.seed(1)
#' m <- finn(
#'   N_species = Nsp, recruits_dbh = 12.9,
#'   regeneration_process = createProcess(
#'     ~ temp + prec, FINN::regeneration_adult,
#'     custom_parameters = list(reg_logK       = rep(log(50), Nsp),
#'                              reg_adult_logx = rep(log(20), Nsp),   # adult from 20 cm
#'                              reg_adult_logA = rep(log(20), Nsp)),  # f = 0.5 at 20 adults/ha
#'     optimizeSpecies = TRUE, optimizeEnv = TRUE))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06, plot_progress = FALSE)
#' # the fitted cap, adult threshold (cm) and half-saturation adult density (per ha)
#' exp(as.numeric(m$reg_logK)); exp(as.numeric(m$reg_adult_logx)); exp(as.numeric(m$reg_adult_logA))
#'
#' @import torch
#' @export
regeneration_adult = function(species, parReg, pred, light, dbh = NULL, trees = NULL, debug = FALSE) {
  needed = c("reg_logK", "reg_adult_logx", "reg_adult_logA")
  missing_par = needed[vapply(needed, function(p) is.null(self[[p]]), logical(1))]
  if(length(missing_par) > 0)
    stop("`regeneration_adult` needs the parameter(s) ", paste(missing_par, collapse = ", "),
         ". Add them in createProcess(custom_parameters = list(reg_logK = ..., ",
         "reg_adult_logx = ..., reg_adult_logA = ...)) -- length 1 for a shared value, ",
         "length N_species for one per species.", call. = FALSE)

  if(self$record_raws) {
    self$raw_r = c(self$raw_r,  list(as_array(light$unsqueeze(4))))
  }

  if("matrix" %in% class(pred)) pred = torch::torch_tensor(pred)
  environment = torch::torch_exp(pred) # Environmental inverse link function
  regP = (1 / (1 + torch_exp(-10 * (light - parReg))) - 1 / (1 + torch_exp(10 * parReg))) / (1 - 1 / (1 + torch_exp(10 * (1 - parReg))))
  reg_floor = if(is.null(self$reg_floor)) 0.2 else self$reg_floor
  potential = regP*(environment[,NULL])$`repeat`(c(1, species$shape[2], 1))
  n_sp = potential$shape[3]

  if(is.null(dbh) || is.null(trees)) {
    adults = NULL
    seed = torch::torch_ones_like(potential)
  } else {
    tau = 2.0 # width of the soft adult threshold, in dbh units
    x = torch::torch_exp(self$reg_adult_logx)$expand(n_sp)  # a length-1 init is shared
    A_half = torch::torch_exp(self$reg_adult_logA)$reshape(c(1L, 1L, -1L))
    if(dbh$shape[3] == 0) {
      adults = torch::torch_zeros(potential$shape[1], 1L, n_sp, device = potential$device)
    } else {
      # adult weight of every cohort against its own species' threshold
      w = torch::torch_sigmoid((dbh - x[species]) / tau)
      zero = torch::torch_zeros(potential$shape[1], n_sp, device = potential$device)
      # sum over cohorts AND patches -> [sites, species], then per ha of the site
      adults = aggregate_results(species, list(trees * w), list(zero), sp_max = n_sp)[[1]]
      adults = (adults / (species$shape[2] * self$patch_size_ha))$unsqueeze(2)
    }
    seed = adults / (adults + A_half)
  }

  mean = potential * seed + reg_floor

  # Beverton-Holt cap, as in regeneration_saturation
  K = torch::torch_exp(self$reg_logK)$reshape(c(1L, 1L, -1L))
  mean = K * mean / (K + mean)

  if(debug == TRUE) out = list(regP = regP, adults = adults, seed = seed, mean = mean) else out = mean
  return(out)
}

#' Regeneration with an effect of conspecific adult density
#'
#' A variant of [FINN::regeneration_adult] in which the density of conspecific
#' adults is an EFFECT on recruitment rather than a gate: a species without
#' adults on the site still recruits a share `f0` of its potential, so bare
#' ground can be colonised from an empty stand.
#'
#' The conspecific adult density `A_k` (stems ha^-1) is computed as in
#' [FINN::regeneration_adult]: the stems of species `k` weighted by a soft adult
#' threshold `sigmoid((d - x_k) / tau)` (`tau` = 2 dbh units), pooled over the
#' patches of the site. The adult term is
#' \deqn{f(A_k) = f_{0,k} + (1 - f_{0,k}) \frac{A_k}{A_k + A_{1/2,k}}}{f(A_k) = f0_k + (1 - f0_k) * A_k / (A_k + A_half_k)}
#' with \eqn{f_0 = \mathrm{sigmoid}(\mathrm{reg\_adult\_logitf0}) \in (0, 1)}{f0 = sigmoid(reg_adult_logitf0) in (0, 1)},
#' and the mean recruitment is
#' \deqn{m = r(\mathrm{light}) \, e^{\mathrm{env}} f(A) + \mathrm{reg\_floor},}{m = regP(light) * exp(env) * f(A) + reg_floor,}
#' capped by the Beverton-Holt term of [FINN::regeneration_capped],
#' `K m / (K + m)`.
#'
#' `f0` is the share of the potential recruitment a species gets on a site
#' without conspecific adults -- seed rain from outside the site, advance
#' regeneration. It rises to 1 as adults become abundant. `f0 -> 1` means no
#' adult effect (the process becomes [FINN::regeneration_capped]); `f0 -> 0`
#' recovers the gate of [FINN::regeneration_adult], where a species without
#' adults recruits only at `reg_floor`.
#'
#' The soft threshold counts small stems partly: a 1-cm stem is about 1.5 % of
#' an adult at `x = 9.4` cm, so a low `A_half` gives a non-zero adult term even
#' in a stand of seedlings.
#'
#' The cap `reg_logK` is created by [FINN::finn()] as for
#' [FINN::regeneration_capped] (resolution from `createProcess(cap = )`, prior
#' from `prior`/`prior_lambda`), so it need not be declared. Declare the three
#' adult parameters with `createProcess(custom_parameters = )`; a length-1 init
#' gives one shared value and a length-`N_species` init one value per species:
#' * `reg_adult_logx`: the log adult threshold, `x = exp(reg_adult_logx)` in the
#'   units of `dbh`;
#' * `reg_adult_logA`: the log half-saturation adult density,
#'   `A_half = exp(reg_adult_logA)` stems ha^-1;
#' * `reg_adult_logitf0`: the logit of the share without adults,
#'   `f0 = plogis(reg_adult_logitf0)`.
#'
#' As for [FINN::regeneration_adult], `forward()` passes `dbh` and `trees`
#' because the signature has them. Called without `dbh`/`trees` (as the ALE
#' code in `xAI.R` does), the function returns recruitment with the seed source
#' present (`f = 1`).
#'
#' @inheritParams regeneration_adult
#'
#' @return torch.Tensor mean recruitment, stems ha^-1 step^-1,
#'   `[sites, patches, species]`. With `debug = TRUE`, a list that also holds
#'   `adults` (A, stems ha^-1, `[sites, 1, species]`), `f0` and `seed` (f(A)).
#'
#' @seealso [FINN::regeneration_adult], [FINN::regeneration_capped],
#'   [FINN::createProcess()], [recruitment] for the recruitment options and how
#'   to combine them
#'
#' @examplesIf torch::torch_is_installed()
#' library(data.table)
#' # three sites of the bundled Oregon FIA sample (?fia_obs_dt)
#' sites <- 1:3
#' Nsp   <- max(fia_obs_dt$species)
#' init  <- makeInitCohorts(fia_init_trees[siteID %in% sites], Nspecies = Nsp)
#'
#' FINN.seed(1)
#' m <- finn(
#'   N_species = Nsp, recruits_dbh = 12.9,
#'   regeneration_process = createProcess(
#'     ~ temp + prec, FINN::regeneration_adult_effect,
#'     custom_parameters = list(reg_adult_logx    = rep(log(20), Nsp),   # adult from 20 cm
#'                              reg_adult_logA    = rep(log(20), Nsp),   # half effect at 20 adults/ha
#'                              reg_adult_logitf0 = qlogis(0.2)),        # 20 % without adults, shared
#'     optimizeSpecies = TRUE, optimizeEnv = TRUE))
#' fit(m, data = fia_obs_dt[siteID %in% sites], env = fia_env_dt[siteID %in% sites],
#'     init_cohort = init, epochs = 2L, patches = 4, patch_size = 0.06, plot_progress = FALSE)
#' # the fitted share of the potential recruitment without conspecific adults
#' plogis(as.numeric(m$reg_adult_logitf0))
#'
#' @import torch
#' @export
regeneration_adult_effect = function(species, parReg, pred, light, dbh = NULL, trees = NULL, debug = FALSE) {
  needed = c("reg_logK", "reg_adult_logx", "reg_adult_logA", "reg_adult_logitf0")
  missing_par = needed[vapply(needed, function(p) is.null(self[[p]]), logical(1))]
  if(length(missing_par) > 0)
    stop("`regeneration_adult_effect` needs the parameter(s) ", paste(missing_par, collapse = ", "),
         ". Add them in createProcess(custom_parameters = list(reg_adult_logx = ..., ",
         "reg_adult_logA = ..., reg_adult_logitf0 = ...)) -- length 1 for a shared value, ",
         "length N_species for one per species. reg_logK is created by finn().", call. = FALSE)

  if(self$record_raws) {
    self$raw_r = c(self$raw_r,  list(as_array(light$unsqueeze(4))))
  }

  if("matrix" %in% class(pred)) pred = torch::torch_tensor(pred)
  environment = torch::torch_exp(pred) # Environmental inverse link function
  regP = (1 / (1 + torch_exp(-10 * (light - parReg))) - 1 / (1 + torch_exp(10 * parReg))) / (1 - 1 / (1 + torch_exp(10 * (1 - parReg))))
  reg_floor = if(is.null(self$reg_floor)) 0.2 else self$reg_floor
  potential = regP*(environment[,NULL])$`repeat`(c(1, species$shape[2], 1))
  n_sp = potential$shape[3]

  # share of the potential without conspecific adults; species is the last dimension
  f0 = torch::torch_sigmoid(self$reg_adult_logitf0)$reshape(c(1L, 1L, -1L))

  if(is.null(dbh) || is.null(trees)) {
    adults = NULL
    seed = torch::torch_ones_like(potential)
  } else {
    tau = 2.0 # width of the soft adult threshold, in dbh units
    x = torch::torch_exp(self$reg_adult_logx)$expand(n_sp)  # a length-1 init is shared
    A_half = torch::torch_exp(self$reg_adult_logA)$reshape(c(1L, 1L, -1L))
    if(dbh$shape[3] == 0) {
      adults = torch::torch_zeros(potential$shape[1], 1L, n_sp, device = potential$device)
    } else {
      # adult weight of every cohort against its own species' threshold
      w = torch::torch_sigmoid((dbh - x[species]) / tau)
      zero = torch::torch_zeros(potential$shape[1], n_sp, device = potential$device)
      # sum over cohorts AND patches -> [sites, species], then per ha of the site
      adults = aggregate_results(species, list(trees * w), list(zero), sp_max = n_sp)[[1]]
      adults = (adults / (species$shape[2] * self$patch_size_ha))$unsqueeze(2)
    }
    seed = f0 + (1 - f0) * adults / (adults + A_half)
  }

  mean = potential * seed + reg_floor

  # Beverton-Holt cap, as in regeneration_capped
  K = torch::torch_exp(self$reg_logK)$reshape(c(1L, 1L, -1L))
  mean = K * mean / (K + mean)

  if(debug == TRUE) out = list(regP = regP, adults = adults, f0 = f0, seed = seed, mean = mean) else out = mean
  return(out)
}




#' @rdname regeneration_capped
#' @export
regeneration_saturation = regeneration_capped
# Share the body rather than call through: a process function is rebound to the
# model by set_environment(), and a wrapper that calls regeneration_capped() by
# name resolves it in the namespace, where there is no `self`.
body(regeneration_saturation) = bquote({
  .Deprecated("regeneration_capped", old = "regeneration_saturation")
  .(body(regeneration_capped))
})
