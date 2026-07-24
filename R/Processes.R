#' Calculate the basal area of a tree given the diameter at breast height (dbh)
#'
#' This function calculates the basal area of a tree given the diameter at breast height (dbh).
#'
#' @param dbh torch.Tensor The diameter at breast height of the tree.
#'
#' @return torch.Tensor The basal area of the tree.
#'
#' @examples
#' \dontrun{
#' dbh = torch::torch_tensor(50)
#' basal_area = BA_stem(dbh)
#' print(basal_area)
#' }
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
#' @return A numeric value representing the basal area of the stand in square meters per hectare.
#' @examplesIf torch::torch_is_installed()
#' # BA_stand operates on torch tensors, so this runs only where the torch
#' # backend (libtorch) is available.
#' dbh_vec <- seq(1, 200, 1)
#' trees_vec <- c(0:500, 10^(seq(2, 4, length.out = 20)))
#'
#' # Generate test data for a patch size of 0.1
#' patch_size <- 0.1
#' cohort_df1 <- expand.grid(
#'   trees_ha = trees_vec,
#'   patch_size_ha = patch_size,
#'   dbh = dbh_vec
#' )
#'
#' cohort_df1 <- data.frame(
#'   patchID = 1,
#'   cohortID = 1,
#'   species = 1,
#'   cohort_df1
#' )
#'
#' cohort_df1$siteID <- 1:nrow(cohort_df1)
#' cohort_df1$trees <- round(cohort_df1$trees_ha * patch_size)
#'
#' cohort <- CohortMat$new(obs_df = cohort_df1)
#'
#' cohort_df1$basal_area <- torch::as_array(
#'   BA_stand(cohort$dbh, cohort$trees, patch_size_ha = patch_size))
#'
#' # View the first few rows of the resulting data frame
#' head(cohort_df1)
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
#' All parameters of parHeight from 0 to 1 result in physiologicaly plausible heights.
#' The range from 0.3 to 0.9 results in realistic tree heights.
#' Values of parHeight close to 1 are physiologically almost impossible, below 0.3 is suitable for small tree species and shrubs.
#'
#' \figure{height_plot1.png}{Parameter range}
#'
#' @return A numeric value representing the calculated height of the tree.
#' @examples
#' height(30, 0.5)
#' height(c(30), c(0.5,0.3))
#' height(c(30,20), c(0.5))
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
#' @param dbh torch.Tensor Diameter at breast height for each cohort.
#' @param species torch.Tensor species index for each cohort.
#' @param trees torch.Tensor Number of trees in each cohort.
#' @param parComp torch.Tensor Competition / height-allometry parameters per species.
#' @param h torch.Tensor (Optional) Height of each cohort. Defaults to NULL.
#' @param patch_size_ha numeric Patch size in hectares.
#' @param ba torch.Tensor (Optional) Pre-computed basal area. Defaults to NULL.
#' @param cohortHeights torch.Tensor (Optional) Pre-computed cohort heights. Defaults to NULL.
#' @param n_quantiles integer Number of height quantiles used when `continuous = FALSE`. Defaults to 10.
#' @param continuous logical Use the continuous competition formulation. Defaults to FALSE.
#'
#' @return torch.Tensor Fraction of available light (light) for each cohort.
#' @import torch
#' @examples
#' \dontrun{
#' competition(dbh = torch::torch_tensor(c(10, 15, 20)), species = torch::torch_tensor(c(1, 2, 1)),
#'         trees = 100, parComp = torch::torch_tensor(c(0.3, 0.5)), patch_size_ha = 0.1)
#' }
#' @export
competition = function(dbh, species, trees, parComp, h = NULL, patch_size_ha, ba = NULL, cohortHeights = NULL, n_quantiles = 10, continuous = FALSE){
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
      # Quantiles and their indices
      cohortHeights_quant = torch::torch_quantile(cohortHeights, q = quants, dim = 3)$squeeze(4)$permute(c(2, 3, 1))
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

# TODO:
# Rework the height comparison, the sigmoid function preserves the gradients but they are not great (either -1 or 1)


#' mortality_wo_growth (internal)
#' @noRd
mortality_wo_growth = function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = F) {
  # TODO remove constant part
  # shade = 1-torch_sigmoid((light + (1-parMort[,1][species]) - 1)/(1/10^(1.5 + torch_abs(light-0.5))))

  # Scale steepness towards the edges
  scaled_steepness <- base_steepness / (0.5 - abs(parMort[,1][species] - 0.5))
  shade = 1 - ((1 / (1 + torch::torch_exp(-scaled_steepness * (light - parMort[,1][species]))) - 1 / (1 + torch::torch_exp(scaled_steepness * parMort[,1][species]))) /
         (1 / (1 + torch::torch_exp(-scaled_steepness * (1 - parMort[,1][species]))) - 1 / (1 + torch::torch_exp(scaled_steepness * parMort[,1][species]))))

  environment = torch::torch_sigmoid(pred)
  # gPSize = torch_clamp(0.1*(dbh/torch_clamp((parMort[,2][species]*100), min = 0.00001))$pow(2.3), max = 1.0)
  gPSize = (1-torch::torch_exp(-(dbh / (parMort[,2][species] * 100))))
  # gPSize = torch_sigmoid(gPSize)
  # TODO
  # clamp can lead to vanishing gradients, sigmoid is not perfect but probably better here!
  # shade -> [0,1]
  # gPsize -> [0, 1]
  # environment -> [0, 1]
  # -> raw pred -> [0, 5]
  predM = environment*((shade*gPSize + shade + gPSize)/3)
  # predM = torch_sigmoid((environment*(shade+gPSize) + shade*gPSize + shade + gPSize -1.5)*2 )
  if(debug == TRUE) out = list(shade = shade, light = light, environment = environment, gPSize = gPSize, predM = predM)
  else out = predM
  return(out)
}

#' Mortality
#'
#' @param dbh dbh
#' @param species species
#' @param trees trees
#' @param parMort parMort
#' @param pred predictions
#' @param light available light
#' @param base_steepness numeric Steepness of the shade-response sigmoid. Defaults to 5.
#' @param debug logical If TRUE, return the intermediate components as a list. Defaults to FALSE.
#' @param growth torch.Tensor (Optional) Growth entering the mortality response; defaults to the model's current growth.
#'
#' @return A `torch` tensor of per-cohort mortality probabilities; or, if
#'   `debug = TRUE`, a list of the intermediate components.
#' @export
mortality = function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = F, growth = NULL) {
  if(is.null(growth)) growth = self$g
  if(self$record_raws) {
    self$raw_m = c(self$raw_m,  list(as_array( torch::torch_cat(list(dbh$unsqueeze(4),
                                                                     growth$unsqueeze(4),
                                                                     light$unsqueeze(4),
                                                                     trees$unsqueeze(4),
                                                                     species$unsqueeze(4)$float()), dim = 4)) ))
  }

  environment = torch::torch_sigmoid(pred+growth*parMort[,3][species]+light*parMort[,1][species] +  parMort[,2][species]*(dbh / ( 100)))
  predM = environment
  if(debug == TRUE) out = list(shade = shade, light = light, environment = environment, gPSize = gPSize, predM = predM)
  else out = predM
  return(out)
}



mortality_hybrid = function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = F, growth = NULL) {
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



#' Calculate growth
#'
#' This function calculates growth based on specified parameters.
#'
#' @param dbh torch.Tensor Diameter at breast height.
#' @param species torch.Tensor species of tree.
#' @param parGrowth torch.Tensor Growth parameters.
#' @param pred torch.Tensor Predicted values.
#' @param light torch.Tensor Accumulated Light.
#' @param light_steepness numeric Steepness of the light-response sigmoid. Defaults to 10.
#' @param debug logical If TRUE, return the intermediate components as a list. Defaults to FALSE.
#' @param trees torch.Tensor (Optional) Number of trees per cohort. Defaults to NULL.
#'
#' @return torch.Tensor A tensor representing the forest plot growth.
#'
#' @import torch
#'
#' @export
growth = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F, trees = NULL){

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


growth_hybrid= function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F, trees = NULL) {

  if(self$record_raws) {
    self$raw_g = c(self$raw_g,  list(as_array( torch::torch_cat(list(dbh$unsqueeze(4),
                                                                     light$unsqueeze(4),
                                                                     trees$unsqueeze(4),
                                                                     species$unsqueeze(4)$float()), dim = 4)) ))
  }

  g = (self$nn_growth(dbh = dbh, trees = trees, light = light, species = species, env = pred)- exp(1))$exp()
  return(g)
}


#' Calculate the regeneration of forest patches based on the input parameters
#'
#' This function calculates the regeneration of forest patches based on species information, regeneration parameters, prediction values, and available light.
#'
#' @param species torch.Tensor species information.
#' @param parReg torch.Tensor Regeneration parameters. 0 <= parReg <= 1
#' This parameter denotes the fraction of light needed for a species to regenerate.
#' In general low values for high regeneration and high values for low regeneration.
#' @param pred torch.Tensor Prediction values.
#' @param light torch.Tensor Available light variable for calculation.
#' @param debug logical If TRUE, return the intermediate components as a list. Defaults to FALSE.
#'
#' @return torch.Tensor Regeneration values for forest patches.
#'
#' @import torch
#' @importFrom torch torch_sigmoid
#' @export
regeneration = function(species, parReg, pred, light, debug = F) {
  if(self$record_raws) {
    self$raw_r = c(self$raw_r,  list(as_array(light$unsqueeze(4))))
  }


  if("matrix" %in% class(pred)) pred = torch::torch_tensor(pred)
  environment = torch::torch_exp(pred) # Environmental inverse link function
  regP = (1 / (1 + torch_exp(-10 * (light - parReg))) - 1 / (1 + torch_exp(10 * parReg))) / (1 - 1 / (1 + torch_exp(10 * (1 - parReg))))
  mean = (regP*(environment[,NULL])$`repeat`(c(1, species$shape[2], 1))+0.2)
  #regP = torch_sigmoid((light + (1-parReg) - 1)/1e-3) # TODO masking? better https://pytorch.org/docs/stable/generated/torch.masked_select.html
  if(debug == T) out = list(regP = regP, mean = mean) else out = mean
  return(out)
}

regeneration_hybrid = function(species, parReg, pred, light, debug = F) {

  if(self$record_raws) {
    self$raw_r = c(self$raw_r,  list(as_array(light$unsqueeze(4))))
  }

  # clamp the pre-exp logit to avoid float overflow (-> Inf -> NaN gradients)
  # if the net's raw output drifts up during training; exp(10) ~ 22000/ha is
  # already far beyond any plausible recruitment rate
  r = self$nn_regeneration(light = light, species = species, env = pred)$clamp(max = 10)$exp()
}



