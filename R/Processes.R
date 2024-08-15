#' Calculate the BA_stemsal area of a tree given the diameter at breast height (dbh).
#'
#' This function calculates the BA_stemsal area of a tree given the diameter at breast height (dbh).
#'
#' @param dbh torch.Tensor The diameter at breast height of the tree.
#'
#' @return torch.Tensor The BA_stemsal area of the tree.
#'
#' @examples
#' dbh = torch::torch_tensor(50)
#' BA_stemsal_area = BA_stem(dbh)
#' print(BA_stemsal_area)
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
#' <img src="figures/BA_stand_plot2.png" alt="dbh, trees, basal area" style="max-width:70%;"/>
#'
#' Sensitivity of basal area for different combinations of dbh, number of trees to patch size.
#'
#' <img src="figures/BA_stand_plot1.png" alt="Patch size, dbh, trees, basal area" style="max-width:70%;"/>
#'
#' @return A numeric value representing the basal area of the stand in square meters per hectare.
#' @examples
#' # Example usage
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
#' cohort_df1$basal_area <- torch::as_array(BA_stand(cohort$dbh, cohort$trees, patch_size_ha = patch_size))
#'
#' # View the first few rows of the resulting data frame
#' head(cohort_df1)
#'
#' # Basic plot showing the function of trees and dbh for basal area
#'
#' # only keep rows with basal area <100
#' cohort_df1 <- cohort_df1[cohort_df1$basal_area < 100,]
#'
#' library(ggplot2)
#' ggplot(cohort_df1, aes(x = dbh, y = basal_area, color = trees, group = trees)) +
#'   geom_line() +
#'   ylab("Basal Area (m^2/ha)") +
#'   xlab("Diameter at Breast Height (cm)") +
#'   scale_color_viridis_c(name = "Trees per ha", trans = "log10", option = "magma", direction = -1) +
#'   ggtitle("Basal Area as a Function of Trees and DBH")
#' @export
BA_stand <- function(dbh, trees, patch_size_ha) {
  return((pi * (dbh / 100 / 2)^2 * trees) / patch_size_ha)
}



#' Calculate the height of a tree based on its diameter at breast height and an alometry parameter.
#'
#' @param dbh A numeric value representing the diameter at breast height of the tree in cm.
#' @param parHeight A numeric value representing the species height alometry.
#'
#' @details
#'
#' This function calculates the height of a tree based on the diameter at breast height (dbh) and a parameter parHeight.
#'
#' The height is calculated using the formula:
#' \deqn{height = \left( \exp \left( \frac{(\text{dbh} \times \text{parHeight})}{(\text{dbh} + 100)} \right) - 1 \right) \times 100 + 0.001}
#' where dbh is the diameter at breast height of the tree in cm and parHeight is an alometric species specific parameter.
#'
#' All parameters of parHeight from 0 to 1 result in physiologicaly plausible heights.
#' The range from 0.3 to 0.9 results in realistic tree heights.
#' Values of parHeight close to 1 are physiologically almost impossible, below 0.3 is suitable for small tree species and shrubs.
#'
#' <br>
#' <img src="figures/height_plot1.png" alt="Parameter range" style="max-width:70%;"/>
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


#' Compute the fraction of available light (light) for each cohort BA_stemsed on the given parameters.
#'
#' This function calculates the fraction of available light for each cohort of trees BA_stemsed on their diameter at breast height (dbh), species, number of trees, and global parameters.
#'
#' @param dbh torch.Tensor Diameter at breast height for each cohort.
#' @param species torch.Tensor species index for each cohort.
#' @param trees Number of trees.
#' @param parHeight torch.Tensor global parameters for all species.
#' @param h torch.Tensor (Optional) Height of each cohort. Defaults to NULL.
#' @param minLight float (Optional) Minimum light requirement. Defaults to 50.
#'
#' @return torch.Tensor Fraction of available light (light) for each cohort.
#' @import torch
#' @examples
#' competition(dbh = torch::torch_tensor(c(10, 15, 20)), species = torch::torch_tensor(c(1, 2, 1)),
#'         trees = 100, parHeight = torch::torch_tensor(c(0.3, 0.5)), h = torch::torch_tensor(c(5, 7, 6)), minLight = 40)
#' @export
competition = function(dbh, species, trees, parHeight, h = NULL, minLight = 50., patch_size_ha, ba = NULL, cohortHeights = NULL){


  if(is.null(ba)) ba = BA_stand(dbh = dbh, trees = trees, patch_size_ha = patch_size_ha)
  if(is.null(cohortHeights)) cohortHeights = height(dbh, parHeight[species])$unsqueeze(4)
  if(is.null(h)) {
    h = cohortHeights

    # TODO:
    # Rework the height comparison, the sigmoid function preserves the gradients but they are not great (either -1 or 1)

    ba_height = (ba$unsqueeze_(4)$multiply(((cohortHeights - h$permute(c(1,2, 4, 3)) - 0.1)/1e-2)$sigmoid_() ))$sum(-2) # AUFPASSEN
  }else{
    ba_height = (ba$unsqueeze_(4)$multiply_(((cohortHeights - 0.1)/1e-2)$sigmoid_() ))$sum(-2)
  }
  light = 1.-ba_height/minLight
  light = torch_clamp(light, min = 0)
  return(light)
}



#' Mortality
#'
#' @param dbh dbh
#' @param species species
#' @param trees trees
#' @param parMort parMort
#' @param pred predictions
#' @param light available light
#'
#' @export
mortality = function(dbh, species, trees, parMort, pred, light, base_steepness = 10, debug = F) {
  # TODO remove constant part
  # shade = 1-torch_sigmoid((light + (1-parMort[,1][species]) - 1)/(1/10^(1.5 + torch_abs(light-0.5))))

  # Scale steepness towards the edges
  scaled_steepness <- base_steepness / (0.5 - abs(parMort[,1][species] - 0.5))
  shade = 1 - ((1 / (1 + torch::torch_exp(-scaled_steepness * (light - parMort[,1][species]))) - 1 / (1 + torch::torch_exp(scaled_steepness * parMort[,1][species]))) /
         (1 / (1 + torch::torch_exp(-scaled_steepness * (1 - parMort[,1][species]))) - 1 / (1 + torch::torch_exp(scaled_steepness * parMort[,1][species]))))

  environment = index_species(pred, species)
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
  mort1 = binomial_from_gamma(trees+trees$le(0.5)$float(), predM)*trees$ge(0.5)$float()
  mort2 = mort1 + mort1$round()$detach() - mort1$detach()
  if(debug == TRUE) out = list(shade = shade, light = light, environment = environment, gPSize = gPSize, predM = predM, mort1 = mort1, mort2 = mort2) else out = mort2
  return(out)
}



#' Calculate growth
#'
#' This function calculates growth BA_stemsed on specified parameters.
#'
#' @param dbh torch.Tensor Diameter at breast height.
#' @param species torch.Tensor species of tree.
#' @param parGrowth torch.Tensor Growth parameters.
#' @param parMort torch.Tensor Mortality parameters.
#' @param pred torch.Tensor Predicted values.
#' @param light torch.Tensor Accumulated Light.
#'
#' @return torch.Tensor A tensor representing the forest plot growth.
#'
#' @import torch
#'
#' @export
growth = function(dbh, species, parGrowth, pred, light, light_steepness = 10, debug = F){
  # K = parGrowth[,3][species]
  # K = 0
  # light_steepness = 10

  # shade = torch_sigmoid((light + (1-parGrowth[,1][species]) - 1)/1e-1)
  shade = ((1 / (1 + torch::torch_exp(-light_steepness * (light - parGrowth[,1][species]))) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]))) /
         (1 / (1 + torch::torch_exp(-light_steepness * (1 - parGrowth[,1][species]))) - 1 / (1 + torch::torch_exp(light_steepness * parGrowth[,1][species]))))
  environment = index_species(pred, species)
  # growth = (1.- torch.pow(1.- pred,4.0)) * parGrowth[species,1]
  growth = shade * environment * (torch::torch_exp(-(dbh / (parGrowth[,2][species] * 100))))
  # growth = environment * shade * torch::torch_exp(-0.5 * (log(dbh / (parGrowth[,2][species])*100) / K)^2)
  # growth = parGrowth[species,1]
  # return torch.nn.functional.softplus(growth)
  if(debug == TRUE) out = list(shade = shade, light = light, environment = environment,growth = growth) else out = growth
  return(out)
}

#' Calculate the regeneration of forest patches BA_stemsed on the input parameters.
#'
#' This function calculates the regeneration of forest patches BA_stemsed on species information, regeneration parameters, prediction values, and available light.
#'
#' @param species torch.Tensor species information.
#' @param parReg torch.Tensor Regeneration parameters. 0 <= parReg <= 1
#' This parameter denotes the fraction of light needed for a species to regenerate.
#' In general low values for high regeneration and high values for low regeneration.
#' @param pred torch.Tensor Prediction values.
#' @param light torch.Tensor Available light variable for calculation.
#'
#' @return torch.Tensor Regeneration values for forest patches.
#'
#' @import torch
#' @importFrom torch torch_sigmoid
#' @export
regeneration = function(species, parReg, pred, light, patch_size_ha, debug = F) {
  if("matrix" %in% class(pred)) pred = torch::torch_tensor(pred)
  environment = pred
  regP = (1 / (1 + torch_exp(-10 * (light - parReg))) - 1 / (1 + torch_exp(10 * parReg))) / (1 - 1 / (1 + torch_exp(10 * (1 - parReg))))
  #regP = torch_sigmoid((light + (1-parReg) - 1)/1e-3) # TODO masking? better https://pytorch.org/docs/stable/generated/torch.masked_select.html
  mean = (regP*(environment[,NULL])$`repeat`(c(1, species$shape[2], 1))+0.2)
  regeneration1 = sample_poisson_gaussian(mean*patch_size_ha) # TODO, check if exp or not?! lambda should be always positive!
  regeneration2 = regeneration1 + regeneration1$round()$detach() - regeneration1$detach()
  if(debug == T) out = list(regP = regP, mean = mean, regeneration1 = regeneration1, regeneration2 = regeneration2) else out = regeneration2
  return(out)
}



