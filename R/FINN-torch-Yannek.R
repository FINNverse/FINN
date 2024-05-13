BA_P = function(dbh) {
  # Calculate the basal area of a tree given the diameter at breast height (dbh).
  #
  # Args:
  # dbh (torch.Tensor): The diameter at breast height of the tree.
  #
  # Returns:
  # torch.Tensor: The basal area of the tree.
  #
  # Example:
  # dbh = torch.tensor(50)
  # basal_area = BA_P(dbh)
  # print(basal_area)
  return(pi*(dbh/100./2.)^2.0)
}

# YANNEK
height_P = function(dbh, parGlobal) {
  # Calculate the height of a tree based on its diameter at breast height (dbh) and a parameter (par).
  # Args:
  # dbh (torch.Tensor): The diameter at breast height of the tree.
  # par (torch.Tensor): The parameter used to calculate the height.
  #
  # Returns:
  # torch.Tensor: The calculated height of the tree.
  height = (exp((((dbh * parGlobal) / (dbh+100))))-1)*100 + 0.001
  return(height)
}

# YANNEK
@torch.jit.script
compF_P = function(dbh, Species, nTree, parGlobal, h = NULL, minLigh = 50.){

  # Compute the fraction of available light (AL) for each cohort based on the given parameters.
  #
  # Args:
  # dbh (torch.Tensor): Diameter at breast height for each cohort.
  # Species (torch.Tensor): Species index for each cohort.
  # parGlobal (torch.Tensor): Global parameters for all species.
  # h (Optional[torch.Tensor], optional): Height of each cohort. Defaults to None.
  # minLight (float, optional): Minimum light requirement. Defaults to 50.
  #
  # Returns:
  # torch.Tensor: Fraction of available light (AL) for each cohort.
  ba = (BA_P(dbh)*nTree)/0.1
  cohortHeights = height_P(dbh, parGlobal[Species])$unsqueeze(3)
  if(is.null(h)) {
    h = cohortHeights
    BA_height = (ba$unsqueeze(3)*torch_sigmoid((cohortHeights - h$permute(0,1, 3, 2) - 0.1)/1e-3) )$sum(-2) # AUFPASSEN
  }else{
    BA_height = (ba.unsqueeze(3)*torch.sigmoid((cohortHeights - 0.1)/1e-3)).sum(-2)
  }
  AL = 1.-BA_height/minLight
  AL = torch_clamp(AL, min = 0)
  return(AL)
}

# YANNEK
growthFP = function(dbh, Species, parGlobal, parGrowth, parMort, pred, AL){

  # Calculate the growth of a forest stand based on the given parameters.
  #
  # Args:
  # dbh (torch.Tensor): Diameter at breast height.
  # Species (torch.Tensor): Species of the trees.
  # parGlobal (torch.Tensor): Global parameters.
  # parGrowth: Growth parameters.
  # pred (torch.Tensor): Predicted values.
  #
  # Returns:
  # torch.Tensor: Growth of the forest stand.
  # shade = (((AL**2)*parGrowth[Species,0]).sigmoid()-0.5)*2

  shade = torch_sigmoid((AL + (1-parGrowth[Species,1]) - 1)/1e-1)
  environment = index_species(pred, Species)
  pred = (shade*environment)
  # growth = (1.- torch.pow(1.- pred,4.0)) * parGrowth[Species,1]
  growth = pred/2 * parGrowth[Species,2] * ((parMort[Species,2]-dbh/100) / parMort[Species,2])^(2)
  # growth = parGrowth[Species,1]
  # return torch.nn.functional.softplus(growth)
  return(torch_clamp(growth, min = 0.0))
}

# Yannek
regFP = function(dbh, Species, parGlobal, parReg, pred, AL) {

  # Calculate the regeneration of forest patches based on the input parameters.
  #
  # Args:
  # dbh (torch.Tensor): Diameter at breast height.
  # Species (torch.Tensor): Species information.
  # parGlobal (torch.Tensor): Global parameters.
  # parReg (torch.Tensor): Regression parameters.
  # pred (torch.Tensor): Prediction values.
  #
  # Returns:
  # torch.Tensor: Regeneration values for forest patches.
  #
  test_tens$
  torch_repeat_interleave()
  regP = torch_sigmoid((AL + (1-parReg) - 1)/1e-3)
  environment = pred
  regeneration = sample_poisson_relaxed((regP*environment[,NULL]$repeat(c(1, Species$shape[2], 1)) + 1e-20), 20)
  regeneration = regeneration + regeneration$round()$detach() - regeneration$detach()
  return(regeneration)
}
