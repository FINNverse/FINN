library(testthat)
library(FINN)

compF_P = function(dbh, Species, nTree, parHeight, h = NULL, minLight = 50.){

  ba = (BA_P(dbh)*nTree)/0.1
  cohortHeights = height_P(dbh, parHeight[Species])$unsqueeze(4)
  if(is.null(h)) {
    h = cohortHeights
    BA_height = (ba$unsqueeze(4)*torch_sigmoid((cohortHeights - h$permute(c(1,2, 4, 3)) - 0.1)/1e-3) )$sum(-2) # AUFPASSEN
  }else{
    BA_height = (ba$unsqueeze(4)*torch_sigmoid((cohortHeights - 0.1)/1e-3))$sum(-2)
  }
  AL = 1.-BA_height/minLight
  AL = torch_clamp(AL, min = 0)
  return(AL)
}

cohort = CohortMat$new(
  obs_df = data.frame(
    siteID = c(1,2),
    patchID = c(1,1),
    cohortID = c(1,1),
    species = c(1,1),
    nTree = c(5094,120),
    dbh = c(6.2,60.5)
  ))

dbh = cohort$dbh
Species = cohort$Species
nTree = cohort$nTree

ba = (BA_P(dbh)*nTree)/1
torch::as_array(ba)
