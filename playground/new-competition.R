library(FINN)
library(ggplot2)
library(torch)
library(data.table)

patch_size = 0.1
Nsp = 8
trees_vec = c(0:100)
dbh_vec = dbh = seq(0,100, 1)

cohort_df1 =
  data.frame(
    patchID = 1,
    cohortID = 1,
    expand.grid(
      trees = round(trees_vec*patch_size),
      dbh = dbh_vec,
      species = 1:Nsp
    )
  )
cohort_df1$siteID = 1:nrow(cohort_df1)


cohort = CohortMat$new(obs_df = cohort_df1)
dbh = cohort$dbh
species = cohort$species
trees = cohort$trees

cohort_df1$ba_stand = BA_stand(dbh = cohort_df1$dbh, trees = cohort_df1$trees, patch_size_ha = patch_size)

cohort_df1 = cohort_df1[cohort_df1$ba_stand <= 150,]

cohort = CohortMat$new(obs_df = cohort_df1)
# competition = function(dbh, species, trees, parHeight, h = NULL, minLight = 50., patch_size_ha, ba = NULL, cohortHeights = NULL){
#   if(is.null(ba)) ba = BA_stand(dbh = dbh, trees = trees, patch_size_ha = patch_size_ha)
#   if(is.null(cohortHeights)) cohortHeights = height(dbh, parHeight[species])$unsqueeze(4)
#   if(is.null(h)) {
#     h = cohortHeights
#     ba_height = (ba$unsqueeze_(4)$multiply(((cohortHeights - h$permute(c(1,2, 4, 3)) - 0.1)/1e-2)$sigmoid_() ))$sum(-2) # AUFPASSEN
#   }else{
#     ba_height = (ba$unsqueeze_(4)$multiply_(((cohortHeights - 0.1)/1e-2)$sigmoid_() ))$sum(-2)
#   }
#   light = 1.-ba_height/minLight
#   light = torch_clamp(light, min = 0)
#   return(light)
# }


plot_dt_out <- data.table()
for(i in seq(0,2,0.1)){

  # growth parameters
  parComp = matrix(c(
    runif(Nsp, 0.3, 0.7), # parHeight
    runif(Nsp, i, i) # Competition strength
  ),Nsp, 2)
  parComp = torch_tensor(cbind(parComp[,1],parComp[,2]), requires_grad=TRUE, dtype=torch_float32())

  competition2 = function(dbh, species, trees, parComp, h = NULL, patch_size_ha, ba = NULL, cohortHeights = NULL){
    parHeight = parComp[,1]
    parCompStr = parComp[,2]
    if(is.null(ba)) ba = BA_stand(dbh = dbh, trees = trees, patch_size_ha = patch_size_ha)*parCompStr[species]*0.1
    if(is.null(cohortHeights)) cohortHeights = height(dbh, parHeight[species])$unsqueeze(4)
    if(is.null(h)) {
      h = cohortHeights
      ba_height = (ba$unsqueeze_(4)$multiply(((cohortHeights - h$permute(c(1,2, 4, 3)) - 0.1)/1e-2)$sigmoid_() ))$sum(-2) # AUFPASSEN
    }else{
      ba_height = (ba$unsqueeze_(4)$multiply_(((cohortHeights - 0.1)/1e-2)$sigmoid_() ))$sum(-2)
    }
    light = 1-ba_height
    light = torch_clamp(light, min = 0)
    return(light)
  }



  comp = competition(
    cohort$dbh, cohort$species, cohort$trees,
    parComp = parComp, h=0,
    patch_size_ha = 0.1)


  parHeight = parComp[,1]
  parCompStr = parComp[,2]

  comp2 = competition2(
    cohort$dbh, cohort$species, cohort$trees,
    parComp = parComp, h=0,
    patch_size_ha = 0.1)

  cohort_df1$light = torch::as_array(comp)[,1,1]
  cohort_df1$light2 = torch::as_array(comp2)[,1,1]

  # p1 = ggplot(cohort_df1,aes(y = cut(ba_stand, seq(0,100,10)), x = factor(trees/0.1), fill = light))+
  #   geom_tile()+
  #   ylab("trees/ha [N]")+
  #   xlab("mean stand dbh [cm]")+
  #   guides(fill=guide_legend(title="light [%]\nat forest floor\n(h=0)"))+
  #   ggtitle(i)
  #   # ggtitle("competition",paste0(deparse(competition),collapse = "\n"))

  cohort_dt = data.table(cohort_df1)

  cohort_dt[light2 > 1]
  plot_dt <- cohort_dt[, .(
    light = unique(light),
    light2 = unique(light2),
    ba_stand = unique(ba_stand),
    parCompStr = i
    ), by = .(siteID)]

  plot_dt_out <- rbind(plot_dt_out, plot_dt)

  # p2 = ggplot(cohort_df1,aes(y = factor(trees/0.1), x = dbh, fill = light2))+
  #   geom_tile()+
  #   ylab("trees/ha [N]")+
  #   xlab("mean stand dbh [cm]")+
  #   guides(fill=guide_legend(title="light [%]\nat forest floor\n(h=0)"))+
  #   ggtitle(i)
  #   # ggtitle("competition",paste0(deparse(competition),collapse = "\n"))
  #
  # print(cowplot::plot_grid(p1,p2, ncol = 2))
}

ggplot(plot_dt_out)+
  # geom_line(aes(x = ba_stand, y = light, color = "minlight = 50", linetype = "dashed"), linewidth = 3)+
  geom_line(aes(x = ba_stand, y = light2, color = factor(parCompStr), linetype = "solid"))









