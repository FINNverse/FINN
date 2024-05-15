library(testthat)
library(FINN)
library(ggplot2)

patch_size = 0.1

trees_vec = c(1:99,10^(seq(2,4, length.out = 16)))
dbh_vec = dbh = seq(1,300, 1)

cohort_df1 =
  data.frame(
    patchID = 1,
    cohortID = 1,
    species = 1,
    expand.grid(
      trees = round(trees_vec*patch_size),
      dbh = dbh_vec
    )
  )
cohort_df1$siteID = 1:nrow(cohort_df1)

cohort = CohortMat$new(obs_df = cohort_df1)

dbh = cohort$dbh
species = cohort$species
trees = cohort$trees


comp = competition(cohort$dbh, cohort$species, cohort$trees,
                   parHeight = torch::torch_tensor(0.5), h=0)


cohort_df1$light = torch::as_array(comp)[,1,1]

ggplot(cohort_df1,aes(y = factor(trees/0.1), x = dbh, fill = light))+
  geom_tile()+
  ylab("trees/ha [N]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="light [%]\nat forest floor\n(h=0)"))+
  ggtitle("competition",paste0(deparse(competition),collapse = "\n"))

# parHeight = torch::torch_tensor(0.5)
# minLight = 50
# cohortHeights = height(dbh, parHeight[species])$unsqueeze(4)
# torch::as_array(cohortHeights)
#
# ba = BA_stand(cohort$dbh, cohort$trees)
# torch::as_array(ba)
# ba_height = (ba$unsqueeze(4)*torch::torch_sigmoid((cohortHeights - 0.1)/1e-2))$sum(-2)
# torch::as_array(ba_height)
# (light = 1.-ba_height/minLight)
# (light = torch::torch_clamp(light, min = 0))
# torch::as_array(light)
#
# h = cohortHeights
# ba_height = (ba$unsqueeze(4)*torch::torch_sigmoid((cohortHeights - h$permute(c(1,2, 4, 3)) - 0.1)/1e-2) )$sum(-2) # AUFPASSEN
# torch::as_array(ba_height)
# (light = 1.-ba_height/minLight)
# (light = torch::torch_clamp(light, min = 0))
# torch::as_array(light)

comp = competition(cohort$dbh, cohort$species, cohort$trees,
                   parHeight = torch::torch_tensor(0.5), h=0)


cohort_df1$light = torch::as_array(comp)[,1,1]

ggplot(cohort_df1,aes(y = factor(trees/0.1), x = dbh, fill = light))+
  geom_tile()+
  ylab("trees/ha [N]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="light [%]\nat forest floor\n(h=0)"))+
  ggtitle("competition",paste0(deparse(competition),collapse = "\n"))
