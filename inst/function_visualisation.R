# Load necessary library
library(ggplot2)
library(data.table)
library(torch)
library(FINN)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## create pdf ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# Open the existing PDF file in append mode
# pdf("function_vis.pdf")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## height ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
test_dt =
  data.table(
    expand.grid(
      list(
        dbh = seq(0, 500,1),
        parHeight = seq(0,1,0.1)
      )
    )
  )
test_dt$rowID = 1:nrow(test_dt)

test_dt[, height := height(dbh, parHeight), by = rowID]

p_height <-
  ggplot(test_dt, aes(x = dbh, y = height, group = parHeight,color = parHeight))+
    ggtitle("height",paste0(deparse(height),collapse = "\n"))+
    ylab("height(dbh,parHeight)")+
    xlab("dbh [cm]")+
    geom_line()+
    geom_label(aes(label = parHeight, color = parHeight), test_dt[,.(height = max(height),dbh = max(dbh)),by=parHeight])
print(p_height)

#### parHeight ####
#' in principle all parameters of parHeight from 0 to 1 result in physiologicaly
#' plausible heights.
#' But very high values close to 1 are really on the edge of how tall a tree
#' can become.
#'
#' Recommendation for default parameterization.
#' The default parameterization should cover the range from 0.3 to 0.9
#'

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## BA_stem ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# Generate test data
test_dt <- data.table(
  dbh = seq(0, 500, 1)
)

# Calculate basal area
test_dt[, basal_area := BA_stem(dbh)]

# Plot
p_BA_stem <- ggplot(test_dt, aes(x = dbh, y = basal_area)) +
  ggtitle("BA_stem", paste0(deparse(BA_stem), collapse = "\n")) +
  ylab("BA_stem(dbh)") +
  xlab("dbh [cm]") +
  geom_line()
print(p_BA_stem)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## BA_stand ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
trees_vec = c(0:500,10^(seq(2,4, length.out = 20)))
dbh_vec = dbh = seq(1,200, 1)
patch_size_vec = seq(0.04, 0.2, 0.02)
patch_size_vec = c(0.01,0.05,0.08,seq(0.1,1,0.2),1)

test_dt_out <- data.table()
i=patch_size_vec[1]
for(i in patch_size_vec){
  # Generate test data
  cohort_df1 =
    data.frame(
      patchID = 1,
      cohortID = 1,
      species = 1,
      expand.grid(
        trees_ha = trees_vec,
        patch_size_ha = i,
        dbh = dbh_vec
      )
    )
  cohort_df1$siteID = 1:nrow(cohort_df1)
  cohort_df1$trees = round(cohort_df1$trees_ha*i)
  cohort = CohortMat$new(obs_df = cohort_df1)

  basal_area = BA_stand(cohort$dbh, cohort$trees, patch_size_ha = i)
  light = competition(dbh = cohort$dbh, species = cohort$species,
                      trees = cohort$trees, parHeight = torch::torch_tensor(c(0.5)),
                      patch_size_ha = i,
                      h = 0)
  cohort_df1$basal_area = torch::as_array(basal_area)[,1,1]
  cohort_df1$light = torch::as_array(light)[,1,1]
  cohort_df1$patch_size_ha = i
  test_dt_out <- rbind(test_dt_out,cohort_df1)
}

max_ba = 60
test_dt_out[basal_area > max_ba, basal_area := max_ba,]
p_BA_stand <- ggplot(test_dt_out[trees_ha <=500],aes(y = factor(trees_ha), x = dbh, fill = basal_area))+
  geom_tile()+
  ylab("trees [N/ha]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="basal area\n[m^2/ha]"))+
  ggtitle("BA_stand",paste0(deparse(BA_stand),collapse = "\n"))+
  scale_fill_viridis_c(breaks = seq(0,max_ba,10), labels = c(seq(0,max_ba-10,10),paste0(">",max_ba)))+
  scale_y_discrete(breaks = seq(0,500,100))+
  facet_wrap(~patch_size_ha, labeller = labeller(patch_size_ha = function(x) paste0("patch size = ",x, " ha")))
print(p_BA_stand)


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## competition ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# sensitivity of competition to parHeight (species)
parHeight_vec = c(0.3,0.7)
for(i in parHeight_vec){
  test_dt_out <- data.frame()
  i=patch_size_vec[1]
  # Generate test data
  cohort_df1 =
    data.frame(
      patchID = 1,
      cohortID = 1,
      species = 1,
      expand.grid(
        trees_ha = trees_vec,
        patch_size_ha = 0,
        dbh = dbh_vec
      )
    )
  cohort_df1$siteID = 1:nrow(cohort_df1)
  cohort_df1$trees = round(cohort_df1$trees_ha*i)
  cohort = CohortMat$new(obs_df = cohort_df1)

  }

# sensitivity of competition to h: light availability vs. height

trees_vec = c(0:500,10^(seq(2,4, length.out = 20)))
dbh_vec = dbh = seq(1,200, 1)
patch_size_vec = seq(0.04, 0.2, 0.02)
patch_size_vec = c(0.1)



test_dt_out <- data.table()
# Generate test data
cohort_df1 =
  data.frame(
    patchID = 1,
    cohortID = 1,
    species = 1,
    expand.grid(
      trees_ha = trees_vec,
      patch_size_ha = 0.1,
      dbh = dbh_vec
    )
  )
cohort_df1$siteID = 1:nrow(cohort_df1)
cohort_df1$trees = round(cohort_df1$trees_ha*i)
cohort = CohortMat$new(obs_df = cohort_df1)

basal_area = BA_stand(cohort$dbh, cohort$trees, patch_size_ha = 0.1)
light = competition(dbh = cohort$dbh, species = cohort$species,
                    trees = cohort$trees, parHeight = torch::torch_tensor(parHeight_vec),
                    h = 1, patch_size_ha = 0.1)
cohort_df1$basal_area = torch::as_array(basal_area)[,1,1]
cohort_df1$light = torch::as_array(light)[,1,1]
test_dt_out <- rbind(test_dt_out,cohort_df1)

p_competition1 <- ggplot(test_dt_out[trees_ha <=500],aes(y = factor(trees_ha), x = dbh, fill = light))+
  geom_tile()+
  ylab("trees [N/ha]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="light[%]"))+
  ggtitle("competition",paste0(deparse(competition),collapse = "\n"))+
  # scale_fill_viridis_c(breaks = seq(0,max_ba,10), labels = c(seq(0,max_ba-10,10),paste0(">",max_ba)))+
  scale_y_discrete(breaks = seq(0,500,100))+
  facet_wrap(~patch_size_ha, labeller = labeller(patch_size_ha = function(x) paste0("patch size = ",x, " ha")))
print(p_competition1)


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

light = competition(cohort$dbh, cohort$species, cohort$trees,
                   parHeight = torch::torch_tensor(0.5), h=0, patch_size_ha = 0.1)
cohort_df1$light = torch::as_array(light)[,1,1]

ggplot(cohort_df1,aes(y = factor(trees/0.1), x = dbh, fill = light*100))+
  geom_tile()+
  ylab("trees/ha [N]")+
  xlab("mean stand dbh [cm]")+
  guides(fill=guide_legend(title="light [%]\nat forest floor\n(h=0)"))+
  ggtitle("competition",paste0(deparse(competition),collapse = "\n"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## regeneration ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

patch_size = 0.1
# trees_vec = c(1:10,10^(seq(2,4, length.out = 8)))
trees_vec = c(0:100)
dbh_vec = 80
# dbh_vec = dbh = seq(30,100, 10)

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

light = competition(cohort$dbh, cohort$species, cohort$trees,
                   parHeight = torch::torch_tensor(0.5), h=0, patch_size_ha = 0.1)
cohort_df1$light = torch::as_array(light)[,1,1]
hist(cohort_df1$light)

sim_dt <-
  data.table(
    expand.grid(list(
      parReg = seq(0,1,0.02),
      pred = seq(0,500,5)
    ))
  )

# r = regeneration(species = cohort$species, parReg = sim_dt[i,]$parReg, pred = array(sim_dt[i,]$pred,dim = c(1,1)), light = comp)
# regeneration2 = function(species, parReg, pred, light) {
#   if("matrix" %in% class(pred)) pred = torch::torch_tensor(pred)
#   environment = pred
#   regP = torch_sigmoid((light + (1-parReg) - 1)/1e-3) # TODO masking? better https://pytorch.org/docs/stable/generated/torch.masked_select.html
#   mean = (regP*(environment[,NULL])$`repeat`(c(1, species$shape[2], 1))+0.2 )
#   regeneration1 = sample_poisson_relaxed(mean, num_samples = 1000) # TODO, check if exp or not?! lambda should be always positive!
#   regeneration2 = regeneration1 + regeneration1$round()$detach() - regeneration1$detach()
#   return(list(regP = regP, mean = mean, regeneration1 = regeneration1, regeneration2 = regeneration2))
# }

light1 = 0.5
parReg = seq(0,1,0.001)
plot(parReg,plogis((light1 + (1-parReg) - 1)/1e-2))

i=1
out_dt <- data.table()
for (i in 1:nrow(sim_dt)) {
  # r = regeneration(species = cohort$species, parReg = sim_dt[i,]$parReg, pred = array(sim_dt[i,]$pred,dim = c(1,1)), light = comp)
  r = regeneration(species = cohort$species, parReg = sim_dt[i,]$parReg, pred = array(sim_dt[i,]$pred,dim = c(1,1)), light = light, patch_size_ha = 0.1, debug = T)
  cohort_df1$regP = as_array(r$regP)
  cohort_df1$mean = as_array(r$mean)
  cohort_df1$r1 = as_array(r$regeneration1)
  cohort_df1$r2 = as_array(r$regeneration2)
  cohort_df1$parReg = sim_dt[i,]$parReg
  cohort_df1$pred = sim_dt[i,]$pred
  out_dt <- rbind(out_dt, cohort_df1)
}

ggplot(out_dt, aes(x = pred, y = parReg, fill = mean))+
  geom_tile()+
  facet_wrap(~paste0("light=", cut(light,6)))+
  theme(legend.position = "top")+
  guides(fill = guide_colorbar(barwidth = 10, title.position = "top"))+
  labs(x = "predicted environment effect", y = "parReg (species light requirement for regeneration)", fill = "lambda (mean of poisson)")

ggplot(out_dt, aes(x = pred, y = parReg, fill = r2))+
  geom_tile()+
  facet_wrap(~paste0("light=", cut(light,6)))+
  theme(legend.position = "top")+
  guides(fill = guide_colorbar(barwidth = 10, title.position = "top"))+
  labs(x = "predicted environment effect", y = "parReg (species light requirement for regeneration)", fill = "simulated regeneration per patch")

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## mortality ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

patch_size = 0.1
# trees_vec = c(1:10,10^(seq(2,4, length.out = 8)))
trees_vec = c(50)
dbh_vec = seq(10,300,10)
# dbh_vec = dbh = seq(30,100, 10)

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
#
# cohort_df1b =
#   data.frame(
#     patchID = 1,
#     cohortID = 2,
#     species = 1,
#     expand.grid(
#       trees = round(trees_vec*patch_size),
#       dbh = dbh_vec+1
#     )
#   )
# cohort_df1b$siteID = 1:nrow(cohort_df1b)
#
# cohort_df1 <- data.table(rbind(cohort_df1, cohort_df1b))

cohort = CohortMat$new(obs_df = cohort_df1)

dbh = cohort$dbh
species = cohort$species
trees = cohort$trees

# light1 = competition(cohort$dbh, cohort$species, cohort$trees,
#                     parHeight = torch::torch_tensor(0.5), patch_size_ha = 0.1)
# dim(torch::as_array(light1))
# cohort_df1[["light1"]] = c(torch::as_array(light1)[,1,1],torch::as_array(light1)[,1,2])
# cohort_df1$light2 = torch::as_array(light1)[,1,2]
# hist(cohort_df1$light1)
# hist(cohort_df1$light2)
# hist(cohort_df1[cohortID == 1]$light1)
# hist(cohort_df1$light2)
# plot(cohort_df1$light1, cohort_df1$light2, xlim = range(as_array(light1)), ylim = range(as_array(light1)))
#
# light2 = competition(cohort$dbh, cohort$species, cohort$trees,
#                     parHeight = torch::torch_tensor(0.5), patch_size_ha = 0.1)
# dim(torch::as_array(light1))
# dim(torch::as_array(light2))
# cohort_df1$light2 = torch::as_array(light2)[,1,2]
# hist(cohort_df1$light2)
# hist(cohort_df1$light)

sim_dt <-
  data.table(
    expand.grid(list(
      parMort1 = seq(0,1,0.25),
      parMort2 = seq(0,4,1),
      pred = seq(0,1,0.1),
      light = seq(0,1,0.1)
    ))
  )


i=50
out_dt <- data.table()
for (i in 1:nrow(sim_dt)) {
  if(i==1) cat("\n")
  cat("\r", i, "of", nrow(sim_dt))
  temp_parMort = torch_tensor(cbind(sim_dt[i,]$parMort1,sim_dt[i,]$parMort2), requires_grad=TRUE, dtype=torch_float32())
  temp_pred = torch_tensor(array(sim_dt[i,]$pred,dim = c(length(as_array(cohort$dbh)),1)), requires_grad=TRUE, dtype=torch_float32())
  light = torch_tensor(array(sim_dt[i,]$light,dim = c(length(as_array(cohort$dbh)),1,1)), requires_grad=TRUE, dtype=torch_float32())
  # r = regeneration(species = cohort$species, parReg = sim_dt[i,]$parReg, pred = array(sim_dt[i,]$pred,dim = c(1,1)), light = comp)
  m = mortality(dbh = cohort$dbh, species = cohort$species, trees = cohort$trees, parMort = temp_parMort, pred = temp_pred, light = light, debug = T)
  cohort_df1$shade = as_array(m$shade)
  cohort_df1$dbh = as_array(cohort$dbh)
  cohort_df1$environment = as_array(m$environment)
  cohort_df1$gPSize = as_array(m$gPSize)
  cohort_df1$predM = as_array(m$predM)
  cohort_df1$mort1 = as_array(m$mort1)
  cohort_df1$mort2 = as_array(m$mort2)
  cohort_df1$pred = sim_dt[i,]$pred
  cohort_df1$parMort1 = sim_dt[i,]$parMort1
  cohort_df1$parMort2 = sim_dt[i,]$parMort2
  cohort_df1$light = sim_dt[i,]$light
  out_dt <- rbind(out_dt, cohort_df1)
  if(i==nrow(sim_dt)) cat("\n")
}

ggplot(out_dt[dbh %in% c(seq(30,300,30))], aes(x = parMort2, y = gPSize, color = factor(dbh), group = dbh))+
  geom_line()+
  coord_cartesian(ylim = c(0,2))+
  geom_hline(yintercept = 0.1, col = "black", lty = 3)+
  scale_y_continuous(breaks = seq(0,2,0.1))

ggplot(out_dt, aes(x = factor(parMort1), y = shade))+
  geom_boxplot()+
  facet_wrap(~cohortID)

ggplot(out_dt[,
              .(shade = mean(shade)
              ), by = .(light = cut(light, seq(0, 1, 0.2), ordered_result = T, include.lowest = T), cohortID, parMort1)
              ], aes(
  x = parMort1,
  y = shade,
  color = light)) +
  geom_point() +
  geom_line() +
  facet_grid(light~ cohortID)

ggplot(out_dt, aes(
  x = parMort1,
  y = shade,
  color = light)) +
  geom_point() +
  stat_smooth(se = F) +
  facet_grid(light~pred)

ggplot(out_dt, aes(x = factor(parMort1), y = predM))+
  geom_boxplot()

ggplot(out_dt, aes(x = pred, y = predM))+
  geom_point()

custom_palette <- c("#440154", "#3B528B", "#21918C", "#5DC863", "#FDE725")

ggplot(out_dt[pred == 0.5 & parMort2 %in% c(1,2,3,4) & parMort1 %in% c(0,0.25,0.5,0.75, 1)], aes(x = dbh, y = light, fill = predM))+
  geom_tile()+
  facet_grid(paste0("parMort1=", parMort1)~paste0("parMort2=", parMort2))+
  theme(legend.position = "top")+
  guides(fill = guide_colorbar(barwidth = 10, title.position = "top"))+
  scale_fill_gradientn(colors = custom_palette)
  # scale_fill_viridis(option = "viridis", name = "Value")
  # labs(x = "predicted environment effect", y = "parReg (species light requirement for regeneration)", fill = "mort")
